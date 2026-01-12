#' Unblinded sample size re-estimation for recurrent events
#'
#' Estimates the event rates and dispersion from unblinded interim data
#' and calculates the required sample size to maintain power, assuming the
#' planned treatment effect holds (or using the observed control rate).
#'
#' @param data A data frame containing the unblinded interim data. Must include
#'   columns `events` (number of events), `tte` (total exposure/follow-up time),
#'   and `treatment` (treatment group identifier, e.g., 1 for control, 2 for experimental).
#'   This is typically the output of [cut_data_by_date()].
#' @param ratio Planned allocation ratio (experimental / control). Default is 1.
#' @param lambda1_planning Planned event rate for the control group used in original calculation.
#' @param lambda2_planning Planned event rate for the experimental group used in original calculation.
#' @param rr0 Rate ratio under the null hypothesis (lambda2/lambda1). Default is 1.
#' @param power Target power (1 - beta). Default is 0.8.
#' @param alpha One-sided significance level. Default is 0.025.
#' @param accrual_rate Vector of accrual rates (patients per unit time).
#' @param accrual_duration Vector of durations for each accrual rate. Must be same length
#'   as `accrual_rate`.
#' @param trial_duration Total planned duration of the trial.
#' @param dropout_rate Dropout rate (hazard rate). Default is 0.
#' @param max_followup Maximum follow-up time for any patient. Default is NULL (infinite).
#' @param event_gap Gap duration after each event during which no new events are counted.
#'   Default is NULL (no gap).
#'
#' @return A list containing:
#'   \describe{
#'     \item{n_total_unblinded}{Re-estimated total sample size using unblinded estimates.}
#'     \item{dispersion_unblinded}{Estimated dispersion parameter (k) from unblinded data.}
#'     \item{lambda1_unblinded}{Estimated control event rate from unblinded data.}
#'     \item{lambda2_unblinded}{Estimated experimental event rate from unblinded data.}
#'     \item{info_fraction}{Estimated information fraction at interim (unblinded information / target information).}
#'     \item{unblinded_info}{Estimated statistical information from the unblinded interim data.}
#'     \item{target_info}{Target statistical information required for the planned power.}
#'   }
#'
#' @importFrom MASS glm.nb
#' @importFrom stats qnorm fitted coef vcov
#' @importFrom utils tail
#'
#' @export
unblinded_ssr <- function(data,
                          ratio = 1,
                          lambda1_planning,
                          lambda2_planning,
                          rr0 = 1,
                          power = 0.8,
                          alpha = 0.025,
                          accrual_rate,
                          accrual_duration,
                          trial_duration,
                          dropout_rate = 0,
                          max_followup = NULL,
                          event_gap = NULL) {

  # 1. Estimate parameters from unblinded data
  # We fit a negative binomial model to estimate rates and dispersion
  # We assume a common dispersion parameter across groups as per standard sample size formula
  
  # Ensure treatment is a factor
  data$treatment <- as.factor(data$treatment)
  
  # Fit NB model: events ~ treatment + offset(log(tte))
  # We use the observed rates for re-estimation of nuisance parameters,
  # but typically we might stick to the PLANNED effect size (hazard ratio) 
  # applied to the OBSERVED control rate, or use both observed rates if we trust them.
  # The prompt says "control rate is < the assumed control rate", implying we should use
  # the observed control rate. Standard practice for "maintaining power" often uses
  # observed control rate and PLANNED relative effect to avoid powering for a noisy observed effect.
  # However, let's calculate both observed rates first.
  
  fit <- tryCatch({
    MASS::glm.nb(events ~ treatment + offset(log(tte)), data = data)
  }, error = function(e) {
    stop("Failed to fit negative binomial model to interim data: ", e$message)
  })
  
  # Extract estimates
  k_est <- 1 / fit$theta
  
  # Get rates (exp(coef))
  # Intercept is log(rate_control) (assuming treatment 1 is reference)
  # Slope is log(rate_ratio)
  coefs <- coef(fit)
  lambda1_est <- exp(coefs[1])
  lambda2_est <- exp(coefs[1] + coefs[2])
  
  # For sample size re-estimation, we usually fix the effect size to the planning assumption
  # to avoid overpowering/underpowering due to early random high/low effect estimates.
  # But we use the observed control rate and dispersion.
  # Target effect size:
  target_rr <- lambda2_planning / lambda1_planning
  
  # Re-estimated parameters for calculation
  lambda1_calc <- lambda1_est
  lambda2_calc <- lambda1_est * target_rr # Maintain planned relative effect
  k_calc <- k_est
  
  # 2. Calculate Information
  # We need the information accumulated SO FAR and the TARGET information.
  
  # Information formula for NB (Wald statistic for log rate ratio):
  # I = 1 / (Var(beta_hat))
  # Per sample size formula: V_tilde = (1/mu1 + k)/p1 + (1/mu2 + k)/p2
  # Info = n_total / V_tilde
  
  # Calculate current information from the fitted model
  # The variance of the log rate ratio (coefficient for treatment) is the element (2,2) of vcov
  # if using treatment contrast.
  # vcov(fit) gives covariance of (Intercept, treatment).
  # Var(log(lambda2) - log(lambda1)) = Var(beta_treatment)
  var_log_rr <- vcov(fit)[2, 2]
  current_info <- 1 / var_log_rr
  
  # 3. Calculate Target Information
  # Target info depends on the PLANNED effect size and power.
  # I_target = (z_alpha + z_beta)^2 / (log(lambda1 * rr0 / lambda2))^2
  # Note: This uses the PLANNED rates for the effect size definition (log(rr)).
  # Actually, the denominator is the effect size we want to detect.
  # If we want to detect the PLANNED effect size, we use that.
  
  z_alpha <- stats::qnorm(1 - alpha) # One-sided
  z_beta <- stats::qnorm(power)
  log_rr_target <- log(lambda1_planning * rr0 / lambda2_planning) # Effect size to detect
  # Note: sample_size_nbinom uses log(lambda1 * rr0 / lambda2) in denominator?
  # Let's check sample_size_nbinom.R formula.
  # den <- (log(lambda1 * rr0 / lambda2))^2
  # Yes.
  
  target_info <- (z_alpha + z_beta)^2 / (log_rr_target)^2
  
  # 4. Re-estimate Sample Size
  # We use the sample_size_nbinom function with the NEW nuisance parameters (lambda1_calc, k_calc)
  # and the ORIGINAL target effect size (implied by lambda2_calc).
  
  # We need to solve for n_total such that we get target_info.
  # sample_size_nbinom returns n_total.
  
  res <- sample_size_nbinom(
    lambda1 = lambda1_calc,
    lambda2 = lambda2_calc,
    rr0 = rr0,
    dispersion = k_calc,
    power = power,
    alpha = alpha,
    sided = 1, # Assuming one-sided as per default
    ratio = ratio,
    accrual_rate = accrual_rate,
    accrual_duration = accrual_duration,
    trial_duration = trial_duration,
    dropout_rate = dropout_rate,
    max_followup = max_followup,
    event_gap = event_gap
  )
  
  list(
    n_total_unblinded = res$n_total,
    dispersion_unblinded = k_calc,
    lambda1_unblinded = lambda1_est,
    lambda2_unblinded = lambda2_est, # The observed one, for reference
    info_fraction = current_info / target_info,
    unblinded_info = current_info,
    target_info = target_info
  )
}
