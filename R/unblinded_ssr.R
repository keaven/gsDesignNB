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
#'     \item{fallback}{Character label for which estimator was used
#'       (\code{"ml"} or \code{"mom"}).}
#'   }
#'
#' @details
#' If the maximum likelihood negative binomial fit fails to converge, the
#' function falls back to method-of-moments estimation via [estimate_nb_mom()]
#' rather than erroring out. The observed Fisher information is then computed
#' analytically from the MoM-estimated rates and dispersion using the same
#' subject-level weight formula as [calculate_blinded_info()]. This keeps SSR
#' updates well-defined under extreme overdispersion or sparse interim data.
#'
#' @importFrom MASS glm.nb
#' @importFrom stats qnorm fitted coef vcov
#' @importFrom utils tail
#'
#' @export
#'
#' @examples
#' interim <- data.frame(
#'   events = c(1, 2, 1, 3),
#'   tte = c(0.8, 1.0, 1.2, 0.9),
#'   treatment = c("Control", "Control", "Experimental", "Experimental")
#' )
#' unblinded_ssr(
#'   interim,
#'   ratio = 1,
#'   lambda1_planning = 0.5,
#'   lambda2_planning = 0.3,
#'   power = 0.8,
#'   alpha = 0.025,
#'   accrual_rate = 10,
#'   accrual_duration = 12,
#'   trial_duration = 18
#' )
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
  
  # Ensure treatment is a factor with exactly two levels.
  data$treatment <- droplevels(as.factor(data$treatment))
  if (nlevels(data$treatment) != 2) {
    stop("unblinded_ssr requires interim data with exactly two treatment groups.",
         call. = FALSE)
  }

  # Fit NB model: events ~ treatment + offset(log(tte))
  # We use the observed rates for re-estimation of nuisance parameters,
  # but typically we might stick to the PLANNED effect size (hazard ratio) 
  # applied to the OBSERVED control rate, or use both observed rates if we trust them.
  # The prompt says "control rate is < the assumed control rate", implying we should use
  # the observed control rate. Standard practice for "maintaining power" often uses
  # observed control rate and PLANNED relative effect to avoid powering for a noisy observed effect.
  # However, let's calculate both observed rates first.
  
  fit <- tryCatch(
    suppressWarnings(MASS::glm.nb(events ~ treatment + offset(log(tte)), data = data)),
    error = function(e) NULL
  )
  ml_ok <- !is.null(fit) && isTRUE(fit$converged) && !is.na(fit$theta) &&
    is.finite(fit$theta) && fit$theta > 0

  if (ml_ok) {
    k_est <- 1 / fit$theta
    coefs <- coef(fit)
    lambda1_est <- exp(coefs[1])
    lambda2_est <- exp(coefs[1] + coefs[2])
    fallback <- "ml"
  } else {
    # Fall back to method-of-moments so that genuinely overdispersed data do
    # not collapse to a Poisson/ML failure. Matches the behaviour of
    # mutze_test() and calculate_blinded_info().
    mom <- tryCatch(
      estimate_nb_mom(data, group = "treatment"),
      error = function(e) NULL
    )
    if (is.null(mom) || any(is.na(mom$lambda)) || !is.finite(mom$dispersion)) {
      stop("Failed to fit negative binomial model to interim data, and ",
           "method-of-moments fallback also failed.", call. = FALSE)
    }
    warning(
      "Unblinded NB ML fit did not converge; falling back to method-of-moments ",
      "for nuisance parameter estimation.",
      call. = FALSE
    )
    lev <- levels(data$treatment)
    lambda1_est <- as.numeric(mom$lambda[lev[1]])
    lambda2_est <- as.numeric(mom$lambda[lev[2]])
    k_est <- as.numeric(mom$dispersion)
    fallback <- "mom"
  }
  
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
  
  # Calculate current information. Under ML we read it directly from vcov().
  # Under the MoM fallback we compute the observed Fisher information using
  # the MoM-estimated rates and dispersion and the same formula used
  # throughout the package: I = 1 / (1/W1 + 1/W2),
  # W_g = sum_i mu_{g,i} / (1 + k * mu_{g,i}).
  if (identical(fallback, "ml")) {
    var_log_rr <- vcov(fit)[2, 2]
    current_info <- 1 / var_log_rr
  } else {
    lev <- levels(data$treatment)
    t1 <- data$tte[data$treatment == lev[1]]
    t2 <- data$tte[data$treatment == lev[2]]
    mu1 <- lambda1_est * t1
    mu2 <- lambda2_est * t2
    w1 <- sum(mu1 / (1 + k_est * mu1))
    w2 <- sum(mu2 / (1 + k_est * mu2))
    current_info <- if (w1 > 0 && w2 > 0) 1 / (1 / w1 + 1 / w2) else 0
  }
  
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
    target_info = target_info,
    fallback = fallback
  )
}
