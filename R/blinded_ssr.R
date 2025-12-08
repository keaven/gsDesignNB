#' Blinded Sample Size Re-estimation for Recurrent Events
#'
#' Estimates the blinded dispersion and event rate from aggregated interim data
#' and calculates the required sample size to maintain power, assuming the
#' planned treatment effect holds. This function supports constant rates (Friede &
#' Schmidli 2010) and accommodates future extensions for time-varying rates
#' (Schneider et al. 2013) by using the exposure-adjusted rate.
#'
#' @param data A data frame containing the blinded interim data. Must include
#'   columns `events` (number of events) and `tte` (total exposure/follow-up time).
#'   This is typically the output of [cut_data_by_date()].
#' @param ratio Planned allocation ratio (experimental / control). Default is 1.
#' @param lambda1_planning Planned event rate for the control group used in original calculation.
#' @param lambda2_planning Planned event rate for the experimental group used in original calculation.
#' @param power Target power (1 - beta). Default is 0.8.
#' @param alpha One-sided significance level. Default is 0.025.
#' @param method Method for sample size recalculation. Currently "friede" (Friede & Schmidli 2010)
#'   is implemented, which uses the blinded nuisance parameter estimates.
#' @param accrual_rate Vector of accrual rates (patients per unit time).
#' @param accrual_duration Vector of durations for each accrual rate. Must be same length
#'   as `accrual_rate`.
#'
#' @references
#' Friede, T., & Schmidli, H. (2010). Blinded sample size reestimation with count data: methods and applications.
#' \emph{Statistics in Medicine}, 29(10), 1145-1156. \doi{10.1002/sim.3891}
#'
#' Schneider, S., Schmidli, H., & Friede, T. (2013). Blinded sample size reestimation for recurrent event data with time trends.
#' \emph{Statistics in Medicine}, 32(30), 5448-5457. \doi{10.1002/sim.5920}
#'
#' @param trial_duration Total planned duration of the trial.
#' @param dropout_rate Dropout rate (hazard rate). Default is 0.
#' @param max_followup Maximum follow-up time for any patient. Default is NULL (infinite).
#' @param event_gap Gap duration after each event during which no new events are counted.
#'   Default is NULL (no gap).
#'
#' @return A list containing:
#'   \describe{
#'     \item{n_total_unadjusted}{Original planned total sample size (based on planning parameters).}
#'     \item{n_total_blinded}{Re-estimated total sample size using blinded estimates.}
#'     \item{dispersion_blinded}{Estimated dispersion parameter (k) from blinded data.}
#'     \item{lambda_blinded}{Estimated overall event rate from blinded data.}
#'     \item{info_fraction}{Estimated information fraction at interim (blinded information / target information).}
#'     \item{blinded_info}{Estimated statistical information from the blinded interim data.}
#'     \item{target_info}{Target statistical information required for the planned power.}
#'   }
#' @export
#' @importFrom MASS glm.nb
#' @importFrom stats qnorm fitted
#' @importFrom utils tail
blinded_ssr <- function(data, ratio = 1, lambda1_planning, lambda2_planning,
                        power = 0.8, alpha = 0.025, method = "friede",
                        accrual_rate, accrual_duration, trial_duration,
                        dropout_rate = 0, max_followup = NULL, event_gap = NULL) {
  # Use calculate_blinded_info for parameter estimation and info calculation
  blind_info_res <- calculate_blinded_info(data, ratio, lambda1_planning, lambda2_planning)

  dispersion_est <- blind_info_res$dispersion_blinded
  lambda_est <- blind_info_res$lambda_blinded
  lambda1_new <- blind_info_res$lambda1_adjusted
  lambda2_new <- blind_info_res$lambda2_adjusted
  observed_info <- blind_info_res$blinded_info

  rate_ratio_planning <- lambda2_planning / lambda1_planning

  # Recalculate Sample Size (Friede method)
  # Calculate NEW sample size with observed k and updated lambdas
  res_new <- sample_size_nbinom(
    lambda1 = lambda1_new,
    lambda2 = lambda2_new,
    dispersion = dispersion_est,
    power = power,
    alpha = alpha,
    ratio = ratio,
    accrual_rate = accrual_rate,
    accrual_duration = accrual_duration,
    trial_duration = trial_duration,
    dropout_rate = dropout_rate,
    max_followup = max_followup,
    event_gap = event_gap
  )

  # Targeted Final Information
  z_alpha <- qnorm(1 - alpha)
  z_beta <- qnorm(power)
  target_info <- (z_alpha + z_beta)^2 / (log(rate_ratio_planning))^2

  info_fraction <- observed_info / target_info

  list(
    n_total_blinded = res_new$n_total,
    dispersion_blinded = dispersion_est,
    lambda_blinded = lambda_est,
    lambda1_adjusted = lambda1_new,
    lambda2_adjusted = lambda2_new,
    info_fraction = info_fraction,
    blinded_info = observed_info,
    target_info = target_info
  )
}
