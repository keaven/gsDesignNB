#' Calculate Blinded Statistical Information
#'
#' Estimates the blinded dispersion and event rate from aggregated interim data
#' and calculates the observed statistical information for the log rate ratio,
#' assuming the planned allocation ratio and treatment effect.
#'
#' @param data A data frame containing the blinded interim data. Must include
#'   columns `events` (number of events) and `tte` (total exposure/follow-up time).
#' @param ratio Planned allocation ratio (experimental / control). Default is 1.
#' @param lambda1_planning Planned event rate for the control group.
#' @param lambda2_planning Planned event rate for the experimental group.
#'
#' @return A list containing:
#'   \describe{
#'     \item{blinded_info}{Estimated statistical information.}
#'     \item{dispersion_blinded}{Estimated dispersion parameter (k).}
#'     \item{lambda_blinded}{Estimated overall event rate.}
#'     \item{lambda1_adjusted}{Re-estimated control rate.}
#'     \item{lambda2_adjusted}{Re-estimated experimental rate.}
#'   }
#' @export
#' @importFrom MASS glm.nb
#' @importFrom stats coef
calculate_blinded_info <- function(data, ratio = 1, lambda1_planning, lambda2_planning) {
  df <- as.data.frame(data)
  if (!all(c("events", "tte") %in% names(df))) {
    stop("Data must contain 'events' and 'tte' columns.")
  }

  # 1. Blinded Parameter Estimation
  # Fit Negative Binomial model to pooled data (intercept only)
  fit_blind <- tryCatch(
    suppressWarnings(MASS::glm.nb(events ~ 1 + offset(log(tte)), data = df)),
    error = function(e) NULL
  )

  if (is.null(fit_blind) || is.na(fit_blind$theta)) {
    warning("Negative Binomial fit failed on blinded data. Falling back to Poisson (dispersion = 0).")
    dispersion_est <- 0
    lambda_est <- sum(df$events) / sum(df$tte)
  } else {
    dispersion_est <- 1 / fit_blind$theta
    lambda_est <- exp(coef(fit_blind)[1])
  }

  # 2. Blinded Information Calculation
  p1 <- 1 / (1 + ratio)
  p2 <- ratio / (1 + ratio)

  # Rename hr_planning to rate_ratio_planning
  rate_ratio_planning <- lambda2_planning / lambda1_planning

  lambda1_new <- lambda_est / (p1 + p2 * rate_ratio_planning)
  lambda2_new <- lambda1_new * rate_ratio_planning

  n_current <- nrow(df)
  avg_exposure <- sum(df$tte) / n_current

  mu1_new <- lambda1_new * avg_exposure
  mu2_new <- lambda2_new * avg_exposure

  V_bar_blind <- (1 / mu1_new + dispersion_est) / p1 + (1 / mu2_new + dispersion_est) / p2

  observed_info <- n_current / V_bar_blind

  list(
    blinded_info = observed_info,
    dispersion_blinded = dispersion_est,
    lambda_blinded = lambda_est,
    lambda1_adjusted = lambda1_new,
    lambda2_adjusted = lambda2_new
  )
}
