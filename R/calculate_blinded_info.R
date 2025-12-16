#' Calculate blinded statistical information
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
#' @param event_gap Optional. Gap duration (numeric) to adjust planning rates if provided.
#'   If provided, planning rates are adjusted as lambda / (1 + lambda * gap).
#'
#' @return A list containing:
#'   \describe{
#'     \item{blinded_info}{Estimated statistical information.}
#'     \item{dispersion_blinded}{Estimated dispersion parameter (k).}
#'     \item{lambda_blinded}{Estimated overall event rate.}
#'     \item{lambda1_adjusted}{Re-estimated control rate.}
#'     \item{lambda2_adjusted}{Re-estimated experimental rate.}
#'   }
#'
#' @importFrom MASS glm.nb
#' @importFrom stats coef
#'
#' @export
#'
#' @examples
#' interim <- data.frame(events = c(1, 2, 1, 3), tte = c(0.8, 1.0, 1.2, 0.9))
#' calculate_blinded_info(
#'   interim,
#'   ratio = 1,
#'   lambda1_planning = 0.5,
#'   lambda2_planning = 0.3
#' )
calculate_blinded_info <- function(data, ratio = 1, lambda1_planning, lambda2_planning, event_gap = NULL) {
  df <- as.data.frame(data)
  if (!all(c("events", "tte") %in% names(df))) {
    stop("Data must contain 'events' and 'tte' columns.")
  }

  # Drop non-positive exposure rows (cannot contribute to information)
  df <- df[df$tte > 0, , drop = FALSE]
  if (nrow(df) == 0) {
    return(list(
      blinded_info = 0,
      dispersion_blinded = NA_real_,
      lambda_blinded = NA_real_,
      lambda1_adjusted = NA_real_,
      lambda2_adjusted = NA_real_
    ))
  }

  .blinded_info_from_tte <- function(tte, lambda1, lambda2, dispersion, ratio) {
    p1 <- 1 / (1 + ratio)
    p2 <- ratio / (1 + ratio)

    mu1 <- lambda1 * tte
    mu2 <- lambda2 * tte

    # Fisher information for the log-rate in NB2 with log link and offset(log(tte)):
    # sum(mu / (1 + dispersion * mu))
    w1 <- sum(mu1 / (1 + dispersion * mu1))
    w2 <- sum(mu2 / (1 + dispersion * mu2))
    if (!is.finite(w1) || !is.finite(w2) || w1 <= 0 || w2 <= 0) {
      return(0)
    }

    var_log_rr <- 1 / (p1 * w1) + 1 / (p2 * w2)
    if (!is.finite(var_log_rr) || var_log_rr <= 0) {
      return(0)
    }
    1 / var_log_rr
  }

  # 1. Blinded Parameter Estimation
  # Fit negative binomial model to pooled data (intercept only)
  fit_blind <- tryCatch(
    suppressWarnings(MASS::glm.nb(events ~ 1 + offset(log(tte)), data = df)),
    error = function(e) NULL
  )

  if (is.null(fit_blind) || is.na(fit_blind$theta)) {
    warning("Negative binomial fit failed on blinded data. Falling back to Poisson (dispersion = 0).")
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

  # If event_gap is present, it reduces the effective rates used in variance calculation
  # But lambda_est is already the effective rate (events / exposure), where exposure accounts for gaps if cut_data_by_date used it.
  # However, the planning parameters lambda1_planning and lambda2_planning are likely "calendar" rates (without gap adjustment)
  # UNLESS the user passed effective rates.

  # Assuming lambda1_planning are raw rates:
  if (!is.null(event_gap) && event_gap > 0) {
    # Adjust planning ratio?
    # The ratio lambda2/lambda1 is roughly preserved even with gaps if rates are small,
    # but strictly: lambda_eff = lambda / (1 + lambda * gap)
    # So RR_eff = (lambda2 / (1 + lambda2*gap)) / (lambda1 / (1 + lambda1*gap))

    lambda1_eff_plan <- lambda1_planning / (1 + lambda1_planning * event_gap)
    lambda2_eff_plan <- lambda2_planning / (1 + lambda2_planning * event_gap)
    rate_ratio_planning <- lambda2_eff_plan / lambda1_eff_plan
  }

  lambda1_new <- lambda_est / (p1 + p2 * rate_ratio_planning)
  lambda2_new <- lambda1_new * rate_ratio_planning

  observed_info <- .blinded_info_from_tte(
    tte = df$tte,
    lambda1 = lambda1_new,
    lambda2 = lambda2_new,
    dispersion = dispersion_est,
    ratio = ratio
  )

  list(
    blinded_info = observed_info,
    dispersion_blinded = dispersion_est,
    lambda_blinded = lambda_est,
    lambda1_adjusted = lambda1_new,
    lambda2_adjusted = lambda2_new
  )
}
