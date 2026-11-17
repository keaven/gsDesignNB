#' Convert average follow-up to a constant dropout hazard
#'
#' Calibrates an exponential dropout model when every participant could be
#' followed for the same duration in the absence of dropout.
#'
#' @param mean_followup Positive, finite average observed follow-up time.
#' @param max_followup Positive, finite potential follow-up time for each
#'   participant. Must be at least `mean_followup`.
#'
#' @details
#' For dropout time \eqn{D} with exponential hazard \eqn{\delta} and a
#' common follow-up cap \eqn{F}, the mean observed time is
#' \eqn{E[\min(D,F)] = (1 - \exp(-\delta F))/\delta}.
#' This function solves that equation, returning zero when the mean equals
#' the cap. Both inputs must be scalars in the same time unit.
#'
#' A mean alone does not identify the dropout distribution. This conversion
#' assumes independent exponential dropout and a common potential follow-up;
#' it does not adjust for staggered enrollment with an earlier calendar cutoff.
#' In an NB regression, the offset is the logarithm of each participant's
#' observed exposure, not the logarithm of the overall average.
#'
#' @return A non-negative dropout hazard per unit time, suitable for the
#'   `dropout_rate` argument of [sample_size_nbinom()].
#' @export
#' @examples
#' # Mean 0.9 years, with one year of potential follow-up for everyone
#' delta <- dropout_rate_from_followup(0.9, 1)
#' sample_size_nbinom(
#'   lambda1 = 0.66, lambda2 = 0.528, dispersion = 0.8, power = 0.9,
#'   accrual_rate = 100, accrual_duration = 1, trial_duration = 2,
#'   max_followup = 1, dropout_rate = delta
#' )
dropout_rate_from_followup <- function(mean_followup, max_followup) {
  for (value in list(mean_followup, max_followup)) {
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0) {
      stop("mean_followup and max_followup must be positive finite scalars.")
    }
  }
  if (mean_followup > max_followup) {
    stop("mean_followup must not exceed max_followup.")
  }
  if (mean_followup == max_followup) return(0)

  # Solve for y = delta * mean_followup on [0, 1], so the bracket
  # remains bounded regardless of the units or amount of dropout.
  fraction <- mean_followup / max_followup
  objective <- function(y) {
    if (y == 0) return(1 - fraction)
    x <- y / fraction
    -expm1(-x) / x - fraction
  }
  y <- stats::uniroot(objective, c(0, 1), tol = .Machine$double.eps)$root
  y / mean_followup
}
