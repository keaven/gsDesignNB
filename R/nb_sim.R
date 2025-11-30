#' Simulate Recurrent Events with Fixed Follow-up
#'
#' Simulates recurrent events for a clinical trial with piecewise constant enrollment,
#' exponential failure rates (Poisson process), and piecewise exponential dropout.
#'
#' @param enroll_rate A data frame with columns \code{rate} and \code{duration} defining
#'   the piecewise constant enrollment rates.
#' @param fail_rate A data frame with columns \code{treatment} and \code{rate} defining
#'   the exponential failure rate for each treatment group.
#' @param dropout_rate A data frame with columns \code{treatment}, \code{rate}, and \code{duration}
#'   defining the piecewise constant dropout rates.
#' @param max_followup Numeric. Maximum duration of follow-up for each individual
#'   (relative to their randomization time).
#' @param n Total sample size. If NULL, it is estimated from \code{enroll_rate}.
#'   If provided, enrollment stops when \code{n} subjects are recruited.
#'
#' @return A data frame (tibble) with columns:
#'   \describe{
#'     \item{id}{Subject identifier}
#'     \item{treatment}{Treatment group}
#'     \item{enroll_time}{Time of enrollment relative to trial start}
#'     \item{tte}{Time to event or censoring relative to randomization}
#'     \item{calendar_time}{Calendar time of event or censoring (enroll_time + tte)}
#'     \item{event}{Binary indicator: 1 for event, 0 for censoring}
#'   }
#'   Multiple rows per subject are returned (one for each event, plus one for the final censoring time).
#'
#' @export
#'
#' @import data.table
#' @importFrom stats rexp runif
#' @importFrom utils tail
#' @importFrom simtrial rpwexp_enroll
nb_sim <- function(enroll_rate, fail_rate, dropout_rate = NULL, max_followup = NULL, n = NULL) {
  # 1. Generate Enrollment
  # Simplified implementation of piecewise constant enrollment
  # If n is provided, we simulate until n. If not, we assume enroll_rate defines the full period.

  # Validate inputs
  if (is.null(enroll_rate) || !is.data.frame(enroll_rate)) stop("enroll_rate must be a data frame")
  if (is.null(fail_rate) || !is.data.frame(fail_rate)) stop("fail_rate must be a data frame")
  if (is.null(max_followup)) stop("max_followup must be provided")

  # Calculate total duration and expected N if n not provided
  if (is.null(n)) {
    n <- round(sum(enroll_rate$rate * enroll_rate$duration))
  }

  # Generate enrollment times using simtrial helper
  enroll_times <- simtrial::rpwexp_enroll(n = n, enroll_rate = enroll_rate)
  if (length(enroll_times) < n) {
    stop("Failed to generate sufficient enrollment times. Please check enroll_rate specification.")
  }
  enroll_times <- sort(enroll_times)[seq_len(n)]

  # Assign treatments with approximately balanced allocation
  treatments <- unique(fail_rate$treatment)
  if (length(treatments) == 0) {
    stop("fail_rate must include at least one treatment")
  }
  assigned_trt <- sample(rep(treatments, length.out = n))

  # Prepare data.tables
  dt_subjects <- data.table(
    id = seq_len(n),
    treatment = assigned_trt,
    enroll_time = enroll_times
  )
  dt_fail <- data.table(fail_rate)
  setkey(dt_fail, treatment)
  dt_subjects <- dt_fail[dt_subjects, on = "treatment"]
  setnames(dt_subjects, "rate", "lambda")

  if (any(is.na(dt_subjects$lambda))) {
    stop("Each treatment must have an associated failure rate.")
  }

  dropout_dt <- NULL
  if (!is.null(dropout_rate)) {
    dropout_dt <- data.table(dropout_rate)
  }

  compute_dropout_time <- function(trt) {
    if (is.null(dropout_dt) || nrow(dropout_dt) == 0) {
      return(Inf)
    }
    dr_sub <- dropout_dt[
      is.na(treatment) | treatment == trt
    ]
    if (nrow(dr_sub) == 0) {
      dr_sub <- dropout_dt[is.na(treatment)]
    }
    if (nrow(dr_sub) == 0) {
      return(Inf)
    }
    t_curr <- 0
    for (j in seq_len(nrow(dr_sub))) {
      rate_j <- dr_sub$rate[j]
      dur_j <- dr_sub$duration[j]
      if (!is.na(rate_j) && rate_j > 0) {
        e <- rexp(1, rate_j)
        if (e < dur_j) {
          return(t_curr + e)
        }
      }
      t_curr <- t_curr + dur_j
    }
    last_rate <- tail(dr_sub$rate, 1)
    if (!is.na(last_rate) && last_rate > 0) {
      return(t_curr + rexp(1, last_rate))
    }
    Inf
  }

  simulate_subject <- function(id, treatment, enroll_time, lambda) {
    dropout_time <- compute_dropout_time(treatment)
    end_time <- min(dropout_time, max_followup)
    if (!is.finite(end_time)) {
      end_time <- max_followup
    }
    event_times <- numeric(0)
    if (!is.na(lambda) && lambda > 0 && end_time > 0) {
      cum_t <- 0
      repeat {
        gap <- rexp(1, lambda)
        cum_t <- cum_t + gap
        if (cum_t <= end_time) {
          event_times <- c(event_times, cum_t)
        } else {
          break
        }
      }
    }
    ttes <- c(event_times, end_time)
    data.table(
      id = id,
      treatment = treatment,
      enroll_time = enroll_time,
      tte = ttes,
      calendar_time = enroll_time + ttes,
      event = c(rep(1, length(event_times)), 0)
    )
  }

  result_dt <- dt_subjects[, simulate_subject(id, treatment, enroll_time, lambda), by = id]
  setorder(result_dt, id, calendar_time)
  result_df <- as.data.frame(result_dt)
  class(result_df) <- unique(c("nb_sim_data", class(result_df)))
  result_df
}
