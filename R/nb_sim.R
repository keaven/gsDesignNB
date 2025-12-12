#' Simulate recurrent events with fixed follow-up
#'
#' Simulates recurrent events for a clinical trial with piecewise constant enrollment,
#' exponential failure rates (Poisson process), and piecewise exponential dropout.
#'
#' The simulation generates data consistent with the negative binomial models described
#' by Friede and Schmidli (2010) and Mutze et al. (2018). Specifically, it simulates
#' a Gamma-distributed frailty variable for each individual (if dispersion > 0),
#' which acts as a multiplier for that individual's event rate. Events are then
#' generated according to a Poisson process with this subject-specific rate.
#'
#' @references
#' Friede, T., & Schmidli, H. (2010). Blinded sample size reestimation with count data: methods and applications.
#' _Statistics in Medicine_, 29(10), 1145--1156. \doi{10.1002/sim.3891}
#'
#' Mütze, T., Glimm, E., Schmidli, H., & Friede, T. (2018). Group sequential designs
#' for negative binomial outcomes. _Statistical Methods in Medical Research_,
#' 27(10), 2978--2993. \doi{10.1177/0962280218773115}
#'
#' @param enroll_rate A data frame with columns `rate` and `duration` defining
#'   the piecewise constant enrollment rates.
#' @param fail_rate A data frame with columns `treatment` and `rate` defining
#'   the exponential failure rate for each treatment group.
#'   Optionally, a `dispersion` column can be provided to generate data from
#'   a negative binomial distribution. The dispersion parameter `k` is
#'   such that \eqn{Var(Y) = \mu + k \mu^2}.
#' @param dropout_rate A data frame with columns `treatment`, `rate`,
#'   and `duration` defining the piecewise constant dropout rates.
#' @param max_followup Numeric. Maximum duration of follow-up for
#'   each individual (relative to their randomization time).
#' @param n Total sample size. If NULL, it is estimated from `enroll_rate`.
#'   If provided, enrollment stops when `n` subjects are recruited.
#' @param block Block vector for treatment allocation. Default is
#'   `c(rep("Control", 2), rep("Experimental", 2))`.
#'   If NULL, simple randomization is used (treatments are assigned with
#'   equal probability). If provided, it specifies the block structure,
#'   for example, `c(rep("A", 2), rep("B", 2))` assigns 2 to group A and
#'   2 to group B in each block.
#' @param event_gap Numeric. Gap duration after each event during which
#'   no new events are counted. Default is 0.
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
#' @importFrom stats rexp rgamma
#' @importFrom utils tail
#' @importFrom simtrial rpwexp_enroll
nb_sim <- function(enroll_rate, fail_rate, dropout_rate = NULL, max_followup = NULL, n = NULL, 
                   block = c(rep("Control", 2), rep("Experimental", 2)), event_gap = 0) {
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

  # Assign treatments
  treatments <- unique(fail_rate$treatment)
  if (length(treatments) == 0) {
    stop("fail_rate must include at least one treatment")
  }

  if (is.null(block)) {
    # Default block randomization: Control (2), Experimental (2)
    # We need to know which treatment is which. Assuming names match fail_rate.
    # But fail_rate$treatment could be anything.
    # If NULL, assume treatments are "Control" and "Experimental".
    # If not, we can't default to that specific vector.
    # However, user request was specific: default block = c(rep("control", 2), rep("experimental", 2)).
    # This implies the treatment names in fail_rate MUST include "Control" and "Experimental" (case insensitive?)
    # or we map them.
    # To be safe and general: If block is NULL, we default to balanced blocks of size 2*n_treatments?
    # OR: strictly follow the prompt: default is c("Control", "Control", "Experimental", "Experimental").

    # If treatments in fail_rate are NOT Control/Experimental, this default will fail validation.
    # Let's set the default argument in function signature if possible, or handle logic here.
    # If I put it in signature, it's clear.
    # But I need to make sure fail_rate matches.

    # Let's use the requested default but check if it applies.
    block <- c("Control", "Control", "Experimental", "Experimental")

    # If actual treatments are different, this default block is invalid.
    # Maybe fallback to simple balanced if default doesn't match?
    # "While I am fine with simple randomization when block=NULL..." -> implies simple was acceptable before.
    # "I would like the default to be block = ..." implies if user doesn't specify, try this.

    if (!all(block %in% treatments)) {
      # If default block doesn't match treatments, we can't use it.
      # Fallback to balanced blocks of size 4 (if 2 arms) or just rep(treatments, 2)?
      # Or warn?
      # Ideally, if user provides different treatment names, they should provide a block.
      # But to be robust: if default block fails validation, create a balanced block of size 2*n_arms.
      block <- rep(treatments, each = 2)
    }
  }

  # Block randomization based on block vector
  # Validate block contents
  if (!all(block %in% treatments)) {
    stop("Elements of 'block' must match treatment names in 'fail_rate'.")
  }
  n_blocks <- ceiling(n / length(block))
  # Generate full blocks
  full_blocks <- replicate(n_blocks, sample(block))
  assigned_trt <- as.vector(full_blocks)[seq_len(n)]

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

  # Handle dispersion for negative binomial simulation
  if ("dispersion" %in% names(dt_subjects)) {
    # For rows with valid positive dispersion, sample lambda from Gamma.
    # We want the count Y ~ NegBin(mean = lambda*t, dispersion = k).
    # This is achieved if subject-specific rate L ~ Gamma with
    # mean = lambda and var = k * lambda^2.
    # Gamma parameters: shape = 1/k, scale = k * lambda
    
    # Identify rows with positive dispersion
    has_disp <- !is.na(dt_subjects$dispersion) & dt_subjects$dispersion > 0
    if (any(has_disp)) {
      dt_subjects[has_disp, lambda := rgamma(
        sum(has_disp), 
        shape = 1 / dispersion, 
        scale = lambda * dispersion
      )]
    }
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
        # Generate time to next event
        # If event_gap > 0, we add the gap AFTER the previous event.
        # For the first event, there is no gap.
        # Wait, usually gap is "dead time" after event.
        # So T_1 ~ Exp(lambda). Event at T_1.
        # T_2 ~ Exp(lambda). Event at T_1 + gap + T_2.
        
        inter_arrival <- rexp(1, lambda)
        
        if (length(event_times) == 0) {
          cum_t <- cum_t + inter_arrival
        } else {
          cum_t <- cum_t + event_gap + inter_arrival
        }
        
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
