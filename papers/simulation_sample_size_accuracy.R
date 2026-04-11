# Supporting survival-validation script: accuracy of Lachin–Foulkes / nSurv() expected events
#
# Part of the gsDesignNB manuscript materials (see papers/README-simulation.md).
# Hosted in this repository as supporting validation of shared timing
# infrastructure; the main package-facing workflow remains negative binomial
# recurrent-event design and simulation in gsDesignNB.
#
# Usage (from gsDesignNB repository root, with gsDesign installed):
#   Rscript papers/simulation_sample_size_accuracy.R
#
# Or from R:
#   source("papers/simulation_sample_size_accuracy.R")
#
# The analytical target is nSurv()$d (total expected events under H1), which
# matches eEvents()-based calculations. This script compares that to the mean
# number of failure events in independent trial simulations under the same
# proportional-hazards + uniform enrollment + exponential dropout model.

suppressPackageStartupMessages({
  if (!requireNamespace("gsDesign", quietly = TRUE)) {
    stop("Install gsDesign first, e.g. install.packages('gsDesign') or devtools::load_all('.')")
  }
  library(gsDesign)
})

#' Simulate total failure events by planned analysis time T (one trial)
#'
#' Competing risks: failure vs dropout; only failures before calendar time T
#' count. Enrollment uniform on [0, sum(R)]. Treatment assignment with
#' randomization ratio experimental:control = `x$ratio` : 1.
#'
#' @param x Object of class `nSurv` from [nSurv()].
#' @param n_subj Integer sample size (default `round(x$n)`).
#' @param seed Integer seed for RNG.
#' @return Integer total event count.
simulate_lf_total_events <- function(x, n_subj = NULL, seed = 1L) {
  if (is.null(n_subj)) {
    n_subj <- round(x$n)
  }
  ratio <- x$ratio
  R <- x$R
  T <- x$T
  lamC <- as.numeric(x$lambdaC)[1L]
  hr <- x$hr
  etaC <- as.numeric(x$etaC)[1L]
  etaE <- as.numeric(x$etaE)[1L]
  accr <- sum(R)
  set.seed(seed)
  enroll <- stats::runif(n_subj, 0, accr)
  is_exp <- stats::rbinom(n_subj, 1L, prob = ratio / (1 + ratio)) == 1L
  lam_fail <- ifelse(is_exp, lamC * hr, lamC)
  eta <- ifelse(is_exp, etaE, etaC)
  Tf <- stats::rexp(n_subj, rate = lam_fail)
  Td <- stats::rexp(n_subj, rate = eta)
  te <- pmin(Tf, Td)
  fail_first <- Tf <= Td
  cal <- enroll + te
  as.integer(sum(fail_first & (cal <= T)))
}

#' Monte Carlo evaluation: mean simulated events vs analytical `x$d`
#'
#' @param x `nSurv` object.
#' @param nsim Number of simulated trials.
#' @param n_subj Per-trial sample size (default `round(x$n)`).
#' @return A list with analytical and simulation summaries.
evaluate_nSurv_event_accuracy <- function(x, nsim = 2000L, n_subj = NULL) {
  if (is.null(n_subj)) {
    n_subj <- round(x$n)
  }
  d_analytical <- unname(x$d)
  vals <- vapply(
    seq_len(nsim),
    function(i) simulate_lf_total_events(x, n_subj = n_subj, seed = i),
    integer(1)
  )
  m <- mean(vals)
  s <- stats::sd(vals)
  bias <- m - d_analytical
  rel_bias <- if (d_analytical > 0) bias / d_analytical else NA_real_
  mc_se <- s / sqrt(nsim)
  z_975 <- stats::qnorm(0.975)
  list(
    analytical_d = d_analytical,
    n_subj = n_subj,
    nsim = nsim,
    mean_simulated_events = m,
    sd_simulated_events = s,
    mc_se_mean = mc_se,
    bias = bias,
    relative_bias = rel_bias,
    ci95_mean = m + c(-1, 1) * z_975 * mc_se
  )
}

#' Optional: thin recurrent events per patient (manuscript Section 2.4)
#'
#' Not used in the default LF first-failure simulation above. If you simulate
#' multiple potential failure times per patient, pass sorted times and keep
#' times separated by at least `delta` after the previous *counted* time.
#'
#' @param times Sorted numeric vector of potential event times for one patient.
#' @param delta Minimum gap (same time units as `times`).
#' @return Sorted vector of counted event times.
thin_interevent_gap <- function(times, delta) {
  if (length(times) == 0L) {
    return(numeric(0))
  }
  out <- times[1L]
  for (i in seq_len(length(times) - 1L) + 1L) {
    if (times[i] >= max(out) + delta) {
      out <- c(out, times[i])
    }
  }
  out
}

# ---- Example designs --------------------------------------------------------

scenario_scalar <- function() {
  nSurv(
    lambdaC = 0.002,
    eta = 5e-4,
    etaE = 5e-4,
    gamma = 10,
    R = 18,
    S = NULL,
    T = 24,
    minfup = 8,
    hr = 0.7,
    hr0 = 1,
    ratio = 1,
    alpha = 0.025,
    beta = 0.1
  )
}

scenario_piecewise <- function() {
  nSurv(
    lambdaC = c(0.003, 0.002),
    eta = 5e-4,
    etaE = 5e-4,
    gamma = c(8, 12),
    R = c(10, 8),
    S = 12,
    T = 30,
    minfup = 6,
    hr = 0.75,
    hr0 = 1,
    ratio = 2,
    alpha = 0.025,
    beta = 0.15
  )
}

# ---- Report -----------------------------------------------------------------

print_accuracy <- function(label, res) {
  cat("\n=== ", label, " ===\n", sep = "")
  cat(sprintf("Analytical E[events] (nSurv $d):     %.4f\n", res$analytical_d))
  cat(sprintf("Monte Carlo mean events (n = %d):   %.4f\n", res$nsim, res$mean_simulated_events))
  cat(sprintf("MC SE of mean:                      %.4f\n", res$mc_se_mean))
  cat(sprintf("95%% CI for E[events] from MC:      [%.4f, %.4f]\n",
              res$ci95_mean[1], res$ci95_mean[2]))
  cat(sprintf("Bias (MC mean - analytical):        %.4f\n", res$bias))
  cat(sprintf("Relative bias:                      %.2f%%\n", 100 * res$relative_bias))
  invisible(res)
}

run_all <- function(nsim = 2000L) {
  cat("gsDesignNB papers: simulation check for gsDesign::nSurv() expected total events\n")
  cat("(see papers/manuscript-gaps-events-sample-size.md, Section 5)\n\n")

  x1 <- scenario_scalar()
  r1 <- evaluate_nSurv_event_accuracy(x1, nsim = nsim)
  print_accuracy("Scenario A: scalar hazards, uniform accrual", r1)

  x2 <- scenario_piecewise()
  # Same competing-risk simulator applies only when hazards are constant per arm;
  # piecewise-by-calendar-time from randomization requires a different simulator.
  # Here we still report analytical d and note validation path.
  cat("\n=== Scenario B: piecewise (note) ===\n")
  cat("Analytical E[events] (nSurv $d): ", x2$d, "\n", sep = "")
  cat(
    "Piecewise-by-time failure hazards need interval-based simulation;\n",
    "compare to eEvents() strata totals or extend simulate_lf_total_events().\n",
    "Scalar Scenario A already validates the core LF event expectation.\n",
    sep = ""
  )
  invisible(list(scalar = r1, piecewise_design = x2))
}

# Run batch report when executed non-interactively (e.g. Rscript); not on source() in a session
if (!interactive()) {
  nsim <- as.integer(Sys.getenv("SIM_NSIM", unset = "2000"))
  run_all(nsim = nsim)
}
