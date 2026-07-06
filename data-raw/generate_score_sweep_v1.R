# Score vs Wald test comparison across scenarios
#
# Factorial grid: lambda1 x k x event_gap (3 x 3 x 3 = 27 scenarios)
# For each scenario:
#   - Power sims (RR = 0.70): n_sims_power trials
#   - Null sims  (RR = 1.00): n_sims_null  trials
# Both Wald and score z/p are recorded for every trial.
#
# Results saved to inst/extdata/score_sweep_results.rds
#
# Usage:
#   Rscript data-raw/generate_score_sweep.R [--debug]
#
# The --debug flag uses 50/200 sims instead of 3500/20000 and runs
# a 2x2x2 subset of scenarios for quick iteration.

library(data.table)
data.table::setDTthreads(1)
devtools::load_all(".")

args <- commandArgs(trailingOnly = TRUE)
debug_mode <- "--debug" %in% args

# ---- Configuration ----
if (debug_mode) {
  cat("=== DEBUG MODE: reduced grid and sim counts ===\n\n")
  n_sims_power <- 50L
  n_sims_null  <- 200L
  lambda1_grid  <- c(0.15, 0.40)
  k_grid        <- c(0.2, 0.5)
  gap_days_grid <- c(0, 30)
} else {
  n_sims_power <- 3500L
  n_sims_null  <- 20000L
  lambda1_grid  <- c(0.15, 0.40, 1.00)
  k_grid        <- c(0.2, 0.5, 1.0)
  gap_days_grid <- c(0, 15, 30)
}

# Fixed across all scenarios
rr_power     <- 0.70
alpha        <- 0.025
power_target <- 0.90
dropout_rate <- 0.1 / 12
max_followup <- 12
trial_duration <- 24
accrual_rel  <- c(1, 2)
accrual_dur  <- c(6, 6)

set.seed(20260428)

# ---- Build scenario grid ----
scenarios <- expand.grid(
  lambda1   = lambda1_grid,
  k         = k_grid,
  gap_days  = gap_days_grid,
  stringsAsFactors = FALSE
)
scenarios$gap <- scenarios$gap_days / 30.42
scenarios$lambda2_power <- scenarios$lambda1 * rr_power

cat(sprintf("Scenarios: %d | Power sims/scenario: %d | Null sims/scenario: %d\n",
            nrow(scenarios), n_sims_power, n_sims_null))
cat(sprintf("Total sims: %s\n\n",
            format(nrow(scenarios) * (n_sims_power + n_sims_null), big.mark = ",")))

# ---- Compute sample size for each scenario ----
cat("Computing sample sizes...\n")
scenarios$n_total <- NA_integer_
scenarios$accrual_rate1 <- NA_real_
scenarios$accrual_rate2 <- NA_real_

for (i in seq_len(nrow(scenarios))) {
  sc <- scenarios[i, ]
  design <- tryCatch(
    sample_size_nbinom(
      lambda1 = sc$lambda1,
      lambda2 = sc$lambda2_power,
      dispersion = sc$k,
      power = power_target,
      alpha = alpha,
      sided = 1,
      accrual_rate = accrual_rel,
      accrual_duration = accrual_dur,
      trial_duration = trial_duration,
      dropout_rate = dropout_rate,
      max_followup = max_followup,
      event_gap = sc$gap
    ),
    error = function(e) NULL
  )
  if (is.null(design)) {
    warning(sprintf("Scenario %d failed: lambda1=%.2f k=%.1f gap=%d",
                    i, sc$lambda1, sc$k, sc$gap_days))
  } else {
    scenarios$n_total[i] <- design$n_total
    scenarios$accrual_rate1[i] <- design$accrual_rate[1]
    scenarios$accrual_rate2[i] <- design$accrual_rate[2]
  }
}

cat("\nScenario grid:\n")
print(scenarios[, c("lambda1", "k", "gap_days", "n_total")])
cat("\n")

# Drop scenarios where sample size failed
scenarios <- scenarios[!is.na(scenarios$n_total), ]

# ---- Simulation helpers ----

run_one_trial <- function(enroll_rate, fail_rate, dropout_rate_df,
                          max_followup, event_gap, trial_duration) {
  sim <- nb_sim(
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    dropout_rate = dropout_rate_df,
    max_followup = max_followup,
    event_gap = event_gap
  )
  cut <- cut_data_by_date(sim, cut_date = trial_duration, event_gap = event_gap)

  w <- tryCatch(mutze_test(cut, test_type = "wald", sided = 1),
                error = function(e) NULL)
  s <- tryCatch(mutze_test(cut, test_type = "score", sided = 1),
                error = function(e) NULL)

  list(
    z_wald   = if (!is.null(w)) w$z else NA_real_,
    p_wald   = if (!is.null(w)) w$p_value else NA_real_,
    z_score  = if (!is.null(s)) s$z else NA_real_,
    p_score  = if (!is.null(s)) s$p_value else NA_real_,
    fallback_wald  = if (!is.null(w)) w$fallback else NA_character_,
    fallback_score = if (!is.null(s)) s$fallback else NA_character_
  )
}

run_scenario <- function(sc, lambda2, n_sims) {
  enroll_rate <- data.frame(
    rate = c(sc$accrual_rate1, sc$accrual_rate2),
    duration = accrual_dur
  )
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(sc$lambda1, lambda2),
    dispersion = c(sc$k, sc$k)
  )
  dropout_rate_df <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(dropout_rate, dropout_rate),
    duration = c(100, 100)
  )

  has_future <- requireNamespace("future.apply", quietly = TRUE)
  if (has_future) {
    res_list <- future.apply::future_lapply(
      seq_len(n_sims),
      function(i) {
        run_one_trial(enroll_rate, fail_rate, dropout_rate_df,
                      max_followup, sc$gap, trial_duration)
      },
      future.seed = TRUE
    )
  } else {
    res_list <- lapply(seq_len(n_sims), function(i) {
      run_one_trial(enroll_rate, fail_rate, dropout_rate_df,
                    max_followup, sc$gap, trial_duration)
    })
  }

  rbindlist(lapply(res_list, as.data.frame), fill = TRUE)
}

# ---- Set up parallelism ----
if (requireNamespace("future", quietly = TRUE) &&
    requireNamespace("future.apply", quietly = TRUE)) {
  n_workers <- max(1L, parallel::detectCores() - 2L)
  # Use multicore (forking) so workers inherit devtools::load_all() environment.
  # Falls back to multisession on Windows where forking is unsupported.
  if (.Platform$OS.type == "unix") {
    future::plan(future::multicore, workers = n_workers)
  } else {
    future::plan(future::multisession, workers = n_workers)
  }
  cat(sprintf("Parallel execution: %d workers (%s)\n\n",
              n_workers,
              if (.Platform$OS.type == "unix") "multicore" else "multisession"))
} else {
  cat("Sequential execution (install future + future.apply for parallel)\n\n")
}

# ---- Run simulations ----
power_results_list <- vector("list", nrow(scenarios))
null_results_list  <- vector("list", nrow(scenarios))

total_start <- proc.time()

for (i in seq_len(nrow(scenarios))) {
  sc <- scenarios[i, ]
  scenario_label <- sprintf("l1=%.2f k=%.1f gap=%dd n=%d",
                            sc$lambda1, sc$k, sc$gap_days, sc$n_total)

  # Power simulation (RR = 0.70)
  cat(sprintf("[%d/%d] Power: %s (%d sims)...",
              i, nrow(scenarios), scenario_label, n_sims_power))
  t0 <- proc.time()
  power_results_list[[i]] <- run_scenario(sc, sc$lambda2_power, n_sims_power)
  power_results_list[[i]]$scenario_id <- i
  elapsed <- (proc.time() - t0)[3]
  cat(sprintf(" %.1fs\n", elapsed))

  # Null simulation (RR = 1.00)
  cat(sprintf("[%d/%d] Null:  %s (%d sims)...",
              i, nrow(scenarios), scenario_label, n_sims_null))
  t0 <- proc.time()
  null_results_list[[i]] <- run_scenario(sc, sc$lambda1, n_sims_null)
  null_results_list[[i]]$scenario_id <- i
  elapsed <- (proc.time() - t0)[3]
  cat(sprintf(" %.1fs\n", elapsed))
}

total_elapsed <- (proc.time() - total_start)[3]
cat(sprintf("\nTotal time: %.1f min\n", total_elapsed / 60))

# ---- Summarize ----
power_raw <- rbindlist(power_results_list)
null_raw  <- rbindlist(null_results_list)

summarize_scenario <- function(dt, alpha_level) {
  dt[, .(
    n_sims        = .N,
    n_valid_wald  = sum(!is.na(p_wald)),
    n_valid_score = sum(!is.na(p_score)),
    reject_wald   = sum(p_wald < alpha_level, na.rm = TRUE),
    reject_score  = sum(p_score < alpha_level, na.rm = TRUE),
    rate_wald     = mean(p_wald < alpha_level, na.rm = TRUE),
    rate_score    = mean(p_score < alpha_level, na.rm = TRUE),
    se_wald       = sqrt(mean(p_wald < alpha_level, na.rm = TRUE) *
                         (1 - mean(p_wald < alpha_level, na.rm = TRUE)) /
                         sum(!is.na(p_wald))),
    se_score      = sqrt(mean(p_score < alpha_level, na.rm = TRUE) *
                         (1 - mean(p_score < alpha_level, na.rm = TRUE)) /
                         sum(!is.na(p_score))),
    mean_z_wald   = mean(z_wald, na.rm = TRUE),
    mean_z_score  = mean(z_score, na.rm = TRUE),
    sd_z_wald     = sd(z_wald, na.rm = TRUE),
    sd_z_score    = sd(z_score, na.rm = TRUE),
    pct_fallback_poisson_wald  = 100 * mean(fallback_wald == "poisson", na.rm = TRUE),
    pct_fallback_mom_wald      = 100 * mean(fallback_wald == "mom", na.rm = TRUE),
    pct_fallback_poisson_score = 100 * mean(fallback_score == "poisson", na.rm = TRUE),
    pct_fallback_mom_score     = 100 * mean(fallback_score == "mom", na.rm = TRUE)
  ), by = scenario_id]
}

power_summary <- summarize_scenario(power_raw, alpha)
null_summary  <- summarize_scenario(null_raw, alpha)

# Merge scenario info
sc_info <- data.frame(
  scenario_id = seq_len(nrow(scenarios)),
  scenarios[, c("lambda1", "k", "gap_days", "n_total")]
)
power_summary <- merge(sc_info, power_summary, by = "scenario_id")
null_summary  <- merge(sc_info, null_summary, by = "scenario_id")

# Subsample z-values for density plots (max 2000 per scenario)
max_z_sample <- 2000L
z_sample_power <- power_raw[, {
  idx <- if (.N > max_z_sample) sample(.N, max_z_sample) else seq_len(.N)
  .SD[idx, .(z_wald, z_score)]
}, by = scenario_id]

z_sample_null <- null_raw[, {
  idx <- if (.N > max_z_sample) sample(.N, max_z_sample) else seq_len(.N)
  .SD[idx, .(z_wald, z_score)]
}, by = scenario_id]

# ---- Save ----
output <- list(
  scenarios      = scenarios,
  power_summary  = as.data.frame(power_summary),
  null_summary   = as.data.frame(null_summary),
  z_sample_power = as.data.frame(z_sample_power),
  z_sample_null  = as.data.frame(z_sample_null),
  config = list(
    n_sims_power = n_sims_power,
    n_sims_null  = n_sims_null,
    rr_power     = rr_power,
    alpha        = alpha,
    power_target = power_target,
    debug_mode   = debug_mode
  )
)

output_path <- file.path("inst", "extdata", "score_sweep_results.rds")
saveRDS(output, file = output_path)
cat(sprintf("\nResults saved to %s\n", output_path))

# ---- Print summary ----
cat("\n=== POWER (RR = 0.70) ===\n")
print(power_summary[, c("lambda1", "k", "gap_days", "n_total",
                         "rate_wald", "se_wald", "rate_score", "se_score")])

cat("\n=== TYPE I ERROR (RR = 1.00) ===\n")
print(null_summary[, c("lambda1", "k", "gap_days", "n_total",
                        "rate_wald", "se_wald", "rate_score", "se_score")])

# ---- Clean up parallel ----
if (requireNamespace("future", quietly = TRUE)) {
  future::plan(future::sequential)
}

cat("\nDone.\n")
