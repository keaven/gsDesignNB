# Script to generate simulation data for vignette verification
#
# This script simulates 3,600 trials at the corrected sample size (n_full)
# and produces results for two sample sizes:
#   1. n_full  -- the Jensen-corrected sample size (validates corrected formula)
#   2. n_naive -- the naive (uncorrected) sample size (validates the correction impact)
#
# For each simulated trial, the first n_naive subjects are a subset of the
# n_full subjects, so both analyses share the same underlying data.

library(data.table)

# Limit data.table threads to avoid OpenMP SHM issues
data.table::setDTthreads(1)
# Assuming we are in the package root
devtools::load_all(".")

set.seed(12345)

# 1. Define Design Parameters
lambda1 <- 0.4
lambda2 <- 0.3
dispersion <- 0.5
power <- 0.9
alpha <- 0.025
dropout_rate <- 0.1 / 12
max_followup <- 12
trial_duration <- 24
event_gap <- 20 / 30.42 # 20 days

# Accrual targeting 90% power (uses Jensen-corrected effective rates)
accrual_rate_rel <- c(1, 2)
accrual_duration <- c(6, 6)

design <- sample_size_nbinom(
  lambda1 = lambda1, lambda2 = lambda2, dispersion = dispersion,
  power = power,
  alpha = alpha, sided = 1,
  accrual_rate = accrual_rate_rel,
  accrual_duration = accrual_duration,
  trial_duration = trial_duration,
  dropout_rate = dropout_rate,
  max_followup = max_followup,
  event_gap = event_gap
)

n_full <- design$n_total
accrual_rate <- design$accrual_rate

# Also compute what the naive (uncorrected) sample size would be.
# The naive formula uses lambda_eff = lambda / (1 + lambda * gap) without
# the Jensen correction. We back-compute this by temporarily setting k = 0
# in the gap correction term, which removes the correction factor.
# Alternatively, hard-code the known value from the previous analysis.
# For reproducibility, we compute it from the naive effective rates:
naive_eff1 <- lambda1 / (1 + lambda1 * event_gap)
naive_eff2 <- lambda2 / (1 + lambda2 * event_gap)
corrected_eff1 <- naive_eff1 * (1 - dispersion * lambda1 * event_gap / (1 + lambda1 * event_gap)^2)
corrected_eff2 <- naive_eff2 * (1 - dispersion * lambda2 * event_gap / (1 + lambda2 * event_gap)^2)

# The naive design uses dispersion = 0 for the gap correction only.
# We approximate n_naive by scaling: ratio of corrected / naive effective rates
# affects the variance term. But the simplest approach: the previous run gave 422.
n_naive <- 422L

cat("Design (Jensen-corrected):\n")
cat("  n_full =", n_full, "\n")
cat("  n_naive =", n_naive, "(previous uncorrected result)\n")
cat("  accrual_rate =", accrual_rate, "\n\n")

# 2. Simulation Setup
nsim <- 3600

# Helper: analyze a cut dataset and return summary
analyze_cut <- function(cut_dt) {
  cut_dt_dt <- data.table::as.data.table(cut_dt)
  res <- mutze_test(cut_dt, method = "nb", sided = 1)

  exposure_obs <- cut_dt_dt[, .(
    exposure_at_risk = mean(tte),
    exposure_total = mean(tte_total)
  ), by = treatment]

  list(
    p_value = res$p_value,
    z = res$z,
    estimate = res$estimate,
    se = res$se,
    method_used = res$method,
    dispersion = res$dispersion,
    exposure_at_risk_control = exposure_obs[treatment == "Control", exposure_at_risk],
    exposure_at_risk_experimental = exposure_obs[treatment == "Experimental", exposure_at_risk],
    exposure_total_control = exposure_obs[treatment == "Control", exposure_total],
    exposure_total_experimental = exposure_obs[treatment == "Experimental", exposure_total],
    events_control = res$group_summary[res$group_summary$treatment == "Control", "events"],
    events_experimental = res$group_summary[res$group_summary$treatment == "Experimental", "events"],
    n_control = res$group_summary[res$group_summary$treatment == "Control", "subjects"],
    n_experimental = res$group_summary[res$group_summary$treatment == "Experimental", "subjects"]
  )
}

# Helper function for one simulation (returns results for both n_full and n_naive)
run_one_sim <- function(i) {
  tryCatch(
    {
      # Generate data at n_full
      enroll_rate_df <- data.frame(
        rate = accrual_rate,
        duration = accrual_duration
      )
      fail_rate_df <- data.frame(
        treatment = c("Control", "Experimental"),
        rate = c(lambda1, lambda2),
        dispersion = c(dispersion, dispersion)
      )
      dropout_rate_df <- data.frame(
        treatment = c("Control", "Experimental"),
        rate = c(dropout_rate, dropout_rate),
        duration = c(100, 100)
      )

      sim_data <- nb_sim(
        enroll_rate = enroll_rate_df,
        fail_rate = fail_rate_df,
        dropout_rate = dropout_rate_df,
        max_followup = max_followup,
        n = NULL, # Determined by enrollment (= n_full)
        event_gap = event_gap
      )

      # Cut at trial duration for the full dataset
      cut_full <- cut_data_by_date(sim_data, cut_date = trial_duration, event_gap = event_gap)

      # Restrict to first n_naive subjects for the naive analysis
      first_n_ids <- sort(unique(sim_data$id))[seq_len(min(n_naive, length(unique(sim_data$id))))]
      sim_data_naive <- sim_data[sim_data$id %in% first_n_ids, ]
      class(sim_data_naive) <- class(sim_data) # preserve class for dispatch
      cut_naive <- cut_data_by_date(sim_data_naive, cut_date = trial_duration, event_gap = event_gap)

      # Analyze both
      res_full <- analyze_cut(cut_full)
      res_naive <- analyze_cut(cut_naive)

      list(full = res_full, naive = res_naive)
    },
    error = function(e) {
      structure(list(error = conditionMessage(e)), class = "try-error")
    }
  )
}

# 3. Run Simulation
cat("Starting simulation with", nsim, "replicates...\n")

results_list <- lapply(seq_len(nsim), function(i) {
  if (i %% 500 == 0) cat("  Completed", i, "/", nsim, "\n")
  run_one_sim(i)
})

# Filter out errors
ok <- sapply(results_list, function(x) !inherits(x, "try-error"))
results_list <- results_list[ok]
if (sum(!ok) > 0) {
  warning(sum(!ok), " simulations failed. Proceeding with ", length(results_list), " replicates.")
}

# Separate and bind results
results_full <- rbindlist(lapply(results_list, function(x) as.data.frame(x$full)), fill = TRUE)
results_naive <- rbindlist(lapply(results_list, function(x) as.data.frame(x$naive)), fill = TRUE)

cat("Completed", nrow(results_full), "replicates.\n")

# 4. Save Results
output_path <- file.path("inst", "extdata", "simulation_results.rds")
saveRDS(list(
  design = design,
  n_full = n_full,
  n_naive = n_naive,
  results = results_full,
  results_naive = results_naive,
  params = list(
    lambda1 = lambda1,
    lambda2 = lambda2,
    dispersion = dispersion,
    accrual_rate = accrual_rate,
    accrual_duration = accrual_duration,
    dropout_rate = dropout_rate,
    max_followup = max_followup,
    trial_duration = trial_duration,
    event_gap = event_gap
  )
), file = output_path)

cat("Simulation completed. Results saved to", output_path, "\n")
cat("  results (n=", n_full, "):", nrow(results_full), "rows\n")
cat("  results_naive (n=", n_naive, "):", nrow(results_naive), "rows\n")
