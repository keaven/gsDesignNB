# Script to generate simulation data for vignette verification

library(parallel)
library(data.table)
library(future.apply)

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

# Accrual targeting 90% power
# We provide relative rates (1:2) and the function scales them to achieve power
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

# Extract calculated absolute accrual rates
accrual_rate <- design$accrual_rate

# 2. Simulation Setup
nsim <- 3600
n_cores <- future::availableCores()
if (is.na(n_cores) || n_cores < 1) n_cores <- 1

# Helper function for one simulation
run_one_sim <- function(i) {
  tryCatch(
    {
      # Generate data
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
        duration = c(100, 100) # Long duration
      )

      sim_data <- nb_sim(
        enroll_rate = enroll_rate_df,
        fail_rate = fail_rate_df,
        dropout_rate = dropout_rate_df,
        max_followup = max_followup,
        n = NULL, # Determined by enrollment
        event_gap = event_gap
      )

      # Cut data at trial duration (administrative censoring)
      cut_dt <- cut_data_by_date(sim_data, cut_date = trial_duration, event_gap = event_gap)
      cut_dt_dt <- data.table::as.data.table(cut_dt)

      # Analyze with mutze_test (defaults to sided=1 as updated)
      res <- mutze_test(cut_dt, method = "nb", sided = 1)

      # Calculate observed exposure metrics by treatment group
      # tte = exposure at risk (time available for events, gaps subtracted)
      # tte_total = total follow-up (calendar time on study)
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
    },
    error = function(e) {
      structure(list(error = conditionMessage(e)), class = "try-error")
    }
  )
}

# 3. Run Simulation
cat("Starting simulation with", nsim, "replicates on", n_cores, "workers (future multisession)...\n")

# Use future multisession to avoid fork/OpenMP issues
plan(multisession, workers = n_cores)
results_list <- future_lapply(
  seq_len(nsim),
  run_one_sim,
  future.seed = TRUE
)

# Filter out errors
results_list <- results_list[sapply(results_list, function(x) !inherits(x, "try-error"))]
if (length(results_list) < nsim) {
  warning("Some simulations failed. Proceeding with ", length(results_list), " replicates.")
}

# Bind results - convert to data.table properly
results_dt <- rbindlist(lapply(results_list, as.data.frame), fill = TRUE)

# 4. Save Results
output_path <- file.path("inst", "extdata", "simulation_results.rds")
saveRDS(list(
  design = design,
  results = results_dt,
  params = list(
    lambda1 = lambda1,
    lambda2 = lambda2,
    dispersion = dispersion,
    accrual_rate = accrual_rate,
    accrual_duration = accrual_duration,
    dropout_rate = dropout_rate,
    max_followup = max_followup,
    trial_duration = trial_duration
  )
), file = output_path)

cat("Simulation completed. Results saved to", output_path, "\n")
