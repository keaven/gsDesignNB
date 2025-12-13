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
dropout_rate <- 0.1
max_followup <- 12
trial_duration <- 24
event_gap <- 30 / 365.25 # 30 days gap

# Piecewise accrual: Ramp up
# Let's target N ~ 200
# First calculate required N with fixed accrual to get an idea
design_prelim <- sample_size_nbinom(
  lambda1 = lambda1, lambda2 = lambda2, dispersion = dispersion,
  power = power, alpha = alpha, sided = 1,
  accrual_rate = 10, accrual_duration = 12, # dummy
  trial_duration = trial_duration,
  dropout_rate = dropout_rate,
  max_followup = max_followup,
  event_gap = event_gap,
  method = "friede"
)

# Target N total
N_target <- design_prelim$n_total
# Create a piecewise accrual that achieves this N in, say, 12 months
# Rate 1 for first 6 months, Rate 2 for next 6 months
# R1 * 6 + R2 * 6 = N
# Let R2 = 2 * R1 -> 3 * R1 * 6 = N -> 18 * R1 = N -> R1 = N / 18
R1 <- N_target / 18
R2 <- 2 * R1

accrual_rate <- c(R1, R2)
accrual_duration <- c(6, 6)
# Ensure total accrual duration <= trial duration
# 12 <= 24. OK.

# Recalculate design with specific accrual to get exact power/info
design <- sample_size_nbinom(
  lambda1 = lambda1, lambda2 = lambda2, dispersion = dispersion,
  power = NULL, # Calculate power for fixed N (implied by rates)
  alpha = alpha, sided = 1,
  accrual_rate = accrual_rate,
  accrual_duration = accrual_duration,
  trial_duration = trial_duration,
  dropout_rate = dropout_rate,
  max_followup = max_followup,
  event_gap = event_gap,
  method = "friede"
)

# 2. Simulation Setup
nsim <- 10000
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

      # Analyze
      # We need to analyze using the gap-adjusted method if possible, or just standard NB?
      # The sample size formula assumes we analyze the counts.
      # If we use standard NB analysis on counts, does it account for gap?
      # Usually, gap reduces exposure.
      # We should calculate exposure at risk.

      # Calculate exposure at risk for each subject
      # Total follow-up time minus gap time for events
      # But cut_data_by_date might handle this?
      # cut_data_by_date returns the cut data.
      # We need to aggregate.

      # Aggregate by subject
      # For each subject, count events and calculate exposure.
      # Exposure = (End - Start) - (Events * Gap)?
      # Or does cut_data_by_date return intervals?
      # nb_sim returns rows per event.

      # Let's look at cut_data_by_date output structure.
      # It returns same structure as nb_sim.

      # We need to calculate total exposure at risk.
      # For a subject with k events:
      # T_total = min(T_cens, T_death, T_admin)
      # T_risk = T_total - k * gap (approx, or exact depending on when gap falls)
      # Actually, if gap is dead time, we just subtract it.

      # Let's assume simple approximation: T_risk = T_total - (N_events * gap)
      # But we must ensure T_risk >= 0.

      # Better: Sum of inter-arrival times?
      # No, we have TTE.

      # Let's use the simple aggregation:
      # Max tte per subject is the follow-up time.
      # Number of events is sum(event).

      analysis_dt <- cut_dt_dt[, .(
        events = sum(event),
        followup = max(tte),
        treatment = unique(treatment)
      ), by = id]

      # Adjust exposure for gap
      analysis_dt[, exposure := followup - events * event_gap]
      # Ensure non-negative (should be if simulation is correct)
      analysis_dt[exposure < 0, exposure := 0]

      # Fit NB model
      # glm.nb with offset(log(exposure))

      fit <- tryCatch(
        {
          MASS::glm.nb(events ~ treatment + offset(log(exposure)), data = analysis_dt)
        },
        error = function(e) NULL
      )

      if (is.null(fit)) {
        return(NULL)
      }

      coefs <- summary(fit)$coefficients
      # Treatment effect (Experimental vs Control)
      # Assuming factor levels: Control is ref?
      # We need to check levels.
      # nb_sim assigns "Control", "Experimental". Alphabetically Control comes first.

      est <- coefs["treatmentExperimental", "Estimate"]
      se <- coefs["treatmentExperimental", "Std. Error"]
      p_val <- coefs["treatmentExperimental", "Pr(>|z|)"]

      # One-sided p-value (assuming H1: Experimental < Control, i.e. est < 0)
      # If est < 0, p = p_val / 2. If est > 0, p = 1 - p_val / 2.
      p_one_sided <- if (est < 0) p_val / 2 else 1 - p_val / 2

      list(
        estimate = est,
        se = se,
        p_value = p_one_sided,
        exposure_control = mean(analysis_dt[treatment == "Control"]$exposure),
        exposure_experimental = mean(analysis_dt[treatment == "Experimental"]$exposure)
      )
    },
    error = function(e) NULL
  )
}

# 2. Simulation Setup
nsim <- 10000
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
        n = NULL # Determined by enrollment
      )

      # Cut data at trial duration (administrative censoring)
      cut_dt <- cut_data_by_date(sim_data, cut_date = trial_duration)
      cut_dt_dt <- data.table::as.data.table(cut_dt)

      # Analyze
      res <- mutze_test(cut_dt, method = "nb")

      # Calculate observed exposure
      exposure_obs <- cut_dt_dt[, .(exposure = mean(tte)), by = treatment]

      list(
        p_value = res$p_value,
        z = res$z,
        estimate = res$estimate,
        se = res$se,
        converged = !is.null(res$model),
        exposure_control = exposure_obs[treatment == "Control", exposure],
        exposure_experimental = exposure_obs[treatment == "Experimental", exposure],
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

# Bind results
results_dt <- rbindlist(results_list)

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
