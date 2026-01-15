# Generate group sequential simulation data for caching
# This produces the results used in group-sequential-simulation.Rmd
library(gsDesignNB)
library(gsDesign)
library(data.table)

set.seed(54321)

# Design parameters
lambda1 <- 1.5 / 12  # Control rate
lambda2 <- 1.0 / 12  # Experimental rate
dispersion <- 0.5
event_gap_val <- 20 / 30.4375
dropout_rate_val <- -log(0.95) / 12
analysis_times <- c(10, 18, 24)

cat("Step 1: Sample size calculation...\n")
nb_ss <- sample_size_nbinom(
  lambda1 = lambda1, lambda2 = lambda2, dispersion = dispersion,
  power = 0.9, alpha = 0.025, sided = 1,
  accrual_rate = 1, accrual_duration = 12,
  trial_duration = 24, dropout_rate = dropout_rate_val,
  max_followup = 12, event_gap = event_gap_val
)
cat("  n_total =", nb_ss$n_total, "\n")

cat("\nStep 2: Create GS design...\n")
gs_nb <- gsNBCalendar(
  x = nb_ss, k = 3, test.type = 4,
  analysis_times = analysis_times,
  sfu = sfLinear, # Linear spending function for upper bound
  sfupar = c(.5, .5), # Identity function
  sfl = sfHSD, # HSD spending for lower bound
  sflpar = -8, # Conservative futility bound
  usTime = c(.1, .18, 1), # Upper spending timing
  lsTime = NULL, # Spending based on information

) |> gsDesignNB::toInteger()
cat("  n per arm:", gs_nb$n.I, "\n")
cat("  timing:", gs_nb$timing, "\n")

n_target <- ceiling(nb_ss$n_total)
enroll_rate <- data.frame(rate = n_target / 12, duration = 12)
fail_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(lambda1, lambda2),
  dispersion = c(dispersion, dispersion)
)
dropout_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(-log(0.95), -log(0.95)),
  duration = c(100, 100)
)

n_sims <- 3600
cat("\nStep 3: Running", n_sims, "simulations...\n")
sim_res <- sim_gs_nbinom(
  n_sims = n_sims,
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  dropout_rate = dropout_rate,
  max_followup = 12,
  event_gap = event_gap_val,
  cuts = lapply(analysis_times, function(x) list(planned_calendar = x)),
  design = gs_nb,
  n_target = n_target
)
cat("  Simulations complete.\n")

cat("\nStep 4: Checking bounds...\n")
results <- check_gs_bound(sim_res, gs_nb, info_scale = "blinded")
cat("  Bounds checked.\n")

summary_gs <- summarize_gs_sim(results)
cat("\nStep 5: Overall Summary\n")
cat("  Power:", round(summary_gs$power * 100, 1), "%\n")
cat("  Futility:", round(summary_gs$futility * 100, 1), "%\n")

# Compute summary statistics at each analysis
dt <- as.data.table(results)

# Use trimmed mean for blinded_info to avoid extreme outliers
trimmed_mean <- function(x, trim = 0.01) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x, trim = trim)
}

analysis_summary <- dt[, .(
  n_enrolled = mean(n_enrolled, na.rm = TRUE),
  n_ctrl = mean(n_ctrl, na.rm = TRUE),
  n_exp = mean(n_exp, na.rm = TRUE),
  exposure_at_risk_ctrl = mean(exposure_at_risk_ctrl, na.rm = TRUE),
  exposure_at_risk_exp = mean(exposure_at_risk_exp, na.rm = TRUE),
  exposure_total_ctrl = mean(exposure_total_ctrl, na.rm = TRUE),
  exposure_total_exp = mean(exposure_total_exp, na.rm = TRUE),
  events_ctrl = mean(events_ctrl, na.rm = TRUE),
  events_exp = mean(events_exp, na.rm = TRUE),
  events_total = mean(events_total, na.rm = TRUE),
  blinded_info = trimmed_mean(blinded_info, trim = 0.01),
  unblinded_info = trimmed_mean(unblinded_info, trim = 0.01),
  z_stat_mean = mean(z_stat, na.rm = TRUE),
  z_stat_sd = sd(z_stat, na.rm = TRUE),
  cross_upper = mean(cross_upper, na.rm = TRUE),
  cross_lower = mean(cross_lower, na.rm = TRUE)
), by = analysis]

cat("\n=== Summary by Analysis ===\n")
print(analysis_summary)

# Save data for vignette
save_list <- list(
  nb_ss = nb_ss,
  gs_nb = gs_nb,
  sim_res = sim_res,
  results = results,
  summary_gs = summary_gs,
  analysis_summary = as.data.frame(analysis_summary),
  n_sims = n_sims,
  params = list(
    lambda1 = lambda1,
    lambda2 = lambda2,
    dispersion = dispersion,
    event_gap = event_gap_val,
    dropout_rate = dropout_rate_val,
    analysis_times = analysis_times
  )
)

saveRDS(save_list, "inst/extdata/gs_simulation_results.rds")
cat("\nResults saved to inst/extdata/gs_simulation_results.rds\n")
