# Generate targeted SSR sensitivity results comparing Wald- and score-sized
# starting designs under the same score-test final analysis.
#
# This is intentionally smaller than the main 90-scenario SSR production grid.
# It supports the recommendation to compare starting sample-size rules for an
# SSR protocol, without replacing the main SSR operating-characteristic cache.
#
# Example:
#   Rscript data-raw/generate_ssr_sizing_sensitivity.R
#   GSDESIGNNB_SSR_SENS_N_NULL=2000 GSDESIGNNB_SSR_SENS_N_POWER=1000 \
#     Rscript data-raw/generate_ssr_sizing_sensitivity.R

suppressPackageStartupMessages({
  library(data.table)
  library(gsDesign)
  library(future)
})

data.table::setDTthreads(1)

project_root <- if (file.exists("DESCRIPTION")) "." else
  if (file.exists("../DESCRIPTION")) ".." else "."

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(project_root, quiet = TRUE)
} else {
  library(gsDesignNB)
}

env_int <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (!is.finite(value) || value < 1L) default else value
}

n_sims_null <- env_int("GSDESIGNNB_SSR_SENS_N_NULL", 5000L)
n_sims_power <- env_int("GSDESIGNNB_SSR_SENS_N_POWER", 3000L)
workers_default <- max(1L, future::availableCores() - 1L)
workers <- env_int("GSDESIGNNB_WORKERS", workers_default)

extdata_dir <- file.path(project_root, "inst", "extdata")
dir.create(extdata_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(extdata_dir, "ssr_sizing_sensitivity.rds")

lambda1_plan <- 0.15
rr_plan <- 0.7
lambda2_plan <- lambda1_plan * rr_plan
k_plan <- 0.5
power_plan <- 0.9
alpha_plan <- 0.025
accrual_dur_plan <- 12
max_followup <- 12
trial_dur_plan <- accrual_dur_plan + max_followup
event_gap_val <- 0
analysis_times_plan <- c(9, 14, 24)
analysis_test_type <- "score"

min_if_futility <- 0.3
max_enrollment_frac_for_ia2 <- 1.00
min_months_to_close_for_adapt <- 2
analysis_lag_months <- 2

make_design <- function(sizing) {
  fixed <- sample_size_nbinom(
    lambda1 = lambda1_plan,
    lambda2 = lambda2_plan,
    dispersion = k_plan,
    power = power_plan,
    alpha = alpha_plan,
    accrual_rate = 30,
    accrual_duration = accrual_dur_plan,
    trial_duration = trial_dur_plan,
    max_followup = max_followup,
    event_gap = event_gap_val,
    test_type = sizing
  )

  gs <- fixed |>
    gsNBCalendar(
      k = 3,
      test.type = 4,
      alpha = alpha_plan,
      sfu = sfHSD,
      sfupar = -2,
      sfl = sfCauchy,
      sflpar = c(0.4, 0.75, 0.56, 0.63),
      analysis_times = analysis_times_plan
    ) |>
    gsDesignNB::toInteger()

  list(
    sizing = sizing,
    fixed = fixed,
    gs = gs,
    fixed_n = fixed$n_total,
    gs_n = gs$n_total[gs$k],
    n_max = 2 * gs$n_total[gs$k],
    accrual_rate = gs$n_total[gs$k] / accrual_dur_plan
  )
}

designs <- list(
  Wald = make_design("wald"),
  Score = make_design("score")
)

make_enroll_rate <- function(accrual_rate, n_max) {
  data.frame(rate = accrual_rate, duration = n_max / accrual_rate)
}

make_fail_rate <- function(rr_true) {
  data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(lambda1_plan, lambda1_plan * rr_true),
    dispersion = k_plan
  )
}

dropout_rate_sim <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0, 0),
  duration = c(100, 100)
)

run_cell <- function(design_label, rr_true, n_sims, ignore_futility, seed) {
  d <- designs[[design_label]]
  message(sprintf(
    "[%s] %s-sized start | RR=%.2f | n=%d | planned GS N=%d",
    Sys.time(), design_label, rr_true, n_sims, d$gs_n
  ))

  sim <- sim_ssr_nbinom(
    n_sims = n_sims,
    enroll_rate = make_enroll_rate(d$accrual_rate, d$n_max),
    fail_rate = make_fail_rate(rr_true),
    dropout_rate = dropout_rate_sim,
    max_followup = max_followup,
    design = d$gs,
    n_max = d$n_max,
    min_if_futility = min_if_futility,
    max_enrollment_frac_for_adapt = max_enrollment_frac_for_ia2,
    min_months_to_close_for_adapt = min_months_to_close_for_adapt,
    analysis_lag_months = analysis_lag_months,
    event_gap = event_gap_val,
    ignore_futility = ignore_futility,
    test_type = analysis_test_type,
    metadata = list(
      starting_sizing = tolower(design_label),
      analysis_test = analysis_test_type,
      rr_true = rr_true,
      lambda1_plan = lambda1_plan,
      k_plan = k_plan,
      fixed_n = d$fixed_n,
      gs_n = d$gs_n
    ),
    seed = seed
  )

  dt <- as.data.table(sim$trial_results)
  sm <- as.data.table(
    summarize_ssr_sim(
      dt,
      by = c(
        "starting_sizing", "analysis_test", "rr_true",
        "lambda1_plan", "k_plan", "fixed_n", "gs_n", "strategy"
      )
    )$trial_summary
  )
  sm[, `:=`(
    mcse = sqrt(rejection_rate * (1 - rejection_rate) / n_sims),
    ignore_futility = ignore_futility
  )]
  sm[]
}

message("=== SSR sizing sensitivity ===")
message("workers: ", workers)
message("n_sims_null: ", n_sims_null)
message("n_sims_power: ", n_sims_power)
message("output_path: ", output_path)
message(sprintf(
  "Wald GS N = %d; Score GS N = %d",
  designs$Wald$gs_n, designs$Score$gs_n
))

old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)
if (workers > 1L) {
  future::plan(future::multisession, workers = workers)
} else {
  future::plan(future::sequential)
}

t0 <- Sys.time()
cells <- list(
  run_cell("Wald", rr_true = 1.0, n_sims = n_sims_null,
           ignore_futility = TRUE, seed = 310001L),
  run_cell("Score", rr_true = 1.0, n_sims = n_sims_null,
           ignore_futility = TRUE, seed = 310002L),
  run_cell("Wald", rr_true = rr_plan, n_sims = n_sims_power,
           ignore_futility = FALSE, seed = 310003L),
  run_cell("Score", rr_true = rr_plan, n_sims = n_sims_power,
           ignore_futility = FALSE, seed = 310004L)
)
runtime_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

summary_dt <- rbindlist(cells, use.names = TRUE, fill = TRUE)
setorder(summary_dt, rr_true, starting_sizing, strategy)

design_summary <- rbindlist(lapply(designs, function(d) {
  data.table(
    starting_sizing = d$sizing,
    fixed_n = d$fixed_n,
    gs_n = d$gs_n,
    n_max = d$n_max,
    accrual_rate = d$accrual_rate
  )
}), use.names = TRUE)

saveRDS(
  list(
    summary = as.data.frame(summary_dt),
    design_summary = as.data.frame(design_summary),
    runtime_seconds = runtime_seconds,
    generated_at = as.character(Sys.time()),
    settings = list(
      lambda1_plan = lambda1_plan,
      rr_plan = rr_plan,
      lambda2_plan = lambda2_plan,
      k_plan = k_plan,
      power_plan = power_plan,
      alpha_plan = alpha_plan,
      accrual_dur_plan = accrual_dur_plan,
      max_followup = max_followup,
      trial_dur_plan = trial_dur_plan,
      event_gap_val = event_gap_val,
      analysis_times_plan = analysis_times_plan,
      analysis_test_type = analysis_test_type,
      n_sims_null = n_sims_null,
      n_sims_power = n_sims_power,
      workers = workers
    )
  ),
  output_path
)

message(sprintf("Saved: %s", output_path))
message(sprintf("Runtime: %.2f minutes", runtime_seconds / 60))
print(summary_dt)
