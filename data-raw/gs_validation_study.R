#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
production <- "--production" %in% args
debug_mode <- "--debug" %in% args
if (production == debug_mode) {
  stop("Specify exactly one of --debug or --production.", call. = FALSE)
}

output_arg <- grep("^--output-dir=", args, value = TRUE)
if (length(output_arg) > 1L) {
  stop("Specify --output-dir at most once.", call. = FALSE)
}
output_dir <- if (length(output_arg) == 1L) {
  sub("^--output-dir=", "", output_arg)
} else if (production) {
  file.path("paper", "gs-validation-results")
} else {
  file.path(tempdir(), "gs-validation-debug")
}

devtools::load_all(".", quiet = TRUE)
data.table::setDTthreads(1L)

scenario_grid <- data.frame(
  scenario = c(
    "reference", "low-rate-ramp", "high-dispersion", "low-rate-gap",
    "moderate-gap", "low-rate-front-loaded", "high-rate-front-loaded",
    "mixed-stress"
  ),
  lambda1 = c(0.5, 0.15, 0.5, 0.15, 0.5, 0.15, 1.0, 0.3),
  dispersion = c(0.2, 0.2, 1.0, 1.0, 0.2, 1.0, 0.2, 0.5),
  gap_days = c(0, 0, 0, 30, 30, 0, 15, 15),
  dropout_annual = c(0, 0.1, 0.1, 0, 0.1, 0.1, 0, 0.05),
  accrual_pattern = c(
    "constant", "ramp", "ramp", "constant", "constant", "front_loaded",
    "front_loaded", "ramp"
  ),
  stringsAsFactors = FALSE
)

if (debug_mode) {
  scenario_grid <- scenario_grid[c(1, 8), ]
  n_sims_null <- 50L
  n_sims_alternative <- 50L
} else {
  n_sims_null <- 20000L
  n_sims_alternative <- 5000L
}

alpha <- 0.025
target_power <- 0.90
rate_ratio <- 0.70
trial_duration <- 24
max_followup <- 12
event_gap_months_per_day <- 1 / 30.4375
calendar_times <- c(10, 18, 24)
timing_methods <- c("information_fraction", "fixed_calendar")
base_seed <- 20261106L

accrual_shape <- function(pattern) {
  switch(
    pattern,
    constant = list(rate = 1, duration = 12),
    ramp = list(rate = c(1, 2), duration = c(6, 6)),
    front_loaded = list(rate = c(2, 1), duration = c(6, 6)),
    stop("Unknown accrual pattern: ", pattern, call. = FALSE)
  )
}

build_design <- function(fixed_design, timing_method) {
  common <- list(
    x = fixed_design,
    k = 3,
    test.type = 1,
    alpha = alpha,
    beta = 1 - target_power,
    sfu = gsDesign::sfHSD,
    sfupar = -4
  )
  design <- if (timing_method == "information_fraction") {
    do.call(gsDesignNB::gsNB, common)
  } else {
    do.call(gsDesignNB::gsNBCalendar, c(common, list(analysis_times = calendar_times)))
  }
  gsDesignNB::toInteger(design)
}

run_cell <- function(scenario, design, timing_method, truth, n_sims, seed) {
  shape <- accrual_shape(scenario$accrual_pattern)
  n_target <- as.integer(ceiling(utils::tail(design$n_total, 1)))
  exposure_units <- sum(shape$rate * shape$duration)
  enroll_rate <- data.frame(
    rate = shape$rate * n_target / exposure_units,
    duration = shape$duration
  )
  lambda2 <- if (truth == "null") scenario$lambda1 else scenario$lambda1 * rate_ratio
  dropout_hazard <- -log1p(-scenario$dropout_annual) / 12
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(scenario$lambda1, lambda2),
    dispersion = scenario$dispersion
  )
  dropout_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(dropout_hazard, dropout_hazard),
    duration = c(100, 100)
  )
  cuts <- lapply(design$T, function(time) list(planned_calendar = time))

  sim_results <- gsDesignNB::sim_gs_nbinom(
    n_sims = n_sims,
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    dropout_rate = dropout_rate,
    max_followup = max_followup,
    event_gap = scenario$gap_days * event_gap_months_per_day,
    cuts = cuts,
    design = design,
    n_target = n_target,
    test_type = "score",
    seed = seed
  )
  expected_information <- vapply(design$T, function(time) {
    gsDesignNB::compute_info_at_time(
      analysis_time = time,
      accrual_rate = enroll_rate$rate,
      accrual_duration = enroll_rate$duration,
      lambda1 = scenario$lambda1,
      lambda2 = lambda2,
      dispersion = scenario$dispersion,
      ratio = design$ratio,
      dropout_rate = dropout_hazard,
      event_gap = scenario$gap_days * event_gap_months_per_day,
      max_followup = max_followup
    )
  }, numeric(1))
  information_metrics <- c(
    "info_unblinded_ml", "info_blinded_ml",
    "info_unblinded_mom", "info_blinded_mom"
  )
  look_rows <- lapply(seq_len(design$k), function(analysis_index) {
    look <- sim_results[sim_results$analysis == analysis_index, ]
    metric_rows <- lapply(information_metrics, function(metric) {
      estimates <- look[[metric]]
      estimates <- estimates[is.finite(estimates)]
      n_valid <- length(estimates)
      mean_information <- if (n_valid > 0L) mean(estimates) else NA_real_
      mcse <- if (n_valid > 1L) stats::sd(estimates) / sqrt(n_valid) else NA_real_
      margin <- if (n_valid > 1L) stats::qt(0.975, n_valid - 1L) * mcse else NA_real_
      data.frame(
        scenario = scenario$scenario,
        timing_method = timing_method,
        truth = truth,
        analysis = analysis_index,
        estimator = metric,
        n_sims = nrow(look),
        n_valid = n_valid,
        n_failed = nrow(look) - n_valid,
        planned_time = design$T[analysis_index],
        mean_analysis_time = mean(look$analysis_time),
        planned_enrolled = design$n_total[analysis_index],
        mean_enrolled = mean(look$n_enrolled),
        planned_information = design$n.I[analysis_index],
        expected_information = expected_information[analysis_index],
        mean_information = mean_information,
        information_ratio = mean_information / expected_information[analysis_index],
        mcse = mcse,
        ci_lower = mean_information - margin,
        ci_upper = mean_information + margin
      )
    })
    do.call(rbind, metric_rows)
  })
  do.call(rbind, look_rows)
}
look_rows <- list()
row_id <- 0L
cell_dir <- file.path(output_dir, "cells")
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)

for (scenario_index in seq_len(nrow(scenario_grid))) {
  scenario <- scenario_grid[scenario_index, ]
  shape <- accrual_shape(scenario$accrual_pattern)
  dropout_hazard <- -log1p(-scenario$dropout_annual) / 12
  fixed_design <- gsDesignNB::sample_size_nbinom(
    lambda1 = scenario$lambda1,
    lambda2 = scenario$lambda1 * rate_ratio,
    dispersion = scenario$dispersion,
    power = target_power,
    alpha = alpha,
    sided = 1,
    accrual_rate = shape$rate,
    accrual_duration = shape$duration,
    trial_duration = trial_duration,
    dropout_rate = dropout_hazard,
    max_followup = max_followup,
    event_gap = scenario$gap_days * event_gap_months_per_day,
    test_type = "score"
  )

  for (timing_index in seq_along(timing_methods)) {
    timing_method <- timing_methods[timing_index]
    design <- build_design(fixed_design, timing_method)
    truth_values <- c("null", "alternative")
    for (truth_index in seq_along(truth_values)) {
      truth <- truth_values[truth_index]
      n_sims <- if (truth == "null") n_sims_null else n_sims_alternative
      seed <- base_seed + scenario_index * 1000L + timing_index * 100L + truth_index
      cell_file <- file.path(
        cell_dir,
        sprintf("information-v1__%s__%s__%s__%d.rds",
                scenario$scenario, timing_method, truth, n_sims)
      )
      if (file.exists(cell_file)) {
        message(sprintf(
          "Scenario %s | %s | %s | cached", scenario$scenario, timing_method, truth
        ))
        cell <- readRDS(cell_file)
      } else {
        message(
          sprintf(
            "Scenario %s | %s | %s | %d replicates",
            scenario$scenario, timing_method, truth, n_sims
          )
        )
        cell <- run_cell(
          scenario, design, timing_method, truth, n_sims, seed
        )
        saveRDS(cell, paste0(cell_file, ".tmp"))
        file.rename(paste0(cell_file, ".tmp"), cell_file)
      }
      row_id <- row_id + 1L
      look_rows[[row_id]] <- cell
    }
  }
}

information_summary <- do.call(rbind, look_rows)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  information_summary,
  file.path(output_dir, "information-summary.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    scenario_grid = scenario_grid,
    information_summary = information_summary,
    n_sims_null = n_sims_null,
    n_sims_alternative = n_sims_alternative,
    seed = base_seed
  ),
  file.path(output_dir, "gs-information-validation.rds")
)
capture.output(
  utils::sessionInfo(),
  file = file.path(output_dir, "session-info.txt")
)
print(information_summary)
message("Results saved under: ", normalizePath(output_dir))