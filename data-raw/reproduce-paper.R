#!/usr/bin/env Rscript

# Reproduce the lightweight computations shown in paper/gsDesignNB-paper.qmd.
# Run from the package root: Rscript data-raw/reproduce-paper.R

args <- commandArgs(trailingOnly = TRUE)
output_arg <- grep("^--output-dir=", args, value = TRUE)

if (length(output_arg) > 1L) {
  stop("Specify --output-dir at most once.", call. = FALSE)
}

output_dir <- if (length(output_arg) == 1L) {
  sub("^--output-dir=", "", output_arg)
} else {
  file.path("paper", "replication-output")
}

if (!requireNamespace("gsDesignNB", quietly = TRUE)) {
  stop(
    "Install gsDesignNB first, then run this script from the package root.",
    call. = FALSE
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_table <- function(x, name) {
  utils::write.csv(
    as.data.frame(x),
    file.path(output_dir, name),
    row.names = FALSE
  )
}

# Worked design example from the manuscript (Section 2.2).
worked_fixed <- gsDesignNB::sample_size_nbinom(
  lambda1 = 0.5,
  lambda2 = 0.35,
  dispersion = 0.5,
  power = 0.90,
  alpha = 0.025,
  accrual_rate = c(12, 24),
  accrual_duration = c(6, 6),
  trial_duration = 24,
  dropout_rate = 0.10 / 12,
  max_followup = 12,
  event_gap = 20 / 30.4375,
  test_type = "score"
)

worked_gs <- worked_fixed |>
  gsDesignNB::gsNBCalendar(
    k = 3,
    test.type = 4,
    alpha = 0.025,
    sfu = gsDesign::sfHSD,
    sfupar = -2,
    sfl = gsDesign::sfCauchy,
    sflpar = c(0.4, 0.75, 0.56, 0.63),
    analysis_times = c(9, 14, 24)
  ) |>
  gsDesignNB::toInteger()

write_table(
  data.frame(
    quantity = c(
      "fixed_design_sample_size",
      "group_sequential_final_sample_size",
      "target_information",
      "expected_final_counted_events",
      "average_final_exposure",
      "test_type"
    ),
    value = c(
      round(worked_fixed$n_total),
      round(worked_gs$n_total[worked_gs$k]),
      round(worked_gs$n.I[worked_gs$k], 1),
      round(worked_gs$events[worked_gs$k], 1),
      round(worked_gs$exposure[worked_gs$k], 1),
      worked_fixed$test_type
    )
  ),
  "worked-design-summary.csv"
)

write_table(
  data.frame(
    look = c("IA 1", "IA 2", "Final"),
    month = worked_gs$T,
    information_fraction = round(worked_gs$timing, 3),
    cumulative_n = round(worked_gs$n_total),
    expected_events = round(worked_gs$events, 1),
    upper_z = round(worked_gs$upper$bound, 3),
    lower_z = round(worked_gs$lower$bound, 3)
  ),
  "worked-design-boundaries.csv"
)

# Table 2: closed-form Jensen-correction magnitudes.
jensen_grid <- expand.grid(
  lambda = c(1.0, 0.5, 0.3),
  dispersion = c(1.0, 0.5),
  gap_months = c(1.0, 0.5)
)
jensen_grid$correction_percent <- with(
  jensen_grid,
  round(100 * dispersion * lambda * gap_months /
    (1 + lambda * gap_months)^2, 1)
)
write_table(jensen_grid, "jensen-correction-table.csv")

# Examples 1 and 2 in the reproducible-workflows section. Example 2 depends
# on the design and rate objects created here, so both examples are run together.
enroll_rate <- data.frame(rate = 12, duration = 6)
fail_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0.6, 0.42),
  dispersion = c(0.4, 0.4)
)
dropout_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0.05, 0.05),
  duration = c(12, 12)
)

fixed_design <- gsDesignNB::sample_size_nbinom(
  lambda1 = 0.6,
  lambda2 = 0.42,
  dispersion = 0.4,
  power = 0.8,
  alpha = 0.025,
  accrual_rate = enroll_rate$rate,
  accrual_duration = enroll_rate$duration,
  trial_duration = 12,
  max_followup = 6
)

gs_design <- gsDesignNB::gsNBCalendar(
  fixed_design,
  k = 3,
  test.type = 4,
  analysis_times = c(4, 8, 12)
)

sim_gs <- gsDesignNB::sim_gs_nbinom(
  n_sims = 50,
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  dropout_rate = dropout_rate,
  max_followup = 6,
  n_target = ceiling(utils::tail(gs_design$n_total, 1)),
  design = gs_design,
  cuts = lapply(gs_design$T, function(t) list(planned_calendar = t)),
  seed = 20260413
)

gs_summary <- gsDesignNB::summarize_gs_sim(
  gsDesignNB::check_gs_bound(
    sim_gs,
    gs_design,
    info_col = "info_unblinded_ml"
  )
)
write_table(gs_summary$analysis_summary, "example-1-gs-analysis-summary.csv")

ssr_res <- gsDesignNB::sim_ssr_nbinom(
  n_sims = 50,
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  dropout_rate = dropout_rate,
  max_followup = 6,
  design = gs_design,
  strategies = c("No adaptation", "Blinded SSR", "Unblinded SSR"),
  bound_info = "unblinded_ml",
  seed = 20260413
)

ssr_summary <- gsDesignNB::summarize_ssr_sim(ssr_res)$trial_summary
write_table(
  ssr_summary[, c("strategy", "rejection_rate", "mean_n", "pct_adapted")],
  "example-2-ssr-trial-summary.csv"
)

# The manuscript's production-scale results are stored as package data. The
# generators are intentionally not run here because they require many hours.
cache_files <- c(
  "inst/extdata/simulation_results.rds",
  "inst/extdata/null_simulation_results.rds",
  "inst/extdata/gs_simulation_results.rds",
  "inst/extdata/jensen_broad_sweep_results.rds",
  "inst/extdata/jensen_score_broad_sweep_results.rds",
  "inst/extdata/score_sweep_summary.rds",
  "inst/extdata/ssr_sim_vignette_summary_score.rds",
  "inst/extdata/ssr_sizing_sensitivity.rds"
)
cache_manifest <- data.frame(
  file = cache_files,
  available = file.exists(cache_files),
  bytes = ifelse(file.exists(cache_files), file.info(cache_files)$size, NA_real_)
)
write_table(cache_manifest, "production-cache-manifest.csv")

commands <- c(
  "Production-cache regeneration commands (run from the package root):",
  "Rscript data-raw/generate_simulation_data.R",
  "Rscript data-raw/generate_gs_simulation_data.R",
  "Rscript data-raw/jensen_correction_broad_sweep.R",
  "Rscript data-raw/jensen_score_broad_sweep.R",
  "Rscript data-raw/generate_score_sweep.R",
  "GSDESIGNNB_PRODUCTION_SSR=true GSDESIGNNB_FORCE_SSR_SIM=true Rscript data-raw/generate_ssr_score_outputs.R",
  "Rscript data-raw/generate_ssr_sizing_sensitivity.R",
  "",
  "These production simulations are intentionally excluded from the default run."
)
writeLines(commands, file.path(output_dir, "production-cache-commands.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "session-info.txt"))

message("Paper replication outputs written to: ", normalizePath(output_dir))
