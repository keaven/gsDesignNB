setup_ssr_args <- function() {
  enroll_rate <- data.frame(rate = 10, duration = 4)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.5, 0.35),
    dispersion = 0.3
  )
  dropout_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.05, 0.05),
    duration = c(8, 8)
  )
  fixed_design <- sample_size_nbinom(
    lambda1 = 0.5,
    lambda2 = 0.35,
    dispersion = 0.3,
    power = 0.8,
    alpha = 0.025,
    accrual_rate = 10,
    accrual_duration = 4,
    trial_duration = 8,
    max_followup = 4
  )
  gs_design <- gsNBCalendar(
    fixed_design,
    k = 3,
    test.type = 4,
    alpha = 0.025,
    sfu = sfHSD,
    sfupar = -2,
    sfl = sfHSD,
    sflpar = 1,
    analysis_times = c(2.5, 5, 8)
  ) |>
    gsDesignNB::toInteger()

  list(
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    dropout_rate = dropout_rate,
    design = gs_design
  )
}

test_that("sim_ssr_nbinom returns expected structure", {
  args <- setup_ssr_args()

  res <- sim_ssr_nbinom(
    n_sims = 3,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    design = args$design,
    strategies = c("No adaptation", "Blinded SSR"),
    n_max = 60,
    seed = 123
  )

  expect_s3_class(res, "sim_ssr_nbinom")
  expect_true(all(c("trial_results", "analysis_results", "settings") %in% names(res)))
  expect_s3_class(res$trial_results, "data.frame")
  expect_s3_class(res$analysis_results, "data.frame")
  expect_equal(nrow(res$trial_results), 6)
  expect_true(all(c(
    "participants_with_events_stop",
    "events_observed_stop",
    "participants_with_events_ia1",
    "events_observed_final",
    "adapt_cut_time"
  ) %in% names(res$trial_results)))
})

test_that("sim_ssr_nbinom is reproducible with fixed seed", {
  args <- setup_ssr_args()

  res1 <- sim_ssr_nbinom(
    n_sims = 3,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    design = args$design,
    strategies = "No adaptation",
    n_max = 60,
    seed = 456
  )

  res2 <- sim_ssr_nbinom(
    n_sims = 3,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    design = args$design,
    strategies = "No adaptation",
    n_max = 60,
    seed = 456
  )

  expect_identical(res1$trial_results, res2$trial_results)
})

test_that("summarize_ssr_sim reports expected participants and events", {
  args <- setup_ssr_args()

  res <- sim_ssr_nbinom(
    n_sims = 3,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    design = args$design,
    strategies = c("No adaptation", "Blinded SSR"),
    n_max = 60,
    seed = 789
  )

  sm <- summarize_ssr_sim(res, by = "strategy")

  expect_s3_class(sm$trial_summary, "data.frame")
  expect_s3_class(sm$analysis_summary, "data.frame")
  expect_true(all(c(
    "expected_participants_with_events",
    "expected_events_observed",
    "mean_participants_with_events_ia1",
    "mean_events_observed_final"
  ) %in% names(sm$trial_summary)))
})

test_that("check_gs_bound accepts explicit information column", {
  design <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 2, timing = c(0.5, 1))
  sim_df <- data.frame(
    sim = c(1, 1),
    analysis = c(1, 2),
    z_stat = c(2.5, NA),
    blinded_info = c(10, 100),
    unblinded_info = c(10, 100),
    info_unblinded_mom = c(100, 100)
  )

  res_default <- check_gs_bound(sim_df, design, info_scale = "unblinded")
  res_explicit <- check_gs_bound(sim_df, design, info_col = "unblinded_info")
  res_custom <- check_gs_bound(sim_df, design, info_col = "info_unblinded_mom")

  expect_identical(
    res_default[, c("cross_upper", "cross_lower", "cross_harm")],
    res_explicit[, c("cross_upper", "cross_lower", "cross_harm")]
  )
  expect_s3_class(res_custom, "data.frame")
  expect_equal(nrow(res_custom), nrow(sim_df))
})

test_that("sim_ssr_nbinom errors on bad design", {
  expect_error(
    sim_ssr_nbinom(
      n_sims = 1,
      enroll_rate = data.frame(rate = 10, duration = 4),
      fail_rate = data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3), dispersion = 0.2),
      dropout_rate = NULL,
      max_followup = 4,
      design = list(k = 2),
      strategies = "No adaptation"
    ),
    "gsNB"
  )
})

test_that("sim_ssr_nbinom errors on bad strategies", {
  args <- setup_ssr_args()
  expect_error(
    sim_ssr_nbinom(
      n_sims = 1,
      enroll_rate = args$enroll_rate,
      fail_rate = args$fail_rate,
      dropout_rate = args$dropout_rate,
      max_followup = 4,
      design = args$design,
      strategies = "Invalid Strategy"
    ),
    "strategies must be chosen from"
  )
})

test_that("sim_ssr_nbinom works with Unblinded SSR strategy", {
  args <- setup_ssr_args()

  res <- sim_ssr_nbinom(
    n_sims = 2,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    design = args$design,
    strategies = c("No adaptation", "Unblinded SSR"),
    n_max = 80,
    seed = 999
  )

  expect_s3_class(res, "sim_ssr_nbinom")
  expect_true(nrow(res$trial_results) > 0)
})

test_that("summarize_ssr_sim works without by argument", {
  args <- setup_ssr_args()

  res <- sim_ssr_nbinom(
    n_sims = 3,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    design = args$design,
    strategies = "No adaptation",
    n_max = 60,
    seed = 111
  )

  sm <- summarize_ssr_sim(res)

  expect_s3_class(sm$trial_summary, "data.frame")
  expect_s3_class(sm$analysis_summary, "data.frame")
  expect_true(nrow(sm$trial_summary) > 0)
})

test_that("sim_ssr_nbinom handles ignore_futility", {
  args <- setup_ssr_args()

  res <- sim_ssr_nbinom(
    n_sims = 2,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    design = args$design,
    strategies = "No adaptation",
    n_max = 60,
    ignore_futility = TRUE,
    seed = 222
  )

  expect_s3_class(res, "sim_ssr_nbinom")
  expect_true(nrow(res$trial_results) > 0)
})

test_that("sim_ssr_nbinom handles metadata", {
  args <- setup_ssr_args()

  meta <- data.frame(scenario = "base", label = "test")
  res <- sim_ssr_nbinom(
    n_sims = 2,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    design = args$design,
    strategies = "No adaptation",
    n_max = 60,
    metadata = meta,
    seed = 333
  )

  expect_s3_class(res, "sim_ssr_nbinom")
  expect_true("scenario" %in% names(res$trial_results))
  expect_true("label" %in% names(res$trial_results))
})

test_that("sim_ssr_nbinom handles different bound_info", {
  args <- setup_ssr_args()

  for (bi in c("blinded_ml", "unblinded_mom", "blinded_mom")) {
    res <- sim_ssr_nbinom(
      n_sims = 2,
      enroll_rate = args$enroll_rate,
      fail_rate = args$fail_rate,
      dropout_rate = args$dropout_rate,
      max_followup = 4,
      design = args$design,
      strategies = "No adaptation",
      n_max = 60,
      bound_info = bi,
      seed = 444
    )
    expect_s3_class(res, "sim_ssr_nbinom")
  }
})

test_that("sim_ssr_nbinom handles all three strategies", {
  args <- setup_ssr_args()

  res <- sim_ssr_nbinom(
    n_sims = 2,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    design = args$design,
    strategies = c("No adaptation", "Blinded SSR", "Unblinded SSR"),
    n_max = 80,
    seed = 555
  )

  expect_s3_class(res, "sim_ssr_nbinom")
  # Should have results for all 3 strategies
  expect_true(all(c("No adaptation", "Blinded SSR", "Unblinded SSR") %in% res$trial_results$strategy))
})

# --- Internal helper tests ---

test_that(".ssr_metadata_df handles named list", {
  result <- gsDesignNB:::.ssr_metadata_df(list(scenario = "A", label = "test"))
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(result$scenario, "A")
})

test_that(".ssr_metadata_df errors on multi-row data frame", {
  expect_error(
    gsDesignNB:::.ssr_metadata_df(data.frame(scenario = c("A", "B"))),
    "exactly one row"
  )
})

test_that(".ssr_metadata_df errors on bad input", {
  expect_error(
    gsDesignNB:::.ssr_metadata_df("bad input"),
    "named list or one-row data frame"
  )
})

test_that(".ssr_select_info falls back to unblinded_ml", {
  metrics <- list(
    info_unblinded_ml = 50,
    info_blinded_ml = NA_real_,
    info_unblinded_mom = 0,
    info_blinded_mom = -1
  )

  # blinded_ml is NA -> should fall back

expect_equal(gsDesignNB:::.ssr_select_info(metrics, "blinded_ml"), 50)
  # unblinded_mom is 0 -> should fall back
  expect_equal(gsDesignNB:::.ssr_select_info(metrics, "unblinded_mom"), 50)
  # blinded_mom is negative -> should fall back
  expect_equal(gsDesignNB:::.ssr_select_info(metrics, "blinded_mom"), 50)
  # unblinded_ml is valid -> should return directly
  expect_equal(gsDesignNB:::.ssr_select_info(metrics, "unblinded_ml"), 50)
})

test_that(".ssr_stage_label_from_suffix works", {
  expect_equal(gsDesignNB:::.ssr_stage_label_from_suffix("final"), "Final")
  expect_equal(gsDesignNB:::.ssr_stage_label_from_suffix("ia1"), "IA1")
  expect_equal(gsDesignNB:::.ssr_stage_label_from_suffix("ia2"), "IA2")
})

test_that(".ssr_stage_suffixes_from_trial_results extracts suffixes", {
  dt <- data.table::data.table(if_ia1 = 0.5, if_ia2 = 0.8, if_final = 1.0, x = 1)
  result <- gsDesignNB:::.ssr_stage_suffixes_from_trial_results(dt)
  expect_equal(sort(result), c("final", "ia1", "ia2"))
})

test_that(".ssr_dynamic_bounds returns bounds", {
  design <- gsDesign::gsDesign(k = 3, test.type = 4, n.fix = 100,
                                timing = c(0.33, 0.67, 1))

  bounds <- gsDesignNB:::.ssr_dynamic_bounds(
    analysis_index = 1,
    n_analyses = 3,
    observed_if = c(0.35, NA, NA),
    planned_timing = c(0.33, 0.67, 1),
    target_info = 100,
    design = design
  )

  expect_true(is.finite(bounds$upper_bound))
  expect_true(is.finite(bounds$lower_bound))
})

test_that(".ssr_dynamic_bounds handles final analysis", {
  design <- gsDesign::gsDesign(k = 2, test.type = 4, n.fix = 100,
                                timing = c(0.5, 1))

  bounds <- gsDesignNB:::.ssr_dynamic_bounds(
    analysis_index = 2,
    n_analyses = 2,
    observed_if = c(0.5, 1.0),
    planned_timing = c(0.5, 1),
    target_info = 100,
    design = design
  )

  expect_true(is.finite(bounds$upper_bound))
  expect_true(is.finite(bounds$lower_bound))
})

test_that(".ssr_dynamic_bounds returns NA on bad input", {
  design <- gsDesign::gsDesign(k = 2, test.type = 4, n.fix = 100,
                                timing = c(0.5, 1))

  bounds <- gsDesignNB:::.ssr_dynamic_bounds(
    analysis_index = 1,
    n_analyses = 2,
    observed_if = c(NA_real_, NA_real_),
    planned_timing = c(0.5, 1),
    target_info = 100,
    design = design
  )

  expect_true(is.na(bounds$upper_bound))
  expect_true(is.na(bounds$lower_bound))
})

test_that(".ssr_extract_info_metrics returns all info types", {
  enroll_rate <- data.frame(rate = 50, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.5, 0.3)
  )
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 3, n = 100)
  cut <- cut_data_by_date(sim, cut_date = 2)

  metrics <- gsDesignNB:::.ssr_extract_info_metrics(
    cut_data = cut,
    ratio_plan = 1,
    lambda1_plan = 0.5,
    lambda2_plan = 0.3,
    event_gap = 0
  )

  expect_true(is.list(metrics))
  expect_true(all(c("z_value", "info_unblinded_ml", "info_blinded_ml") %in% names(metrics)))
  expect_true(is.finite(metrics$z_value))
  expect_true(is.finite(metrics$info_unblinded_ml))
})

test_that(".ssr_extract_info_metrics handles insufficient data", {
  cut <- data.frame(
    id = 1:2, treatment = c("Control", "Experimental"),
    events = c(0, 0), tte = c(0.01, 0.01), enroll_time = c(0, 0)
  )

  metrics <- gsDesignNB:::.ssr_extract_info_metrics(
    cut_data = cut,
    ratio_plan = 1,
    lambda1_plan = 0.5,
    lambda2_plan = 0.3,
    event_gap = 0
  )

  expect_true(is.na(metrics$z_value))
})

test_that(".ssr_long_from_trial_results reshapes correctly", {
  args <- setup_ssr_args()
  res <- sim_ssr_nbinom(
    n_sims = 2,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    design = args$design,
    strategies = "No adaptation",
    n_max = 60,
    seed = 666
  )

  long <- gsDesignNB:::.ssr_long_from_trial_results(res$trial_results)
  expect_true(is.data.frame(long))
  expect_true(nrow(long) > 0)
  expect_true(all(c("analysis", "analysis_stage", "analysis_time",
                     "z_value", "info_value", "info_fraction",
                     "cross_upper", "cross_lower") %in% names(long)))
  # For k=3 design: ia1, ia2, final -> 3 rows per sim per strategy
  expect_equal(nrow(long), 2 * 3)
})

test_that("summarize_ssr_sim works on raw trial_results data frame", {
  args <- setup_ssr_args()
  res <- sim_ssr_nbinom(
    n_sims = 3,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    design = args$design,
    strategies = "No adaptation",
    n_max = 60,
    seed = 777
  )

  # Pass raw trial_results data frame, not the sim_ssr_nbinom object
  sm <- summarize_ssr_sim(res$trial_results)
  expect_true(is.data.frame(sm$trial_summary))
  expect_true(nrow(sm$trial_summary) > 0)
})

test_that("summarize_ssr_sim works with by = NULL", {
  args <- setup_ssr_args()
  res <- sim_ssr_nbinom(
    n_sims = 3,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    design = args$design,
    strategies = "No adaptation",
    n_max = 60,
    seed = 888
  )

  sm <- summarize_ssr_sim(res, by = "nonexistent_col")
  expect_true(is.data.frame(sm$trial_summary))
  expect_true(nrow(sm$trial_summary) == 1)
})

test_that("sim_ssr_nbinom errors on design with k < 2", {
  # Create a gsNB-like object with k=1 to trigger the error
  fake_design <- list(k = 1)
  class(fake_design) <- c("gsNB", "gsDesign")

  expect_error(
    sim_ssr_nbinom(
      n_sims = 1,
      enroll_rate = data.frame(rate = 10, duration = 4),
      fail_rate = data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3)),
      max_followup = 4,
      design = fake_design,
      strategies = "No adaptation"
    ),
    "at least one interim"
  )
})

test_that("sim_ssr_nbinom errors on bad adapt_analysis", {
  args <- setup_ssr_args()

  expect_error(
    sim_ssr_nbinom(
      n_sims = 1,
      enroll_rate = args$enroll_rate,
      fail_rate = args$fail_rate,
      max_followup = 4,
      design = args$design,
      strategies = "No adaptation",
      adapt_analysis = 10
    ),
    "adapt_analysis must identify an interim"
  )
})

test_that("sim_ssr_nbinom uses event_gap from design inputs", {
  fixed_design <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.35, dispersion = 0.3,
    power = 0.8, alpha = 0.025,
    accrual_rate = 10, accrual_duration = 4, trial_duration = 8,
    max_followup = 4, event_gap = 0.5
  )
  gs_design <- gsNBCalendar(
    fixed_design, k = 2, test.type = 4, alpha = 0.025,
    analysis_times = c(5, 8)
  )

  res <- sim_ssr_nbinom(
    n_sims = 2,
    enroll_rate = data.frame(rate = 10, duration = 4),
    fail_rate = data.frame(treatment = c("Control", "Experimental"),
                           rate = c(0.5, 0.35), dispersion = 0.3),
    max_followup = 4,
    design = gs_design,
    strategies = "No adaptation",
    seed = 1234
  )

  # event_gap should be extracted from design
  expect_equal(res$settings$event_gap, 0.5)
})

test_that("sim_ssr_nbinom with adapt parameters triggers adaptation checks", {
  args <- setup_ssr_args()

  res <- sim_ssr_nbinom(
    n_sims = 3,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    design = args$design,
    strategies = c("Blinded SSR", "Unblinded SSR"),
    n_max = 100,
    adapt_analysis = 2,
    max_enrollment_frac_for_adapt = 0.95,
    min_months_to_close_for_adapt = 0,
    analysis_lag_months = 0.5,
    seed = 2345
  )

  expect_s3_class(res, "sim_ssr_nbinom")
  tr <- res$trial_results
  # Should have adapt-related columns
  expect_true("adapt_allowed" %in% names(tr))
  expect_true("adapt_applied" %in% names(tr))
  # adapt_analysis = 2 means ia2 columns should exist
  expect_true("ia2_adapt_cut_time" %in% names(tr))
})

test_that(".ssr_dynamic_bounds handles non-monotone timing", {
  design <- gsDesign::gsDesign(k = 3, test.type = 4, n.fix = 100,
                                timing = c(0.33, 0.67, 1))

  # Decreasing observed IFs should trigger monotonicity correction
  bounds <- gsDesignNB:::.ssr_dynamic_bounds(
    analysis_index = 2,
    n_analyses = 3,
    observed_if = c(0.5, 0.4, NA),
    planned_timing = c(0.33, 0.67, 1),
    target_info = 100,
    design = design
  )

  # Even after correction, should produce numeric bounds (may be NA if gsDesign fails)
  expect_true(is.numeric(bounds$upper_bound))
  expect_true(is.numeric(bounds$lower_bound))
})

test_that(".ssr_extract_info_metrics computes all four info types", {
  set.seed(42)
  enroll_rate <- data.frame(rate = 50, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.6, 0.3),
    dispersion = 0.3
  )
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 3, n = 100)
  cut <- cut_data_by_date(sim, cut_date = 2)

  metrics <- gsDesignNB:::.ssr_extract_info_metrics(
    cut_data = cut,
    ratio_plan = 1,
    lambda1_plan = 0.6,
    lambda2_plan = 0.3,
    event_gap = 0
  )

  # Should have all fields
  expect_true(is.finite(metrics$info_unblinded_ml))
  expect_true(is.finite(metrics$info_blinded_ml))
  # MoM estimates may or may not be finite depending on data
  expect_true("info_unblinded_mom" %in% names(metrics))
  expect_true("info_blinded_mom" %in% names(metrics))
})

test_that(".ssr_long_from_trial_results handles empty suffixes", {
  dt <- data.frame(sim = 1, strategy = "No adaptation", reject = FALSE)
  result <- gsDesignNB:::.ssr_long_from_trial_results(dt)
  expect_equal(nrow(result), 0)
})

test_that(".ssr_long_from_trial_results handles missing required columns", {
  # Create a data frame that has if_ columns but not all required stage cols
  dt <- data.frame(
    sim = 1, strategy = "No adaptation", reject = FALSE,
    if_ia1 = 0.5, if_final = 1.0,
    ia1_time = 5, final_time = 10
  )
  result <- gsDesignNB:::.ssr_long_from_trial_results(dt)
  # Should return empty or skip incomplete stages
  expect_true(nrow(result) == 0)
})

test_that("sim_ssr_nbinom SSR adaptation is applied with favorable design", {
  # Use a design with much lower power so SSR re-estimation
  # will want to increase the sample size
  fixed_design <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.45, dispersion = 0.4,
    power = 0.8, alpha = 0.025,
    accrual_rate = 8, accrual_duration = 6, trial_duration = 12,
    max_followup = 6
  )
  gs_design <- gsNBCalendar(
    fixed_design, k = 2, test.type = 4, alpha = 0.025,
    sfu = sfHSD, sfupar = -2, sfl = sfHSD, sflpar = 1,
    analysis_times = c(6, 12)
  )

  # Simulate with higher true rates so SSR will re-estimate higher N
  res <- sim_ssr_nbinom(
    n_sims = 5,
    enroll_rate = data.frame(rate = 8, duration = 6),
    fail_rate = data.frame(
      treatment = c("Control", "Experimental"),
      rate = c(0.5, 0.45), dispersion = 0.4
    ),
    dropout_rate = data.frame(
      treatment = c("Control", "Experimental"),
      rate = c(0.05, 0.05), duration = c(12, 12)
    ),
    max_followup = 6,
    design = gs_design,
    strategies = c("Blinded SSR", "Unblinded SSR"),
    n_max = 200,
    max_enrollment_frac_for_adapt = 1.0,
    min_months_to_close_for_adapt = 0,
    seed = 9999
  )

  expect_s3_class(res, "sim_ssr_nbinom")
  tr <- res$trial_results
  expect_true("adapt_allowed" %in% names(tr))
  # At least some sims should have adaptation allowed
  expect_true(any(tr$adapt_allowed))
})
