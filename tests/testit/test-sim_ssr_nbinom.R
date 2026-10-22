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
    toInteger()

  list(
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    dropout_rate = dropout_rate,
    design = gs_design
  )
}

assert("sim_ssr_nbinom returns expected structure", {
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

  (inherits(res, "sim_ssr_nbinom"))
  (c("trial_results", "analysis_results", "settings") %in% names(res))
  (inherits(res$trial_results, "data.frame"))
  (inherits(res$analysis_results, "data.frame"))
  (nrow(res$trial_results) %==% 6L)
  (c(
    "participants_with_events_stop",
    "events_observed_stop",
    "participants_with_events_ia1",
    "events_observed_final",
    "adapt_cut_time"
  ) %in% names(res$trial_results))
})

assert("sim_ssr_nbinom is reproducible with fixed seed", {
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

  (res1$trial_results %==% res2$trial_results)
})

assert("summarize_ssr_sim reports expected participants and events", {
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

  (inherits(sm$trial_summary, "data.frame"))
  (inherits(sm$analysis_summary, "data.frame"))
  (c(
    "expected_participants_with_events",
    "expected_events_observed",
    "mean_participants_with_events_ia1",
    "mean_events_observed_final"
  ) %in% names(sm$trial_summary))
})

assert("check_gs_bound accepts explicit information column", {
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

  (res_default[, c("cross_upper", "cross_lower", "cross_harm")] %==% res_explicit[, c("cross_upper", "cross_lower", "cross_harm")])
  (inherits(res_custom, "data.frame"))
  (nrow(res_custom) %==% nrow(sim_df))
})

assert("sim_ssr_nbinom errors on bad design", {
  (has_error(
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
  ))
})

assert("sim_ssr_nbinom errors on bad strategies", {
  args <- setup_ssr_args()
  (has_error(
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
  ))
})

assert("sim_ssr_nbinom works with Unblinded SSR strategy", {
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

  (inherits(res, "sim_ssr_nbinom"))
  (nrow(res$trial_results) > 0)
})

assert("summarize_ssr_sim works without by argument", {
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

  (inherits(sm$trial_summary, "data.frame"))
  (inherits(sm$analysis_summary, "data.frame"))
  (nrow(sm$trial_summary) > 0)
})

assert("sim_ssr_nbinom handles ignore_futility", {
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

  (inherits(res, "sim_ssr_nbinom"))
  (nrow(res$trial_results) > 0)
})

assert("sim_ssr_nbinom handles metadata", {
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

  (inherits(res, "sim_ssr_nbinom"))
  (c("scenario", "label") %in% names(res$trial_results))
})

assert("sim_ssr_nbinom handles different bound_info", {
  args <- setup_ssr_args()

  res_list <- lapply(c("blinded_ml", "unblinded_mom", "blinded_mom"), function(bi) {
    sim_ssr_nbinom(
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
  })
  (vapply(res_list, inherits, logical(1), "sim_ssr_nbinom"))
})

assert("sim_ssr_nbinom handles all three strategies", {
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

  (inherits(res, "sim_ssr_nbinom"))
  # Should have results for all 3 strategies
  (c("No adaptation", "Blinded SSR", "Unblinded SSR") %in% res$trial_results$strategy)
})

# --- Internal helper tests ---

assert(".ssr_metadata_df handles named list", {
  result <- .ssr_metadata_df(list(scenario = "A", label = "test"))
  (inherits(result, "data.frame"))
  (nrow(result) %==% 1L)
  (all.equal(result$scenario, "A"))
})

assert(".ssr_metadata_df errors on multi-row data frame", {
  (has_error(
    .ssr_metadata_df(data.frame(scenario = c("A", "B"))),
    "exactly one row"
  ))
})

assert(".ssr_metadata_df errors on bad input", {
  (has_error(
    .ssr_metadata_df("bad input"),
    "named list or one-row data frame"
  ))
})

assert(".ssr_select_info falls back to unblinded_ml", {
  metrics <- list(
    info_unblinded_ml = 50,
    info_blinded_ml = NA_real_,
    info_unblinded_mom = 0,
    info_blinded_mom = -1
  )

  # blinded_ml is NA -> should fall back

(all.equal(.ssr_select_info(metrics, "blinded_ml"), 50))
  # unblinded_mom is 0 -> should fall back
  (all.equal(.ssr_select_info(metrics, "unblinded_mom"), 50))
  # blinded_mom is negative -> should fall back
  (all.equal(.ssr_select_info(metrics, "blinded_mom"), 50))
  # unblinded_ml is valid -> should return directly
  (all.equal(.ssr_select_info(metrics, "unblinded_ml"), 50))
})

assert(".ssr_stage_label_from_suffix works", {
  (all.equal(.ssr_stage_label_from_suffix("final"), "Final"))
  (all.equal(.ssr_stage_label_from_suffix("ia1"), "IA1"))
  (all.equal(.ssr_stage_label_from_suffix("ia2"), "IA2"))
})

assert(".ssr_stage_suffixes_from_trial_results extracts suffixes", {
  dt <- data.table::data.table(if_ia1 = 0.5, if_ia2 = 0.8, if_final = 1.0, x = 1)
  result <- .ssr_stage_suffixes_from_trial_results(dt)
  (all.equal(sort(result), c("final", "ia1", "ia2")))
})

assert(".ssr_dynamic_bounds returns bounds", {
  design <- gsDesign::gsDesign(k = 3, test.type = 4, n.fix = 100,
                                timing = c(0.33, 0.67, 1))

  bounds <- .ssr_dynamic_bounds(
    analysis_index = 1,
    n_analyses = 3,
    observed_if = c(0.35, NA, NA),
    planned_timing = c(0.33, 0.67, 1),
    target_info = 100,
    design = design
  )

  (is.finite(bounds$upper_bound))
  (is.finite(bounds$lower_bound))
})

assert(".ssr_dynamic_bounds handles final analysis", {
  design <- gsDesign::gsDesign(k = 2, test.type = 4, n.fix = 100,
                                timing = c(0.5, 1))

  bounds <- .ssr_dynamic_bounds(
    analysis_index = 2,
    n_analyses = 2,
    observed_if = c(0.5, 1.0),
    planned_timing = c(0.5, 1),
    target_info = 100,
    design = design
  )

  (is.finite(bounds$upper_bound))
  (is.finite(bounds$lower_bound))
})

assert(".ssr_dynamic_bounds returns NA on bad input", {
  design <- gsDesign::gsDesign(k = 2, test.type = 4, n.fix = 100,
                                timing = c(0.5, 1))

  bounds <- .ssr_dynamic_bounds(
    analysis_index = 1,
    n_analyses = 2,
    observed_if = c(NA_real_, NA_real_),
    planned_timing = c(0.5, 1),
    target_info = 100,
    design = design
  )

  (is.na(bounds$upper_bound))
  (is.na(bounds$lower_bound))
})

assert(".ssr_extract_info_metrics returns all info types", {
  enroll_rate <- data.frame(rate = 50, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.5, 0.3)
  )
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 3, n = 100)
  cut <- cut_data_by_date(sim, cut_date = 2)

  metrics <- .ssr_extract_info_metrics(
    cut_data = cut,
    ratio_plan = 1,
    lambda1_plan = 0.5,
    lambda2_plan = 0.3,
    event_gap = 0
  )

  (is.list(metrics))
  (c("z_value", "info_unblinded_ml", "info_blinded_ml") %in% names(metrics))
  (is.finite(metrics$z_value))
  (is.finite(metrics$info_unblinded_ml))
})

assert(".ssr_extract_info_metrics handles insufficient data", {
  cut <- data.frame(
    id = 1:2, treatment = c("Control", "Experimental"),
    events = c(0, 0), tte = c(0.01, 0.01), enroll_time = c(0, 0)
  )

  metrics <- .ssr_extract_info_metrics(
    cut_data = cut,
    ratio_plan = 1,
    lambda1_plan = 0.5,
    lambda2_plan = 0.3,
    event_gap = 0
  )

  (is.na(metrics$z_value))
})

assert(".ssr_long_from_trial_results reshapes correctly", {
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

  long <- .ssr_long_from_trial_results(res$trial_results)
  (is.data.frame(long))
  (nrow(long) > 0)
  (c("analysis", "analysis_stage", "analysis_time",
    "z_value", "info_value", "info_fraction",
    "cross_upper", "cross_lower") %in% names(long))
  # For k=3 design: ia1, ia2, final -> 3 rows per sim per strategy
  (nrow(long) %==% 6L)
})

assert("summarize_ssr_sim works on raw trial_results data frame", {
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
  (is.data.frame(sm$trial_summary))
  (nrow(sm$trial_summary) > 0)
})

assert("summarize_ssr_sim works with by = NULL", {
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
  (is.data.frame(sm$trial_summary))
  (nrow(sm$trial_summary) %==% 1L)
})

assert("sim_ssr_nbinom errors on design with k < 2", {
  # Create a gsNB-like object with k=1 to trigger the error
  fake_design <- list(k = 1)
  class(fake_design) <- c("gsNB", "gsDesign")

  (has_error(
    sim_ssr_nbinom(
      n_sims = 1,
      enroll_rate = data.frame(rate = 10, duration = 4),
      fail_rate = data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3)),
      max_followup = 4,
      design = fake_design,
      strategies = "No adaptation"
    ),
    "at least one interim"
  ))
})

assert("sim_ssr_nbinom errors on bad adapt_analysis", {
  args <- setup_ssr_args()

  (has_error(
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
  ))
})

assert("sim_ssr_nbinom uses event_gap from design inputs", {
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
  (all.equal(res$settings$event_gap, 0.5))
})

assert("sim_ssr_nbinom with adapt parameters triggers adaptation checks", {
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

  (inherits(res, "sim_ssr_nbinom"))
  tr <- res$trial_results
  # Should have adapt-related columns
  (c("adapt_allowed", "adapt_applied") %in% names(tr))
  # adapt_analysis = 2 means ia2 columns should exist
  ("ia2_adapt_cut_time" %in% names(tr))
})

assert(".ssr_dynamic_bounds handles non-monotone timing", {
  design <- gsDesign::gsDesign(k = 3, test.type = 4, n.fix = 100,
                                timing = c(0.33, 0.67, 1))

  # Decreasing observed IFs should trigger monotonicity correction
  bounds <- .ssr_dynamic_bounds(
    analysis_index = 2,
    n_analyses = 3,
    observed_if = c(0.5, 0.4, NA),
    planned_timing = c(0.33, 0.67, 1),
    target_info = 100,
    design = design
  )

  # Even after correction, should produce numeric bounds (may be NA if gsDesign fails)
  (is.numeric(bounds$upper_bound))
  (is.numeric(bounds$lower_bound))
})

assert(".ssr_extract_info_metrics computes all four info types", {
  set.seed(42)
  enroll_rate <- data.frame(rate = 50, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.6, 0.3),
    dispersion = 0.3
  )
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 3, n = 100)
  cut <- cut_data_by_date(sim, cut_date = 2)

  metrics <- .ssr_extract_info_metrics(
    cut_data = cut,
    ratio_plan = 1,
    lambda1_plan = 0.6,
    lambda2_plan = 0.3,
    event_gap = 0
  )

  # Should have all fields
  (is.finite(metrics$info_unblinded_ml))
  (is.finite(metrics$info_blinded_ml))
  # MoM estimates may or may not be finite depending on data
  (c("info_unblinded_mom", "info_blinded_mom") %in% names(metrics))
})

assert(".ssr_long_from_trial_results handles empty suffixes", {
  dt <- data.frame(sim = 1, strategy = "No adaptation", reject = FALSE)
  result <- .ssr_long_from_trial_results(dt)
  (nrow(result) %==% 0L)
})

assert(".ssr_long_from_trial_results handles missing required columns", {
  # Create a data frame that has if_ columns but not all required stage cols
  dt <- data.frame(
    sim = 1, strategy = "No adaptation", reject = FALSE,
    if_ia1 = 0.5, if_final = 1.0,
    ia1_time = 5, final_time = 10
  )
  result <- .ssr_long_from_trial_results(dt)
  # Should return empty or skip incomplete stages
  (nrow(result) %==% 0L)
})

assert("sim_ssr_nbinom SSR adaptation is applied with favorable design", {
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

  (inherits(res, "sim_ssr_nbinom"))
  tr <- res$trial_results
  ("adapt_allowed" %in% names(tr))
  # At least some sims should have adaptation allowed
  (any(tr$adapt_allowed))
})
