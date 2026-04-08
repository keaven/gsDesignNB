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
