test_that("sim_gs_nbinom runs basic simulation", {
  set.seed(123)
  enroll_rate <- data.frame(rate = 10, duration = 3)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.6, 0.4),
    dispersion = 0.2
  )
  dropout_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.05, 0.05),
    duration = c(6, 6)
  )
  design <- sample_size_nbinom(
    lambda1 = 0.6, lambda2 = 0.4, dispersion = 0.2, power = 0.8,
    accrual_rate = enroll_rate$rate, accrual_duration = enroll_rate$duration,
    trial_duration = 6
  )
  cuts <- list(
    list(planned_calendar = 2),
    list(planned_calendar = 4)
  )
  sim_results <- sim_gs_nbinom(
    n_sims = 2,
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    dropout_rate = dropout_rate,
    max_followup = 4,
    n_target = 30,
    design = design,
    cuts = cuts,
    seed = TRUE
  )

  expect_true(is.data.frame(sim_results))

  # Expected columns
  expected_cols <- c("sim", "analysis", "analysis_time", "n_enrolled",
                     "n_ctrl", "n_exp", "events_total", "z_stat",
                     "blinded_info", "unblinded_info",
                     "info_unblinded_ml", "info_blinded_ml",
                     "info_unblinded_mom", "info_blinded_mom")
  expect_true(all(expected_cols %in% names(sim_results)))

  # 2 sims x 2 analyses = 4 rows
  expect_equal(nrow(sim_results), 4)

  # Sim IDs should be 1 and 2
  expect_equal(sort(unique(sim_results$sim)), c(1, 2))

  # Analysis IDs should be 1 and 2
  expect_equal(sort(unique(sim_results$analysis)), c(1, 2))

  # Analysis times should be non-decreasing within each sim
  for (s in unique(sim_results$sim)) {
    sub <- sim_results[sim_results$sim == s, ]
    expect_true(all(diff(sub$analysis_time) >= 0))
  }
})

test_that("sim_gs_nbinom errors without design", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.5, 0.3),
    dispersion = 0.1
  )

  expect_error(
    sim_gs_nbinom(
      n_sims = 1,
      enroll_rate = enroll_rate,
      fail_rate = fail_rate,
      max_followup = 2,
      cuts = list(list(planned_calendar = 1))
    ),
    "design object must be provided"
  )
})

test_that("sim_gs_nbinom errors without cuts or analysis_times", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.5, 0.3),
    dispersion = 0.1
  )
  design <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 1, trial_duration = 3
  )

  expect_error(
    sim_gs_nbinom(
      n_sims = 1,
      enroll_rate = enroll_rate,
      fail_rate = fail_rate,
      max_followup = 2,
      design = design
    ),
    "Either analysis_times or cuts must be provided"
  )
})

test_that("sim_gs_nbinom accepts analysis_times", {
  set.seed(42)
  enroll_rate <- data.frame(rate = 10, duration = 3)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.6, 0.4),
    dispersion = 0.2
  )
  design <- sample_size_nbinom(
    lambda1 = 0.6, lambda2 = 0.4, dispersion = 0.2, power = 0.8,
    accrual_rate = 10, accrual_duration = 3, trial_duration = 6
  )

  # Use analysis_times instead of cuts
  sim_results <- sim_gs_nbinom(
    n_sims = 2,
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    max_followup = 4,
    n_target = 30,
    design = design,
    analysis_times = c(2, 4),
    seed = TRUE
  )

  expect_true(is.data.frame(sim_results))
  expect_equal(nrow(sim_results), 4)
})

test_that("sim_gs_nbinom is reproducible with seed", {
  enroll_rate <- data.frame(rate = 10, duration = 2)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.6, 0.4),
    dispersion = 0.2
  )
  design <- sample_size_nbinom(
    lambda1 = 0.6, lambda2 = 0.4, dispersion = 0.2, power = 0.8,
    accrual_rate = 10, accrual_duration = 2, trial_duration = 4
  )
  cuts <- list(list(planned_calendar = 2), list(planned_calendar = 4))

  set.seed(42)
  res1 <- sim_gs_nbinom(
    n_sims = 3, enroll_rate = enroll_rate, fail_rate = fail_rate,
    max_followup = 3, n_target = 20, design = design, cuts = cuts, seed = TRUE
  )

  set.seed(42)
  res2 <- sim_gs_nbinom(
    n_sims = 3, enroll_rate = enroll_rate, fail_rate = fail_rate,
    max_followup = 3, n_target = 20, design = design, cuts = cuts, seed = TRUE
  )

  expect_identical(res1, res2)
})
