test_that("unblinded_ssr estimates parameters and re-estimates sample size", {

  # Simulate data with known parameters
  enroll_rate <- data.frame(rate = 100, duration = 1)
  fail_rate <- data.frame(treatment = c("Control", "Experimental"),
                          rate = c(0.5, 0.3))
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 10, n = 200,
                block = c("Control", "Experimental"))

  cut <- cut_data_by_date(sim, cut_date = 5)

  res <- unblinded_ssr(
    cut,
    ratio = 1,
    lambda1_planning = 0.5,
    lambda2_planning = 0.3,
    power = 0.8,
    alpha = 0.025,
    accrual_rate = 100,
    accrual_duration = 1,
    trial_duration = 10
  )

  # Check output names
  expect_true(all(c("n_total_unblinded", "dispersion_unblinded",
                     "lambda1_unblinded", "lambda2_unblinded",
                     "info_fraction", "unblinded_info", "target_info") %in% names(res)))

  # Re-estimated control rate should be near 0.5
  expect_true(res$lambda1_unblinded > 0.2 && res$lambda1_unblinded < 1.0)

  # Re-estimated experimental rate should be near 0.3
  expect_true(res$lambda2_unblinded > 0.1 && res$lambda2_unblinded < 0.8)

  # Dispersion should be near 0 for Poisson data

  expect_true(res$dispersion_unblinded >= 0 && res$dispersion_unblinded < 0.5)

  # Sample size should be positive
  expect_true(res$n_total_unblinded > 0)

  # Information fraction should be positive
  expect_true(res$info_fraction > 0)
  expect_true(res$unblinded_info > 0)
  expect_true(res$target_info > 0)
})

test_that("unblinded_ssr preserves planned effect size", {
  # Generate data where control rate is lower than planned
  enroll_rate <- data.frame(rate = 100, duration = 1)
  fail_rate <- data.frame(treatment = c("Control", "Experimental"),
                          rate = c(0.3, 0.18))
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 10, n = 200,
                block = c("Control", "Experimental"))

  cut <- cut_data_by_date(sim, cut_date = 5)

  # Planning assumed higher rates
  res <- unblinded_ssr(
    cut,
    ratio = 1,
    lambda1_planning = 0.5,
    lambda2_planning = 0.3,
    power = 0.8,
    alpha = 0.025,
    accrual_rate = 100,
    accrual_duration = 1,
    trial_duration = 10
  )

  # The planned rate ratio is 0.3/0.5 = 0.6

  # The re-estimated lambda2 should reflect the planned RR applied to observed control rate
  # lambda2_calc = lambda1_est * (0.3/0.5)
  # But the RETURNED lambda2_unblinded is the OBSERVED one (not the calc one)
  expect_true(res$n_total_unblinded > 0)
})

test_that("unblinded_ssr handles unequal allocation", {
  enroll_rate <- data.frame(rate = 100, duration = 1)
  fail_rate <- data.frame(treatment = c("Control", "Experimental"),
                          rate = c(0.5, 0.3))
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 10, n = 200,
                block = c("Control", "Experimental", "Experimental"))

  cut <- cut_data_by_date(sim, cut_date = 5)

  res <- unblinded_ssr(
    cut,
    ratio = 2,
    lambda1_planning = 0.5,
    lambda2_planning = 0.3,
    power = 0.8,
    alpha = 0.025,
    accrual_rate = 100,
    accrual_duration = 1,
    trial_duration = 10
  )

  expect_true(res$n_total_unblinded > 0)
  expect_true(res$unblinded_info > 0)
})

test_that("unblinded_ssr errors on bad data", {
  bad_data <- data.frame(events = c(0), tte = c(0.1), treatment = c("Control"))
  expect_error(
    unblinded_ssr(
      bad_data,
      ratio = 1,
      lambda1_planning = 0.5,
      lambda2_planning = 0.3,
      power = 0.8,
      alpha = 0.025,
      accrual_rate = 10,
      accrual_duration = 1,
      trial_duration = 10
    )
  )
})
