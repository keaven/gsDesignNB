test_that("mutze_test runs on cut data", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3))
  dropout_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.1, 0.05), duration = c(10, 10))
  sim <- nb_sim(enroll_rate, fail_rate, dropout_rate, max_followup = 1.5, n = 40)
  cut <- cut_data_by_date(sim, cut_date = 1)
  res <- mutze_test(cut)
  expect_true(is.list(res))
  expect_true(all(c("estimate", "se", "z", "p_value", "rate_ratio") %in% names(res)))
  expect_true(is.finite(res$estimate))
  expect_true(nrow(res$group_summary) == 2)
})

test_that("mutze_test handles poisson option", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3))
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 1, n = 20)
  cut <- cut_data_by_date(sim, cut_date = 0.8)
  res <- mutze_test(cut, method = "poisson")
  expect_true(grepl("Poisson", res$method))
})

test_that("mutze_test errors on missing columns", {
  bad_data <- data.frame(x = 1, y = 2, z = 3)
  expect_error(mutze_test(bad_data), "treatment.*events.*tte")
})

test_that("mutze_test errors with no positive follow-up", {
  df <- data.frame(treatment = c("A", "B"), events = c(0, 0), tte = c(0, 0))
  expect_error(mutze_test(df), "No rows with positive follow-up")
})

test_that("mutze_test errors with non-two groups", {
  df <- data.frame(treatment = c("A", "A", "A"), events = c(1, 2, 3), tte = c(1, 1, 1))
  expect_error(mutze_test(df), "exactly two treatment groups")
})

test_that("mutze_test drops zero-tte rows", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3))
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 1, n = 20)
  cut <- cut_data_by_date(sim, cut_date = 0.8)

  # Add a zero-tte row
  extra <- data.frame(
    id = 999, treatment = "Control", enroll_time = 0,
    tte = 0, tte_total = 0, events = 0L
  )
  cut_with_zero <- rbind(cut, extra)

  res <- mutze_test(cut_with_zero)
  expect_true(is.finite(res$estimate))
})

test_that("mutze_test two-sided p-value", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3))
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 1, n = 40)
  cut <- cut_data_by_date(sim, cut_date = 0.8)
  res <- mutze_test(cut, sided = 2)
  expect_equal(res$sided, 2)
  expect_true(res$p_value >= 0 && res$p_value <= 1)
})

test_that("mutze_test falls back to Poisson on near-Poisson data", {
  # Very high theta (near-Poisson) should trigger fallback
  # nb_sim generates Poisson data (k=0), so theta will be very large
  enroll_rate <- data.frame(rate = 30, duration = 1)
  fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3))
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 2, n = 60)
  cut <- cut_data_by_date(sim, cut_date = 1.5)

  res <- mutze_test(cut, poisson_threshold = 1000)
  # Should get a valid result even if it fell back
  expect_true(is.finite(res$z))
})

test_that("print.mutze_test works", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3))
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 1, n = 20)
  cut <- cut_data_by_date(sim, cut_date = 0.8)
  res <- mutze_test(cut)

  expect_output(print(res), "Mutze Test Results")
  expect_output(print(res), "Rate Ratio")
  expect_output(print(res), "Group Summary")
})
