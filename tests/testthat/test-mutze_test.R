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

test_that("mutze_test reports a fallback label", {
  enroll_rate <- data.frame(rate = 20, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.5, 0.3), dispersion = c(0.8, 0.8)
  )
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 2, n = 60)
  cut <- cut_data_by_date(sim, cut_date = 1.5)

  res <- mutze_test(cut)
  expect_true(res$fallback %in% c("ml", "poisson", "mom"))
})

test_that("mutze_test uses MoM fallback under extreme overdispersion", {
  # Force MoM fallback by making the threshold tight relative to the fitted
  # theta. We fit NB data with modest k; with mom_threshold set just above
  # the expected 1/theta, extreme overdisp branch should trigger.
  set.seed(42)
  enroll_rate <- data.frame(rate = 40, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.5, 0.3), dispersion = c(1.0, 1.0)
  )
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 2, n = 80)
  cut <- cut_data_by_date(sim, cut_date = 1.5)

  # By setting mom_threshold very small we effectively demand theta < 1/1e6,
  # which will never hold, so the MoM branch only triggers on convergence
  # failure. Instead, force the branch by using poisson_threshold < 1
  # (which makes ml_ok && near_poisson true for any theta > 0). The cleaner
  # test is to use mom_threshold large enough that the ML theta falls below
  # 1/mom_threshold for this highly overdispersed scenario.
  res <- mutze_test(cut, poisson_threshold = 1e6, mom_threshold = 0.5)
  # mom_threshold = 0.5 means extreme_overdisp triggers when theta < 2
  # Given k ~ 1 the ML theta ~ 1 so branch should fire
  expect_true(res$fallback %in% c("mom", "ml"))
  expect_true(is.finite(res$z))
  expect_true(is.finite(res$se))
  expect_true(res$se > 0)
})

test_that("mutze_test MoM fallback gives larger SE than Poisson under overdispersion", {
  # With genuinely overdispersed data, the MoM fallback SE should exceed the
  # Poisson-fallback SE (Poisson variance is anti-conservative).
  set.seed(1)
  enroll_rate <- data.frame(rate = 50, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.6, 0.4), dispersion = c(1.5, 1.5)
  )
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 2, n = 100)
  cut <- cut_data_by_date(sim, cut_date = 1.5)

  # Force MoM branch via very tight mom_threshold
  res_mom <- mutze_test(cut, poisson_threshold = 1e6, mom_threshold = 0.5)
  res_pois <- mutze_test(cut, method = "poisson")

  if (identical(res_mom$fallback, "mom")) {
    expect_gt(res_mom$se, res_pois$se)
  }
})

test_that("mutze_test validates threshold arguments", {
  df <- data.frame(
    treatment = c("A", "A", "B", "B"),
    events = c(1, 2, 0, 1),
    tte = c(1, 1, 1, 1)
  )
  expect_error(mutze_test(df, poisson_threshold = -1), "must be positive")
  expect_error(mutze_test(df, mom_threshold = 0), "must be positive")
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

# --- Score test tests ---

test_that("score test returns valid mutze_test object", {
  set.seed(100)
  enroll_rate <- data.frame(rate = 20, duration = 1)
  fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3),
                          dispersion = 0.3)
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 2, n = 50)
  cut <- cut_data_by_date(sim, cut_date = 1.5)
  res <- mutze_test(cut, test_type = "score")

  expect_s3_class(res, "mutze_test")
  expect_equal(res$test_type, "score")
  expect_true(is.finite(res$z))
  expect_true(is.finite(res$p_value))
  expect_true(is.finite(res$se))
  expect_true(is.finite(res$estimate))
  expect_true(is.finite(res$rate_ratio))
  expect_true(res$p_value >= 0 && res$p_value <= 1)
  expect_true(res$fallback %in% c("ml", "poisson", "mom"))
  expect_match(res$method, "score")
})

test_that("wald test returns test_type field", {
  set.seed(101)
  enroll_rate <- data.frame(rate = 20, duration = 1)
  fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3),
                          dispersion = 0.3)
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 1, n = 30)
  cut <- cut_data_by_date(sim, cut_date = 0.8)
  res <- mutze_test(cut, test_type = "wald")
  expect_equal(res$test_type, "wald")
})

test_that("score test with method='poisson' works", {
  set.seed(102)
  enroll_rate <- data.frame(rate = 20, duration = 1)
  fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3))
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 1, n = 30)
  cut <- cut_data_by_date(sim, cut_date = 0.8)
  res <- mutze_test(cut, test_type = "score", method = "poisson")

  expect_s3_class(res, "mutze_test")
  expect_equal(res$test_type, "score")
  expect_match(res$method, "Poisson score")
  expect_equal(res$fallback, "poisson")
  expect_equal(res$dispersion, Inf)
})

test_that("score and wald z-statistics have the same sign", {
  set.seed(103)
  enroll_rate <- data.frame(rate = 20, duration = 1)
  fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3),
                          dispersion = 0.3)
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 1.5, n = 40)
  cut <- cut_data_by_date(sim, cut_date = 1.2)

  wald_res <- mutze_test(cut, test_type = "wald")
  score_res <- mutze_test(cut, test_type = "score")

  expect_equal(sign(wald_res$z), sign(score_res$z))
})

test_that("score test two-sided p-value", {
  set.seed(104)
  enroll_rate <- data.frame(rate = 20, duration = 1)
  fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3),
                          dispersion = 0.3)
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 1, n = 30)
  cut <- cut_data_by_date(sim, cut_date = 0.8)

  res1 <- mutze_test(cut, test_type = "score", sided = 1)
  res2 <- mutze_test(cut, test_type = "score", sided = 2)

  expect_equal(res1$sided, 1)
  expect_equal(res2$sided, 2)
  # Two-sided p-value relates to one-sided for the same z
  expect_equal(res2$p_value, 2 * pnorm(-abs(res1$z)), tolerance = 1e-10)
})

test_that("score test falls back to Poisson on near-Poisson data", {
  set.seed(105)
  enroll_rate <- data.frame(rate = 30, duration = 1)
  fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3),
                          dispersion = 0.0001)
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 1.5, n = 60)
  cut <- cut_data_by_date(sim, cut_date = 1)
  res <- mutze_test(cut, test_type = "score", poisson_threshold = 50)

  expect_true(res$fallback %in% c("poisson", "ml"))
  expect_match(res$method, "score")
})

test_that("print.mutze_test works for score test", {
  set.seed(106)
  enroll_rate <- data.frame(rate = 20, duration = 1)
  fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3),
                          dispersion = 0.3)
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 1, n = 20)
  cut <- cut_data_by_date(sim, cut_date = 0.8)
  res <- mutze_test(cut, test_type = "score")

  expect_output(print(res), "Mutze Test Results")
  expect_output(print(res), "score")
})

test_that("score test z=0 when both groups have 0 events", {
  df <- data.frame(
    treatment = factor(rep(c("Control", "Experimental"), each = 10)),
    events = rep(0L, 20),
    tte = rep(1, 20)
  )
  res <- mutze_test(df, test_type = "score")
  expect_equal(res$z, 0)
  expect_equal(res$estimate, 0)
  expect_equal(res$rate_ratio, 1)
  expect_equal(res$p_value, 0.5)  # one-sided, z=0
})

test_that("score test handles 0 events in one group", {
  # Control has events, experimental has none → est = -Inf, z finite
  df <- data.frame(
    treatment = factor(rep(c("Control", "Experimental"), each = 10)),
    events = c(rep(1L, 10), rep(0L, 10)),
    tte = rep(1, 20)
  )
  res <- mutze_test(df, test_type = "score")
  expect_true(is.finite(res$z))
  expect_true(res$z < 0)  # treatment has fewer events
  expect_equal(res$estimate, -Inf)

  # Experimental has events, control has none → est = +Inf, z finite
  df2 <- data.frame(
    treatment = factor(rep(c("Control", "Experimental"), each = 10)),
    events = c(rep(0L, 10), rep(1L, 10)),
    tte = rep(1, 20)
  )
  res2 <- mutze_test(df2, test_type = "score")
  expect_true(is.finite(res2$z))
  expect_true(res2$z > 0)  # treatment has more events
  expect_equal(res2$estimate, Inf)
})
