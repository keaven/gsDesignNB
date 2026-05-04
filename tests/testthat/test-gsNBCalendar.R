test_that("gsNBCalendar creates valid gsNB object", {
  # Create sample size object
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )

  # Create group sequential design
  gs_design <- gsNBCalendar(nb_ss, k = 3, test.type = 4, analysis_times = c(12, 18, 24))

  # Check class inheritance
  expect_s3_class(gs_design, "gsNB")
  expect_s3_class(gs_design, "gsDesign")
  expect_s3_class(gs_design, "sample_size_nbinom_result")

  # Check that nb_design is preserved
  expect_identical(gs_design$nb_design, nb_ss)

  # Check sample size vectors
  expect_length(gs_design$n1, 3)
  expect_length(gs_design$n2, 3)
  expect_length(gs_design$n_total, 3)

  # Cumulative sample sizes should increase
  expect_true(all(diff(gs_design$n1) > 0))
  expect_true(all(diff(gs_design$n2) > 0))

  # Final sample sizes should match original (approximately, due to inflation)
  expect_equal(
    gs_design$n1[3] + gs_design$n2[3],
    gs_design$n_total[3]
  )
})

test_that("gsNBCalendar rejects invalid input", {
  # Not a sample_size_nbinom_result object
  expect_error(
    gsNBCalendar(list(n_total = 100)),
    "must be an object of class 'sample_size_nbinom_result'"
  )
})

test_that("gsNBCalendar respects allocation ratio", {
  # Create sample size with 2:1 allocation
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    ratio = 2,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )

  gs_design <- gsNBCalendar(nb_ss, k = 2, analysis_times = c(12, 24))

  # Check ratio is preserved
  expect_equal(gs_design$n2[2] / gs_design$n1[2], 2)
})

test_that("gsNBCalendar works with different test types", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )

  # One-sided
  gs1 <- gsNBCalendar(nb_ss, k = 2, test.type = 1, analysis_times = c(12, 24))
  expect_s3_class(gs1, "gsNB")

  # Two-sided symmetric
  gs2 <- gsNBCalendar(nb_ss, k = 2, test.type = 2, analysis_times = c(12, 24))
  expect_s3_class(gs2, "gsNB")

  # Two-sided asymmetric non-binding
  gs4 <- gsNBCalendar(nb_ss, k = 2, test.type = 4, analysis_times = c(12, 24))
  expect_s3_class(gs4, "gsNB")
})

test_that("gsNBCalendar works with custom spending functions", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )

  # O'Brien-Fleming-like spending
  gs_design <- gsNBCalendar(
    nb_ss,
    k = 3,
    sfu = gsDesign::sfHSD,
    sfupar = -4,
    sfl = gsDesign::sfHSD,
    sflpar = -2,
    analysis_times = c(12, 18, 24)
  )

  expect_s3_class(gs_design, "gsNB")
})

test_that("summary.gsNB returns gsNBsummary object", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  gs_design <- gsNBCalendar(nb_ss, k = 3, analysis_times = c(12, 18, 24))

  s <- summary(gs_design)
  expect_s3_class(s, "gsNBsummary")
  expect_true(is.character(s))

  s_collapsed <- paste(s, collapse = " ")
  expect_true(nchar(s_collapsed) > 0)

  # Should mention key terms
  expect_true(grepl("negative binomial", s_collapsed, ignore.case = TRUE))
  expect_true(grepl("0.5", s_collapsed))
  expect_true(grepl("0.3", s_collapsed))
})

test_that("print.gsNBsummary prints without error", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  gs_design <- gsNBCalendar(nb_ss, k = 2, analysis_times = c(12, 24))
  s <- summary(gs_design)

  expect_output(print(s))
})

test_that("toInteger.gsNB rounds final sample size", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  gs_design <- gsNBCalendar(nb_ss, k = 2, analysis_times = c(12, 24))

  gs_int <- toInteger(gs_design)

  expect_s3_class(gs_int, "gsNB")

  # Final sample size should be integer (within rounding)
  k <- gs_int$k
  expect_equal(gs_int$n_total[k], round(gs_int$n_total[k]))

  # n1 + n2 should equal n_total
  expect_equal(gs_int$n1[k] + gs_int$n2[k], gs_int$n_total[k])
})

test_that("toInteger.gsNB respects ratio", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    ratio = 2,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  gs_design <- gsNBCalendar(nb_ss, k = 2, analysis_times = c(12, 24))

  gs_int <- toInteger(gs_design)

  # Final total should be divisible by (ratio + 1) = 3
  expect_equal(gs_int$n_total[2] %% 3, 0)
})

test_that("toInteger.gsNB preserves calendar enrollment at interim analyses", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 12, trial_duration = 24
  )
  analysis_times <- c(10, 18, 24)
  gs_design <- gsNBCalendar(nb_ss, k = 3, analysis_times = analysis_times)

  gs_int <- toInteger(gs_design)

  expected_n <- gs_int$n_total[gs_int$k] *
    pmin(analysis_times, nb_ss$accrual_duration) / nb_ss$accrual_duration

  expect_equal(gs_int$n_total, expected_n, tolerance = 1e-8)
  expect_true(all(diff(gs_int$n.I) > 0))
})

test_that("toInteger.gsNB preserves piecewise accrual shape", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.35, dispersion = 0.5, power = 0.9,
    accrual_rate = c(12, 24), accrual_duration = c(6, 6),
    trial_duration = 24, dropout_rate = 0.10 / 12,
    max_followup = 12, event_gap = 20 / 30.4375,
    test_type = "score"
  )
  analysis_times <- c(9, 14, 24)
  gs_design <- gsNBCalendar(nb_ss, k = 3, analysis_times = analysis_times)

  gs_int <- toInteger(gs_design)

  information_scale <- gs_design$n.I[gs_design$k] / (1 / nb_ss$variance)
  rounding_scale <- gs_int$n_total[gs_int$k] / gs_design$n_total[gs_design$k]
  expected_n <- vapply(
    analysis_times,
    function(t) {
      sample_size_nbinom(
        lambda1 = nb_ss$inputs$lambda1,
        lambda2 = nb_ss$inputs$lambda2,
        rr0 = nb_ss$inputs$rr0,
        dispersion = nb_ss$inputs$dispersion,
        power = NULL,
        alpha = nb_ss$inputs$alpha,
        sided = nb_ss$inputs$sided,
        ratio = nb_ss$inputs$ratio,
        accrual_rate = nb_ss$accrual_rate * information_scale * rounding_scale,
        accrual_duration = nb_ss$accrual_duration,
        trial_duration = t,
        dropout_rate = nb_ss$inputs$dropout_rate,
        max_followup = nb_ss$inputs$max_followup,
        event_gap = nb_ss$inputs$event_gap
      )$n_total
    },
    numeric(1)
  )

  expect_equal(gs_int$n_total, expected_n, tolerance = 1e-8)
  expect_lt(gs_int$n_total[gs_int$k], 1.05 * gs_design$n_total[gs_design$k])
})

test_that("toInteger.gsDesign dispatches correctly", {
  gs <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 2)
  gs_int <- toInteger(gs)
  expect_s3_class(gs_int, "gsDesign")
})

test_that("compute_info_at_time returns positive value", {
  info <- compute_info_at_time(
    analysis_time = 12,
    accrual_rate = 10,
    accrual_duration = 10,
    lambda1 = 0.5,
    lambda2 = 0.3,
    dispersion = 0.1
  )

  expect_true(is.numeric(info))
  expect_true(info > 0)
})

test_that("compute_info_at_time increases with analysis time", {
  base_args <- list(
    accrual_rate = 10,
    accrual_duration = 10,
    lambda1 = 0.5,
    lambda2 = 0.3,
    dispersion = 0.1
  )

  info_early <- do.call(compute_info_at_time, c(list(analysis_time = 6), base_args))
  info_late <- do.call(compute_info_at_time, c(list(analysis_time = 12), base_args))

  expect_true(info_late > info_early)
})

test_that("gsNBCalendar errors without analysis_times", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )

  expect_error(
    gsNBCalendar(nb_ss, k = 3),
    "analysis_times must be provided"
  )
})

test_that("gsNBCalendar errors when analysis_times length != k", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )

  expect_error(
    gsNBCalendar(nb_ss, k = 3, analysis_times = c(12, 24)),
    "analysis_times must have length k"
  )
})
