test_that("gsNBCalendar creates valid gsNB object", {
  # Create sample size object
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )

  # Create group sequential design
  gs_design <- gsNBCalendar(nb_ss, k = 3, test.type = 4)

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

  gs_design <- gsNBCalendar(nb_ss, k = 2)

  # Check ratio is preserved
  expect_equal(gs_design$n2[2] / gs_design$n1[2], 2)
})

test_that("gsNBCalendar works with different test types", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )

  # One-sided
  gs1 <- gsNBCalendar(nb_ss, k = 2, test.type = 1)
  expect_s3_class(gs1, "gsNB")

  # Two-sided symmetric
  gs2 <- gsNBCalendar(nb_ss, k = 2, test.type = 2)
  expect_s3_class(gs2, "gsNB")

  # Two-sided asymmetric non-binding
  gs4 <- gsNBCalendar(nb_ss, k = 2, test.type = 4)
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
    sflpar = -2
  )

  expect_s3_class(gs_design, "gsNB")
})
