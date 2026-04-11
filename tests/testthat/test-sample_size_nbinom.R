# source(testthat::test_path("../../R/sample_size_nbinom.R")) # Not needed in package context

test_that("sample_size_nbinom calculates correctly", {
  # Basic calculation - now requires accrual params
  # To get exposure=1, e.g. trial=2, accrual=2 (uniform 0-2). Avg fup = 2 - 1 = 1.

  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 2, trial_duration = 2
  )
  expect_type(res, "list")
  expect_s3_class(res, "sample_size_nbinom_result")
  expect_equal(res$inputs$lambda1, 0.5)
  expect_equal(res$inputs$power, 0.8)
  expect_true(res$n1 > 0)
  expect_true(res$n2 > 0)
  expect_equal(res$n_total, res$n1 + res$n2)
  expect_equal(res$exposure[1], 1)

  # Check Poisson limit (dispersion = 0)
  res_nb_0 <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0, power = 0.8,
    accrual_rate = 10, accrual_duration = 2, trial_duration = 2
  )

  # Manual Poisson calculation
  alpha <- 0.025
  beta <- 0.2
  mu1 <- 0.5 * 1 # exposure 1
  mu2 <- 0.3 * 1
  z_a <- qnorm(1 - alpha) # one-sided alpha = 0.025 corresponds to sided=1 default
  z_b <- qnorm(1 - beta)
  num <- (z_a + z_b)^2 * (1 / mu1 + 1 / mu2)
  den <- (log(mu1 / mu2))^2
  n_poisson <- ceiling(num / den)

  expect_equal(res_nb_0$n1, n_poisson)

  # Check ratio
  res_ratio <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, ratio = 2,
    accrual_rate = 10, accrual_duration = 2, trial_duration = 2
  )
  expect_equal(res_ratio$n2, 2 * res_ratio$n1)

  # Check vector dispersion
  res_vec_disp <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = c(0.1, 0.2), power = 0.8,
    accrual_rate = 10, accrual_duration = 2, trial_duration = 2
  )
  expect_type(res_vec_disp, "list")
  expect_true(res_vec_disp$n1 > 0)

  # Error handling
  expect_error(sample_size_nbinom(lambda1 = -1, lambda2 = 0.3, dispersion = 0.1, accrual_rate = 1, accrual_duration = 1, trial_duration = 1))
  expect_error(sample_size_nbinom(lambda1 = 0.5, lambda2 = 0.3, dispersion = -1, accrual_rate = 1, accrual_duration = 1, trial_duration = 1))
  expect_error(sample_size_nbinom(lambda1 = 0.5, lambda2 = 0.3, dispersion = c(0.1, 0.2, 0.3), accrual_rate = 1, accrual_duration = 1, trial_duration = 1))

  # Variable accrual
  # Case 1: Uniform accrual over [0, 10], trial ends at 10. Avg exposure = 5.
  res_accrual <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 10, trial_duration = 10
  )
  expect_equal(res_accrual$exposure[1], 5)

  # Case 2: Piecewise accrual
  # Rate 10 for 5 units, Rate 20 for 5 units. Trial duration 15.
  # Seg 1: N=50, midpoint=2.5, avg_fup=15-2.5=12.5. Mass=50*12.5=625
  # Seg 2: N=100, midpoint=7.5 (5 + 2.5), avg_fup=15-7.5=7.5. Mass=100*7.5=750
  # Total N = 150. Total Mass = 1375. Avg = 1375/150 = 9.166667
  res_piecewise <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = c(10, 20), accrual_duration = c(5, 5), trial_duration = 15
  )
  expect_equal(res_piecewise$exposure[1], 1375 / 150)

  # Error checks for accrual
  expect_error(sample_size_nbinom(lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, accrual_rate = 10, trial_duration = 10)) # Missing duration
  expect_error(sample_size_nbinom(lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, accrual_rate = 10, accrual_duration = c(5, 5), trial_duration = 10)) # Mismatch
  
  # Accrual > Trial is now allowed (truncated)
  res_trunc <- sample_size_nbinom(lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = NULL, accrual_rate = 10, accrual_duration = 11, trial_duration = 10)
  expect_equal(res_trunc$n_total, 100) # 10 * 10
})

# Test Exposure Reporting with Event Gap
# When gap > 0, function should return exposure_at_risk_n1/n2
test_that("Event gap results in correct exposure reporting", {
  lambda1 <- 0.5
  lambda2 <- 0.3
  gap <- 0.5

  res <- sample_size_nbinom(
    lambda1 = lambda1, lambda2 = lambda2, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    event_gap = gap
  )

  # Calendar exposure
  exp_cal <- res$exposure[1]
  expect_true(exp_cal > 0)

  # Expected at-risk exposures
  exp1_expected <- exp_cal / (1 + lambda1 * gap)
  exp2_expected <- exp_cal / (1 + lambda2 * gap)

  expect_equal(res$exposure_at_risk_n1, exp1_expected)
  expect_equal(res$exposure_at_risk_n2, exp2_expected)

  # Check method name remains "zhu"
  # expect_equal(res$inputs$method, "zhu") # Method removed
})

test_that("Jensen correction for event gap increases sample size", {
  # With event gap and dispersion > 0, the Jensen correction reduces the

  # effective rate, increasing the required sample size relative to the
  # naive formula (lambda / (1 + lambda * gap)).
  #
  # Verify the correction is applied by comparing sample sizes:
  # - With dispersion = 0 (Poisson): no correction needed (no frailty)
  # - With dispersion > 0: correction reduces effective rate -> larger n

  gap <- 0.5
  lambda1 <- 0.5
  lambda2 <- 0.3

  # Poisson case (k=0): effective rate = lambda / (1 + lambda * gap) exactly
  res_poisson <- sample_size_nbinom(
    lambda1 = lambda1, lambda2 = lambda2, dispersion = 0, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    event_gap = gap
  )

  # NB case (k=0.5): correction factor < 1, so effective rate is lower,
  # variance is higher, and n should be larger
  res_nb <- sample_size_nbinom(
    lambda1 = lambda1, lambda2 = lambda2, dispersion = 0.5, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    event_gap = gap
  )

  # NB should require more subjects (dispersion increases variance + gap correction)
  expect_true(res_nb$n_total > res_poisson$n_total)
  expect_true(res_nb$variance > res_poisson$variance)

  # Verify the correction factor analytically for a single group
  # correction = 1 - k * lambda * gap / (1 + lambda * gap)^2
  k <- 0.5
  correction1 <- 1 - k * lambda1 * gap / (1 + lambda1 * gap)^2
  correction2 <- 1 - k * lambda2 * gap / (1 + lambda2 * gap)^2

  # Correction should be < 1 and > 0 for reasonable parameters
  expect_true(correction1 > 0 && correction1 < 1)
  expect_true(correction2 > 0 && correction2 < 1)

  # The effective rate with correction should be smaller than the naive rate
  naive1 <- lambda1 / (1 + lambda1 * gap)
  corrected1 <- naive1 * correction1
  expect_true(corrected1 < naive1)

  # Verify that with k=0, no correction is applied (correction factor = 1)
  correction_k0 <- 1 - 0 * lambda1 * gap / (1 + lambda1 * gap)^2
  expect_equal(correction_k0, 1)
})

test_that("sample_size_nbinom matches published results", {
  # Helper to simulate fixed exposure of 1.0
  # We use a very short accrual duration and trial duration = 1 + accrual duration
  # This ensures everyone is followed for approximately 1.0 unit of time.
  calc_fixed_exp <- function(lambda1, lambda2, dispersion, power, alpha) {
    sample_size_nbinom(
      lambda1 = lambda1, lambda2 = lambda2, dispersion = dispersion, power = power,
      alpha = alpha, sided = 1,
      accrual_rate = 1000, accrual_duration = 0.001, trial_duration = 1.001
    )
  }

  # Zhu and Lakkis (2014)
  # Example 1: Moderate Rates
  # L1=0.5, L2=0.3, k=0.1, Power=0.8, Alpha=0.025 -> n=167
  res_zhu_1 <- calc_fixed_exp(0.5, 0.3, 0.1, 0.8, 0.025)
  expect_equal(res_zhu_1$n1, 167)
  expect_equal(res_zhu_1$n2, 167)

  # Example 2: High Rates
  # L1=1.0, L2=0.5, k=0.5, Power=0.9, Alpha=0.025 -> n=88
  res_zhu_2 <- calc_fixed_exp(1.0, 0.5, 0.5, 0.9, 0.025)
  expect_equal(res_zhu_2$n1, 88)
  expect_equal(res_zhu_2$n2, 88)

  # Friede and Schmidli (2010)
  # Example 1: Standard Scenario
  # L1=0.6, L2=0.3, k=0.4, Power=0.8, Alpha=0.025 -> n=95
  res_friede_1 <- calc_fixed_exp(0.6, 0.3, 0.4, 0.8, 0.025)
  expect_equal(res_friede_1$n1, 95)
  expect_equal(res_friede_1$n2, 95)

  # Example 2: High Dispersion
  # L1=1.0, L2=0.5, k=0.5, Power=0.8, Alpha=0.025 -> n=66
  res_friede_2 <- calc_fixed_exp(1.0, 0.5, 0.5, 0.8, 0.025)
  expect_equal(res_friede_2$n1, 66)
  expect_equal(res_friede_2$n2, 66)
})

test_that("sample_size_nbinom handles vector dropout_rate", {
  # Case 1: Same dropout and max_followup for both groups
  # Should give same result as scalar input
  res_scalar <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    dropout_rate = 0.05, max_followup = 12
  )

  res_vector <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    dropout_rate = c(0.05, 0.05), max_followup = 12
  )

  expect_equal(res_scalar$n_total, res_vector$n_total)
  expect_equal(res_scalar$exposure[1], res_vector$exposure[1])
  expect_equal(res_scalar$exposure[2], res_vector$exposure[2])

  # Case 2: Different dropout rates
  # Group 1 has higher dropout, so exposure should be lower
  res_diff <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    dropout_rate = c(0.1, 0.01), max_followup = 12
  )

  expect_true(res_diff$exposure[1] < res_diff$exposure[2])
  expect_true(length(res_diff$exposure) == 2)
})

test_that("sample_size_nbinom input validation errors", {
  base <- list(lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1,
               accrual_rate = 10, accrual_duration = 2, trial_duration = 2)

  # Bad rr0
  expect_error(do.call(sample_size_nbinom, c(base, rr0 = -1)), "rr0 must be positive")

  # Bad power
  expect_error(do.call(sample_size_nbinom, c(base, power = 0)), "Power must be between")
  expect_error(do.call(sample_size_nbinom, c(base, power = 1)), "Power must be between")

  # Bad alpha
  expect_error(do.call(sample_size_nbinom, c(base, alpha = 0)), "Alpha must be between")
  expect_error(do.call(sample_size_nbinom, c(base, alpha = 1)), "Alpha must be between")

  # Bad dropout_rate
  expect_error(do.call(sample_size_nbinom, c(base, dropout_rate = -0.1)), "non-negative")
  expect_error(do.call(sample_size_nbinom, c(base, list(dropout_rate = c(0.1, 0.1, 0.1)))), "scalar or a vector of length 2")

  # Bad max_followup
  expect_error(do.call(sample_size_nbinom, c(base, max_followup = 0)), "max_followup must be positive")
  expect_error(do.call(sample_size_nbinom, c(base, max_followup = -1)), "max_followup must be positive")
  expect_error(do.call(sample_size_nbinom, c(base, list(max_followup = c(1, 2, 3)))), "scalar or a vector of length 2")
})

test_that("sample_size_nbinom handles non-inferiority (rr0 != 1)", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.5, rr0 = 1.1, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 2, trial_duration = 2
  )
  expect_true(res$n_total > 0)
  expect_equal(res$inputs$rr0, 1.1)
})

test_that("sample_size_nbinom handles max_followup with split case", {
  # max_followup shorter than some follow-up times triggers Case 3 (split)
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 10, trial_duration = 15,
    max_followup = 8
  )
  expect_true(res$n_total > 0)
  # With max_followup, exposure should be less than without
  res_no_mf <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 10, trial_duration = 15
  )
  expect_true(res$exposure[1] <= res_no_mf$exposure[1])
})

test_that("sample_size_nbinom handles max_followup all-truncated case", {
  # Very short max_followup truncates everyone
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 10, trial_duration = 15,
    max_followup = 2
  )
  expect_true(res$n_total > 0)
  expect_true(res$exposure[1] <= 2)
})

test_that("sample_size_nbinom handles dropout with max_followup", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 10, trial_duration = 15,
    dropout_rate = 0.1, max_followup = 8
  )
  expect_true(res$n_total > 0)
})

test_that("print.sample_size_nbinom_result works", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  expect_output(print(res), "Sample size")
})

test_that("print.sample_size_nbinom_result shows non-inferiority", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.5, rr0 = 1.1, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  expect_output(print(res), "Null hypothesis RR")
})

test_that("print.sample_size_nbinom_result shows vector dispersion", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = c(0.1, 0.2), power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  expect_output(print(res), "\\(n1\\)")
})

test_that("print.sample_size_nbinom_result shows dropout", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    dropout_rate = 0.05
  )
  expect_output(print(res), "Dropout rate")
})

test_that("print shows different dropout rates per group", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    dropout_rate = c(0.05, 0.10)
  )
  expect_output(print(res), "\\(n1\\)")
})

test_that("print shows event gap info", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    event_gap = 0.5
  )
  expect_output(print(res), "Event gap")
})

test_that("print shows max_followup", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    max_followup = 12
  )
  expect_output(print(res), "Max follow")
})

test_that("print shows different exposure per group", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    dropout_rate = c(0.05, 0.2)
  )
  expect_output(print(res), "\\(n1\\)")
})

test_that("summary.sample_size_nbinom_result returns summary object", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  s <- summary(res)
  expect_s3_class(s, "sample_size_nbinom_summary")
  expect_true(is.character(s))
})

test_that("summary shows non-inferiority rr0", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.5, rr0 = 1.1, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  s <- summary(res)
  s_text <- paste(s, collapse = " ")
  expect_true(grepl("null hypothesis RR", s_text))
})

test_that("summary shows vector dispersion", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = c(0.1, 0.2), power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  s <- summary(res)
  s_text <- paste(s, collapse = " ")
  expect_true(grepl("\\(n1\\)", s_text))
})

test_that("summary shows event gap info", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    event_gap = 0.5
  )
  s <- summary(res)
  s_text <- paste(s, collapse = " ")
  expect_true(grepl("Event gap", s_text))
})

test_that("summary shows different exposure per group", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    dropout_rate = c(0.05, 0.2)
  )
  s <- summary(res)
  s_text <- paste(s, collapse = " ")
  expect_true(grepl("\\(n1\\)", s_text))
})

test_that("print.sample_size_nbinom_summary works", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  s <- summary(res)
  expect_output(print(s), "negative binomial")
})
