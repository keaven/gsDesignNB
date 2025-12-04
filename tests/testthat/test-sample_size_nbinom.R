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
  expect_equal(res$exposure, 1)

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

  # Error handling
  expect_error(sample_size_nbinom(lambda1 = -1, lambda2 = 0.3, dispersion = 0.1, accrual_rate = 1, accrual_duration = 1, trial_duration = 1))
  expect_error(sample_size_nbinom(lambda1 = 0.5, lambda2 = 0.3, dispersion = -1, accrual_rate = 1, accrual_duration = 1, trial_duration = 1))

  # Variable accrual
  # Case 1: Uniform accrual over [0, 10], trial ends at 10. Avg exposure = 5.
  res_accrual <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 10, trial_duration = 10
  )
  expect_equal(res_accrual$exposure, 5)

  # Case 2: Piecewise accrual
  # Rate 10 for 5 units, Rate 20 for 5 units. Trial duration 15.
  # Seg 1: N=50, midpoint=2.5, avg_fup=15-2.5=12.5. Mass=50*12.5=625
  # Seg 2: N=100, midpoint=7.5 (5 + 2.5), avg_fup=15-7.5=7.5. Mass=100*7.5=750
  # Total N = 150. Total Mass = 1375. Avg = 1375/150 = 9.166667
  res_piecewise <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = c(10, 20), accrual_duration = c(5, 5), trial_duration = 15
  )
  expect_equal(res_piecewise$exposure, 1375 / 150)

  # Error checks for accrual
  expect_error(sample_size_nbinom(lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, accrual_rate = 10, trial_duration = 10)) # Missing duration
  expect_error(sample_size_nbinom(lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, accrual_rate = 10, accrual_duration = c(5, 5), trial_duration = 10)) # Mismatch
  expect_error(sample_size_nbinom(lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, accrual_rate = 10, accrual_duration = 11, trial_duration = 10)) # Accrual > Trial

  # Test Friede method
  res_friede <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8, method = "friede",
    accrual_rate = 10, accrual_duration = 2, trial_duration = 2
  )
  expect_type(res_friede, "list")
  expect_s3_class(res_friede, "sample_size_nbinom_result")
  expect_true(res_friede$n1 > 0)
  expect_true(res_friede$n2 > 0)

  # Friede and Zhu methods usually give similar but slightly different results
  # because variance estimation details differ slightly or are equivalent under specific conditions
  res_zhu <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8, method = "zhu",
    accrual_rate = 10, accrual_duration = 2, trial_duration = 2
  )

  # In this implementation:
  # Zhu: V = (1/mu1 + k) + 1/r * (1/mu2 + k)
  # Friede: V_bar = (1/mu1 + k)/p1 + (1/mu2 + k)/p2
  # If r=1, p1=0.5, p2=0.5.
  # Friede V_bar = 2(1/mu1 + k) + 2(1/mu2 + k) = 2 * [ (1/mu1+k) + (1/mu2+k) ]
  # Zhu V = (1/mu1+k) + (1/mu2+k)
  # Zhu n1 = (z+z)^2 * V / log^2
  # Friede n_total = (z+z)^2 * V_bar / log^2 = (z+z)^2 * 2*V / log^2
  # Friede n1 = n_total/2 = (z+z)^2 * V / log^2
  # So for r=1 they should be identical.
  expect_equal(res_zhu$n1, res_friede$n1)
  expect_equal(res_zhu$n2, res_friede$n2)

  # Test unequal allocation r=2
  # p1=1/3, p2=2/3.
  # Friede V_bar = 3(1/mu1+k) + 1.5(1/mu2+k)
  # Zhu V = (1/mu1+k) + 0.5(1/mu2+k)
  # Friede n_total = C * V_bar. n1 = C * V_bar / 3 = C * ( (1/mu1+k) + 0.5(1/mu2+k) ) = C * Zhu_V
  # So they should still be identical for n1.
  res_zhu_r2 <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, ratio = 2, method = "zhu",
    accrual_rate = 10, accrual_duration = 2, trial_duration = 2
  )
  res_friede_r2 <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, ratio = 2, method = "friede",
    accrual_rate = 10, accrual_duration = 2, trial_duration = 2
  )

  expect_equal(res_zhu_r2$n1, res_friede_r2$n1)

  # Error check for unknown method
  expect_error(sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, method = "unknown",
    accrual_rate = 10, accrual_duration = 2, trial_duration = 2
  ))
})
