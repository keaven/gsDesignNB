
test_that("sample_size_nbinom truncates accrual correctly", {
  # Case 1: Trial duration < Accrual duration (scalar)
  # Accrual: 10/month for 10 months. Trial ends at 5 months.
  # Expected: Accrual 10/month for 5 months. Total N = 50.
  res1 <- sample_size_nbinom(
    lambda1 = 0.1, lambda2 = 0.1, dispersion = 0.1, power = NULL,
    accrual_rate = 10, accrual_duration = 10, trial_duration = 5
  )
  
  expect_equal(res1$n_total, 50)
  expect_equal(res1$accrual_duration, 5)
  
  # Case 2: Trial duration < Accrual duration (vector)
  # Accrual: 10/month for 5 months, then 20/month for 5 months.
  # Trial ends at 8 months.
  # Expected: 10/month for 5 months (50) + 20/month for 3 months (60) = 110.
  res2 <- sample_size_nbinom(
    lambda1 = 0.1, lambda2 = 0.1, dispersion = 0.1, power = NULL,
    accrual_rate = c(10, 20), accrual_duration = c(5, 5), trial_duration = 8
  )
  
  expect_equal(res2$n_total, 110)
  expect_equal(res2$accrual_duration, c(5, 3))
  
  # Case 3: Trial duration exactly at a breakpoint
  # Trial ends at 5 months.
  # Expected: 10/month for 5 months (50).
  res3 <- sample_size_nbinom(
    lambda1 = 0.1, lambda2 = 0.1, dispersion = 0.1, power = NULL,
    accrual_rate = c(10, 20), accrual_duration = c(5, 5), trial_duration = 5
  )
  
  expect_equal(res3$n_total, 50)
  expect_equal(res3$accrual_duration, 5)
  
  # Case 4: Trial duration > Total accrual duration
  # Trial ends at 12 months.
  # Expected: Full accrual. 10*5 + 20*5 = 150.
  res4 <- sample_size_nbinom(
    lambda1 = 0.1, lambda2 = 0.1, dispersion = 0.1, power = NULL,
    accrual_rate = c(10, 20), accrual_duration = c(5, 5), trial_duration = 12
  )
  
  expect_equal(res4$n_total, 150)
  expect_equal(res4$accrual_duration, c(5, 5))
})
