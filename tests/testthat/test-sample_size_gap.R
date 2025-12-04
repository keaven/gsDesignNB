
test_that("sample_size_nbinom adjusts for event_gap", {
  # Fixed design parameters for testing
  acc_r <- 1
  acc_d <- 1
  trial_d <- 2 # Exposure approx 1 if no dropout?
  
  # No gap
  res_no_gap <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.5, dispersion = 0.1, power = 0.8, 
    accrual_rate = acc_r, accrual_duration = acc_d, trial_duration = trial_d
  )
  
  # Use different lambdas for power calc
  l1 <- 0.5
  l2 <- 0.3
  gap <- 0.1
  
  res1 <- sample_size_nbinom(
    lambda1 = l1, lambda2 = l2, dispersion = 0, power = 0.8,
    accrual_rate = acc_r, accrual_duration = acc_d, trial_duration = trial_d
  )
  
  res2 <- sample_size_nbinom(
    lambda1 = l1, lambda2 = l2, dispersion = 0, power = 0.8,
    accrual_rate = acc_r, accrual_duration = acc_d, trial_duration = trial_d,
    event_gap = gap
  )
  expect_equal(res2$inputs$event_gap, gap)
  
  expect_true(res2$n_total > res1$n_total)
  
  # Manually check effective rates
  eff_l1 <- l1 / (1 + l1 * gap)
  eff_l2 <- l2 / (1 + l2 * gap)
  
  res3 <- sample_size_nbinom(
    lambda1 = eff_l1, lambda2 = eff_l2, dispersion = 0, power = 0.8,
    accrual_rate = acc_r, accrual_duration = acc_d, trial_duration = trial_d
  )
  
  # res2 should be smaller than res3 (effect size logic)
  expect_true(res2$n_total < res3$n_total)
  
  # Calculated events should match
  mu1_eff <- eff_l1 * res2$exposure
  expect_equal(res2$events_n1, res2$n1 * mu1_eff)
})
