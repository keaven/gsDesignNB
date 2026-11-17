assert("mean follow-up recovers a known exponential dropout hazard", {
  for (cap in c(1, 12)) {
    for (hazard in c(0.001, 0.1, 1, 10)) {
      average <- -expm1(-hazard * cap) / hazard
      (all.equal(dropout_rate_from_followup(average, cap), hazard,
                 tolerance = 1e-9))
    }
  }
  (dropout_rate_from_followup(1, 1) == 0)
  (all.equal(dropout_rate_from_followup(10.8, 12),
             dropout_rate_from_followup(0.9, 1) / 12))
  delta <- dropout_rate_from_followup(0.9, 1)
  design <- sample_size_nbinom(
    lambda1 = 0.66, lambda2 = 0.528, dispersion = 0.8, power = 0.9,
    accrual_rate = 100, accrual_duration = 1, trial_duration = 2,
    max_followup = 1, dropout_rate = delta
  )
  (all.equal(design$exposure, c(0.9, 0.9)))
})

assert("follow-up conversion rejects invalid inputs", {
  for (bad in list(NULL, numeric(), NA_real_, NaN, Inf, 0, -1,
                   c(0.8, 0.9), "0.9", TRUE)) {
    (has_error(dropout_rate_from_followup(bad, 1)))
    (has_error(dropout_rate_from_followup(0.9, bad)))
  }
  (has_error(dropout_rate_from_followup(1.1, 1)))
})

assert("constant-exposure sizing reproduces the PASS published example", {
  # NCSS, Tests for the Ratio of Two Negative Binomial Rates, Example 2:
  # https://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/PASS/Tests_for_the_Ratio_of_Two_Negative_Binomial_Rates.pdf
  # PASS M3 uses pooled null rates; here this is test_type = "score".
  design <- sample_size_nbinom(
    lambda1 = 0.66, lambda2 = 0.66 * 0.8, dispersion = 0.8,
    power = 0.9, alpha = 0.05, sided = 2, test_type = "score",
    accrual_rate = 100, accrual_duration = 1, trial_duration = 2,
    max_followup = 0.9
  )
  (design$n1 == 1131)
  (design$n2 == 1131)
  # Independently check the true-rate variance (PASS M2 / package Wald).
  v <- (1 / 0.66 + 1 / (0.66 * 0.8)) / 0.9 + 2 * 0.8
  expected <- ceiling((qnorm(0.975) + qnorm(0.9))^2 * v / log(0.8)^2)
  wald <- sample_size_nbinom(
    lambda1 = 0.66, lambda2 = 0.66 * 0.8, dispersion = 0.8,
    power = 0.9, alpha = 0.05, sided = 2, test_type = "wald",
    accrual_rate = 100, accrual_duration = 1, trial_duration = 2,
    max_followup = 0.9
  )
  (wald$n1 == expected)
})
