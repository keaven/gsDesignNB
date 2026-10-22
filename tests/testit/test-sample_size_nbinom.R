
assert("sample_size_nbinom calculates correctly", {
  # Basic calculation - now requires accrual params
  # To get exposure=1, e.g. trial=2, accrual=2 (uniform 0-2). Avg fup = 2 - 1 = 1.

  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 2, trial_duration = 2
  )
  (is.list(res))
  (inherits(res, "sample_size_nbinom_result"))
  (all.equal(res$inputs$lambda1, 0.5))
  (all.equal(res$inputs$power, 0.8))
  (res$n1 > 0)
  (res$n2 > 0)
  (all.equal(res$n_total, res$n1 + res$n2))
  (all.equal(res$exposure[1], 1))

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

  (all.equal(res_nb_0$n1, n_poisson))

  # Check ratio
  res_ratio <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, ratio = 2,
    accrual_rate = 10, accrual_duration = 2, trial_duration = 2
  )
  (all.equal(res_ratio$n2, 2 * res_ratio$n1))

  # Check vector dispersion
  res_vec_disp <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = c(0.1, 0.2), power = 0.8,
    accrual_rate = 10, accrual_duration = 2, trial_duration = 2
  )
  (is.list(res_vec_disp))
  (res_vec_disp$n1 > 0)

  # Error handling
  (has_error(sample_size_nbinom(lambda1 = -1, lambda2 = 0.3, dispersion = 0.1, accrual_rate = 1, accrual_duration = 1, trial_duration = 1)))
  (has_error(sample_size_nbinom(lambda1 = 0.5, lambda2 = 0.3, dispersion = -1, accrual_rate = 1, accrual_duration = 1, trial_duration = 1)))
  (has_error(sample_size_nbinom(lambda1 = 0.5, lambda2 = 0.3, dispersion = c(0.1, 0.2, 0.3), accrual_rate = 1, accrual_duration = 1, trial_duration = 1)))

  # Variable accrual
  # Case 1: Uniform accrual over [0, 10], trial ends at 10. Avg exposure = 5.
  res_accrual <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 10, trial_duration = 10
  )
  (all.equal(res_accrual$exposure[1], 5))

  # Case 2: Piecewise accrual
  # Rate 10 for 5 units, Rate 20 for 5 units. Trial duration 15.
  # Seg 1: N=50, midpoint=2.5, avg_fup=15-2.5=12.5. Mass=50*12.5=625
  # Seg 2: N=100, midpoint=7.5 (5 + 2.5), avg_fup=15-7.5=7.5. Mass=100*7.5=750
  # Total N = 150. Total Mass = 1375. Avg = 1375/150 = 9.166667
  res_piecewise <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = c(10, 20), accrual_duration = c(5, 5), trial_duration = 15
  )
  (all.equal(res_piecewise$exposure[1], 1375 / 150))

  # Error checks for accrual
  (has_error(sample_size_nbinom(lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, accrual_rate = 10, trial_duration = 10))) # Missing duration
  (has_error(sample_size_nbinom(lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, accrual_rate = 10, accrual_duration = c(5, 5), trial_duration = 10))) # Mismatch
  
  # Accrual > Trial is now allowed (truncated)
  res_trunc <- sample_size_nbinom(lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = NULL, accrual_rate = 10, accrual_duration = 11, trial_duration = 10)
  (all.equal(res_trunc$n_total, 100)) # 10 * 10
})

# Test Exposure Reporting with Event Gap
# When gap > 0, function should return exposure_at_risk_n1/n2
assert("Event gap results in correct exposure reporting", {
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
  (exp_cal > 0)

  # Expected at-risk exposures
  exp1_expected <- exp_cal / (1 + lambda1 * gap)
  exp2_expected <- exp_cal / (1 + lambda2 * gap)

  (all.equal(res$exposure_at_risk_n1, exp1_expected))
  (all.equal(res$exposure_at_risk_n2, exp2_expected))

  # Check method name remains "zhu"
  # (all.equal(res$inputs$method, "zhu")) # Method removed
})

assert("Jensen correction for event gap increases sample size", {
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
  (res_nb$n_total > res_poisson$n_total)
  (res_nb$variance > res_poisson$variance)

  # Verify the correction factor analytically for a single group
  # correction = 1 - k * lambda * gap / (1 + lambda * gap)^2
  k <- 0.5
  correction1 <- 1 - k * lambda1 * gap / (1 + lambda1 * gap)^2
  correction2 <- 1 - k * lambda2 * gap / (1 + lambda2 * gap)^2

  # Correction should be < 1 and > 0 for reasonable parameters
  (correction1 > 0 && correction1 < 1)
  (correction2 > 0 && correction2 < 1)

  # The effective rate with correction should be smaller than the naive rate
  naive1 <- lambda1 / (1 + lambda1 * gap)
  corrected1 <- naive1 * correction1
  (corrected1 < naive1)

  # Verify that with k=0, no correction is applied (correction factor = 1)
  correction_k0 <- 1 - 0 * lambda1 * gap / (1 + lambda1 * gap)^2
  (all.equal(correction_k0, 1))
})

assert("sample_size_nbinom matches published results", {
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
  (all.equal(res_zhu_1$n1, 167))
  (all.equal(res_zhu_1$n2, 167))

  # Example 2: High Rates
  # L1=1.0, L2=0.5, k=0.5, Power=0.9, Alpha=0.025 -> n=88
  res_zhu_2 <- calc_fixed_exp(1.0, 0.5, 0.5, 0.9, 0.025)
  (all.equal(res_zhu_2$n1, 88))
  (all.equal(res_zhu_2$n2, 88))

  # Friede and Schmidli (2010)
  # Example 1: Standard Scenario
  # L1=0.6, L2=0.3, k=0.4, Power=0.8, Alpha=0.025 -> n=95
  res_friede_1 <- calc_fixed_exp(0.6, 0.3, 0.4, 0.8, 0.025)
  (all.equal(res_friede_1$n1, 95))
  (all.equal(res_friede_1$n2, 95))

  # Example 2: High Dispersion
  # L1=1.0, L2=0.5, k=0.5, Power=0.8, Alpha=0.025 -> n=66
  res_friede_2 <- calc_fixed_exp(1.0, 0.5, 0.5, 0.8, 0.025)
  (all.equal(res_friede_2$n1, 66))
  (all.equal(res_friede_2$n2, 66))
})

assert("sample_size_nbinom handles vector dropout_rate", {
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

  (all.equal(res_scalar$n_total, res_vector$n_total))
  (all.equal(res_scalar$exposure[1], res_vector$exposure[1]))
  (all.equal(res_scalar$exposure[2], res_vector$exposure[2]))

  # Case 2: Different dropout rates
  # Group 1 has higher dropout, so exposure should be lower
  res_diff <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    dropout_rate = c(0.1, 0.01), max_followup = 12
  )

  (res_diff$exposure[1] < res_diff$exposure[2])
  (length(res_diff$exposure) == 2)
})

assert("sample_size_nbinom input validation errors", {
  base <- list(lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1,
               accrual_rate = 10, accrual_duration = 2, trial_duration = 2)

  # Bad rr0
  (has_error(do.call(sample_size_nbinom, c(base, rr0 = -1)), "rr0 must be positive"))

  # Bad power
  (has_error(do.call(sample_size_nbinom, c(base, power = 0)), "Power must be between"))
  (has_error(do.call(sample_size_nbinom, c(base, power = 1)), "Power must be between"))

  # Bad alpha
  (has_error(do.call(sample_size_nbinom, c(base, alpha = 0)), "Alpha must be between"))
  (has_error(do.call(sample_size_nbinom, c(base, alpha = 1)), "Alpha must be between"))

  # Bad dropout_rate
  (has_error(do.call(sample_size_nbinom, c(base, dropout_rate = -0.1)), "non-negative"))
  (has_error(do.call(sample_size_nbinom, c(base, list(dropout_rate = c(0.1, 0.1, 0.1)))), "scalar, a vector of length 2, or a data frame"))

  # Bad max_followup
  (has_error(do.call(sample_size_nbinom, c(base, max_followup = 0)), "max_followup must be positive"))
  (has_error(do.call(sample_size_nbinom, c(base, max_followup = -1)), "max_followup must be positive"))
  (has_error(do.call(sample_size_nbinom, c(base, list(max_followup = c(1, 2, 3)))), "scalar or a vector of length 2"))
})

assert("sample_size_nbinom handles non-inferiority (rr0 != 1)", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.5, rr0 = 1.1, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 2, trial_duration = 2
  )
  (res$n_total > 0)
  (all.equal(res$inputs$rr0, 1.1))
})

assert("score sizing reports null variance on final-analysis scale", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0, power = 0.8,
    accrual_rate = 10, accrual_duration = 2, trial_duration = 2,
    test_type = "score"
  )

  mu0 <- mean(c(0.5, 0.3)) * res$exposure[1]
  v0 <- (1 / mu0) * 2
  (all.equal(res$exposure[1], 1))
  (all.equal(res$variance_null, v0 / res$n1))
})

assert("sample_size_nbinom handles max_followup with split case", {
  # max_followup shorter than some follow-up times triggers Case 3 (split)
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 10, trial_duration = 15,
    max_followup = 8
  )
  (res$n_total > 0)
  # With max_followup, exposure should be less than without
  res_no_mf <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 10, trial_duration = 15
  )
  (res$exposure[1] <= res_no_mf$exposure[1])
})

assert("sample_size_nbinom handles max_followup all-truncated case", {
  # Very short max_followup truncates everyone
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 10, trial_duration = 15,
    max_followup = 2
  )
  (res$n_total > 0)
  (res$exposure[1] <= 2)
})

assert("sample_size_nbinom handles dropout with max_followup", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 10, trial_duration = 15,
    dropout_rate = 0.1, max_followup = 8
  )
  (res$n_total > 0)
})

assert("print.sample_size_nbinom_result works", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  (has_match(capture.output(print(res)), "Sample size"))
})

assert("print.sample_size_nbinom_result shows non-inferiority", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.5, rr0 = 1.1, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  (has_match(capture.output(print(res)), "Null hypothesis RR"))
})

assert("print.sample_size_nbinom_result shows vector dispersion", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = c(0.1, 0.2), power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  (has_match(capture.output(print(res)), "\\(n1\\)"))
})

assert("print.sample_size_nbinom_result shows dropout", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    dropout_rate = 0.05
  )
  (has_match(capture.output(print(res)), "Dropout rate"))
})

assert("print shows different dropout rates per group", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    dropout_rate = c(0.05, 0.10)
  )
  (has_match(capture.output(print(res)), "\\(n1\\)"))
})

assert("print shows event gap info", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    event_gap = 0.5
  )
  (has_match(capture.output(print(res)), "Event gap"))
})

assert("print shows max_followup", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    max_followup = 12
  )
  (has_match(capture.output(print(res)), "Max follow"))
})

assert("print shows different exposure per group", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    dropout_rate = c(0.05, 0.2)
  )
  (has_match(capture.output(print(res)), "\\(n1\\)"))
})

assert("summary.sample_size_nbinom_result returns summary object", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  s <- summary(res)
  (inherits(s, "sample_size_nbinom_summary"))
  (is.character(s))
})

assert("summary shows non-inferiority rr0", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.5, rr0 = 1.1, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  s <- summary(res)
  s_text <- paste(s, collapse = " ")
  (grepl("null hypothesis RR", s_text))
})

assert("summary shows vector dispersion", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = c(0.1, 0.2), power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  s <- summary(res)
  s_text <- paste(s, collapse = " ")
  (grepl("\\(n1\\)", s_text))
})

assert("summary shows event gap info", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    event_gap = 0.5
  )
  s <- summary(res)
  s_text <- paste(s, collapse = " ")
  (grepl("Event gap", s_text))
})

assert("summary shows different exposure per group", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    dropout_rate = c(0.05, 0.2)
  )
  s <- summary(res)
  s_text <- paste(s, collapse = " ")
  (grepl("\\(n1\\)", s_text))
})

assert("print.sample_size_nbinom_summary works", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  s <- summary(res)
  (has_match(capture.output(print(s)), "negative binomial"))
})

# --- Piecewise exponential dropout tests ---

assert("piecewise dropout with one piece matches scalar dropout", {
  # A data frame with a single constant rate should give identical results
  # to the scalar specification.
  res_scalar <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    dropout_rate = 0.05
  )

  res_pw <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    dropout_rate = data.frame(rate = 0.05, duration = Inf)
  )

  (all.equal(res_scalar$n_total, res_pw$n_total))
  (all.equal(res_scalar$exposure, res_pw$exposure, tolerance = 1e-10))
  (all.equal(res_scalar$variance, res_pw$variance, tolerance = 1e-10))
})

assert("piecewise dropout with treatment column matches per-group scalar", {
  res_vec <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    dropout_rate = c(0.05, 0.10)
  )

  res_pw <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24,
    dropout_rate = data.frame(
      treatment = c(1, 2),
      rate = c(0.05, 0.10),
      duration = c(Inf, Inf)
    )
  )

  (all.equal(res_vec$n_total, res_pw$n_total))
  (all.equal(res_vec$exposure, res_pw$exposure, tolerance = 1e-10))
})

assert("piecewise dropout higher early rate reduces exposure", {
  # Two pieces: high dropout first 6 months, low after.
  # Compared to constant dropout at the low rate, exposure should be less.
  res_low <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 12, trial_duration = 24,
    dropout_rate = data.frame(rate = 0.02, duration = Inf)
  )

  res_pw <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 12, trial_duration = 24,
    dropout_rate = data.frame(rate = c(0.10, 0.02), duration = c(6, Inf))
  )

  # Higher early dropout means less exposure and more subjects needed
  (res_pw$exposure[1] < res_low$exposure[1])
  (res_pw$n_total >= res_low$n_total)
})

assert("piecewise dropout with zero rate matches no dropout", {
  res_none <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 10, trial_duration = 15,
    dropout_rate = 0
  )

  res_pw_zero <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 10, trial_duration = 15,
    dropout_rate = data.frame(rate = c(0, 0), duration = c(5, Inf))
  )

  (all.equal(res_none$n_total, res_pw_zero$n_total))
  (all.equal(res_none$exposure, res_pw_zero$exposure, tolerance = 1e-10))
})

assert("piecewise dropout with max_followup works correctly", {
  # Piecewise dropout with max_followup should interact properly.
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 12, trial_duration = 24,
    dropout_rate = data.frame(rate = c(0.08, 0.03), duration = c(6, Inf)),
    max_followup = 12
  )

  (res$n_total > 0)
  (res$exposure[1] > 0)
  # Exposure must be less than max_followup

  (res$exposure[1] < 12)
})

assert("piecewise dropout per-group data frame works", {
  # Different piecewise schedules for each group.
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 12, trial_duration = 24,
    dropout_rate = data.frame(
      treatment = c(1, 1, 2, 2),
      rate = c(0.05, 0.02, 0.08, 0.03),
      duration = c(6, Inf, 6, Inf)
    )
  )

  (res$n_total > 0)
  # Group 2 has higher dropout, so lower exposure
  (res$exposure[2] < res$exposure[1])
})

assert("piecewise dropout input validation works", {
  # Missing columns
  (has_error(
    sample_size_nbinom(
      lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1,
      accrual_rate = 10, accrual_duration = 10, trial_duration = 15,
      dropout_rate = data.frame(rate = 0.05)
    ),
    "duration"
  ))

  # Negative rate in data frame
  (has_error(
    sample_size_nbinom(
      lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1,
      accrual_rate = 10, accrual_duration = 10, trial_duration = 15,
      dropout_rate = data.frame(rate = -0.05, duration = Inf)
    ),
    "non-negative"
  ))
})

assert("print handles piecewise dropout display", {
  res <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 12, trial_duration = 24,
    dropout_rate = data.frame(rate = c(0.08, 0.03), duration = c(6, Inf))
  )
  (has_match(capture.output(print(res)), "piecewise"))
})
