assert("blinded_ssr estimates dispersion and rate", {
  # Simulate data with known k and rate
  # Control: 0.5, Exp: 0.3. Ratio 1:1.
  # Pooled rate approx 0.4. k = 0.5.

  enroll_rate <- data.frame(rate = 100, duration = 1)
  fail_rate <- data.frame(treatment = c("C", "E"), rate = c(0.5, 0.3))
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 10, n = 200, block = c("C", "E"))

  # Add dispersion (not directly supported in nb_sim which is Poisson,
  # but we can hack it or just test that the function runs on Poisson data where k=0)
  # nb_sim generates Poisson counts (k=0).

  cut <- cut_data_by_date(sim, cut_date = 5)

  # Mock accrual params for blinded_ssr (needed for sample size recalc)
  res <- blinded_ssr(
    cut,
    lambda1_planning = 0.5,
    lambda2_planning = 0.3,
    accrual_rate = 100,
    accrual_duration = 1,
    trial_duration = 10
  )

  # Since data is Poisson, k should be near 0, but mixture effect might inflate it slightly
  # Allow more noise
  (res$dispersion_blinded < 0.5)

  # Pooled rate should be between 0.3 and 0.5
  (res$lambda_blinded > 0.3 && res$lambda_blinded < 0.5)

  # N total should be calculated
  (res$n_total_blinded > 0)

  # Check information output
  (c("blinded_info", "target_info", "info_fraction") %in% names(res))
  (res$blinded_info > 0)
  (res$info_fraction > 0)
})

assert("blinded SSR projects fixed-size final information with consistent event gaps", {
  set.seed(40)
  sim <- nb_sim(
    data.frame(rate = 40, duration = 2),
    data.frame(treatment = c("Control", "Experimental"),
               rate = c(0.8, 0.5), dispersion = 0.4),
    max_followup = 4, n = 80, event_gap = 0.2
  )
  cut <- cut_data_by_date(sim, cut_date = 3, event_gap = 0.2)
  blinded <- calculate_blinded_info(cut, lambda1_planning = 0.8,
                                    lambda2_planning = 0.5, event_gap = 0.2)
  result <- blinded_ssr(
    cut, lambda1_planning = 0.8, lambda2_planning = 0.5,
    accrual_rate = 40, accrual_duration = 2, trial_duration = 6,
    max_followup = 4, event_gap = 0.2
  )
  expected <- sample_size_nbinom(
    lambda1 = blinded$lambda1_adjusted, lambda2 = blinded$lambda2_adjusted,
    dispersion = blinded$dispersion_blinded, power = NULL,
    accrual_rate = 40, accrual_duration = 2, trial_duration = 6,
    max_followup = 4, event_gap = 0.2
  )
  (all.equal(result$blinded_info, blinded$blinded_info))
  (all.equal(result$info_projected, 1 / expected$variance))
  (result$fallback %==% blinded$fallback)
  planned <- list(lambda1_adjusted = 0.8, lambda2_adjusted = 0.5,
                  dispersion_blinded = 0.4)
  projected <- .gs_nb_projected_info(planned, 40, 2, 6, 0, 4, 0.2, 1,
                                    analysis_time = 1, n_observed = 40)
  reference <- sample_size_nbinom(
    lambda1 = 0.8, lambda2 = 0.5, dispersion = 0.4, power = NULL,
    accrual_rate = 40, accrual_duration = 2, trial_duration = 6,
    max_followup = 4, event_gap = 0.2
  )
  (all.equal(projected, 1 / reference$variance))
})
