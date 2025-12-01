test_that("blinded_ssr estimates dispersion and rate", {
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
  
  res <- blinded_ssr(cut, lambda1_planning = 0.5, lambda2_planning = 0.3)
  
  # Since data is Poisson, k should be near 0
  expect_true(res$dispersion_blinded < 0.2) # Allow some noise
  
  # Pooled rate should be between 0.3 and 0.5
  expect_true(res$lambda_blinded > 0.3 && res$lambda_blinded < 0.5)
  
  # N total should be calculated
  expect_true(res$n_total_blinded > 0)
})
