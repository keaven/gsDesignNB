test_that("calculate_blinded_info uses subject-level tte (not just an average)", {
  # Access internal helper via ::: (defined inside calculate_blinded_info())
  # To keep the test stable, we test the exported function behavior by
  # crafting data with deterministic tte patterns and large event counts.

  # Two datasets with same mean tte but different distribution
  tte_a <- rep(1, 100)
  tte_b <- c(rep(0.1, 50), rep(1.9, 50))

  # Construct events so glm.nb is well-behaved and lambda_est is similar
  # (events roughly proportional to tte)
  df_a <- data.frame(tte = tte_a, events = round(2 * tte_a + 10))
  df_b <- data.frame(tte = tte_b, events = round(2 * tte_b + 10))

  res_a <- calculate_blinded_info(
    df_a,
    ratio = 1,
    lambda1_planning = 0.5,
    lambda2_planning = 0.3
  )

  res_b <- calculate_blinded_info(
    df_b,
    ratio = 1,
    lambda1_planning = 0.5,
    lambda2_planning = 0.3
  )

  # With the subject-level formula, these should generally differ because
  # the exposure distribution differs, even though mean(tte) is the same.
  expect_true(is.finite(res_a$blinded_info))
  expect_true(is.finite(res_b$blinded_info))
  expect_false(isTRUE(all.equal(res_a$blinded_info, res_b$blinded_info, tolerance = 1e-12)))
})
