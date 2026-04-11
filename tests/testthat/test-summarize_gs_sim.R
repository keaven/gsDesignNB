test_that("summarize_gs_sim returns correct structure", {
  # Create mock data that check_gs_bound would produce
  sim_df <- data.frame(
    sim = c(1, 1, 2, 2, 3, 3),
    analysis = c(1, 2, 1, 2, 1, 2),
    z_stat = c(2.5, NA, -0.5, 2.1, 1.0, 1.5),
    blinded_info = c(40, 80, 40, 80, 40, 80),
    unblinded_info = c(40, 80, 40, 80, 40, 80),
    n_enrolled = c(30, 60, 30, 60, 30, 60),
    events_total = c(12, 25, 10, 22, 11, 20),
    cross_upper = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
    cross_lower = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
  )

  result <- summarize_gs_sim(sim_df)

  expect_true(is.list(result))
  expect_true(all(c("n_sim", "power", "futility", "analysis_summary") %in% names(result)))

  # 3 simulations
 expect_equal(result$n_sim, 3)

  # Power: sims 1 and 2 cross upper
  expect_equal(result$power, 2 / 3)

  # Futility: sim 3 crosses lower but not upper
  expect_equal(result$futility, 1 / 3)

  # Analysis summary should be a data frame with 2 rows
  expect_true(is.data.frame(result$analysis_summary))
  expect_equal(nrow(result$analysis_summary), 2)

  # Check cumulative probability column exists
  expect_true("cum_prob_upper" %in% names(result$analysis_summary))
})

test_that("summarize_gs_sim errors without required columns", {
  bad_df <- data.frame(sim = 1, analysis = 1, z_stat = 1.5)

  expect_error(
    summarize_gs_sim(bad_df),
    "cross_upper.*cross_lower"
  )
})

test_that("summarize_gs_sim handles all successes", {
  sim_df <- data.frame(
    sim = c(1, 1, 2, 2),
    analysis = c(1, 2, 1, 2),
    z_stat = c(3.0, NA, 3.5, NA),
    blinded_info = c(40, 80, 40, 80),
    unblinded_info = c(40, 80, 40, 80),
    n_enrolled = c(30, 60, 30, 60),
    events_total = c(15, 25, 15, 25),
    cross_upper = c(TRUE, FALSE, TRUE, FALSE),
    cross_lower = c(FALSE, FALSE, FALSE, FALSE)
  )

  result <- summarize_gs_sim(sim_df)
  expect_equal(result$power, 1.0)
  expect_equal(result$futility, 0.0)
})

test_that("summarize_gs_sim handles no crossings", {
  sim_df <- data.frame(
    sim = c(1, 1),
    analysis = c(1, 2),
    z_stat = c(0.5, 1.0),
    blinded_info = c(40, 80),
    unblinded_info = c(40, 80),
    n_enrolled = c(30, 60),
    events_total = c(10, 20),
    cross_upper = c(FALSE, FALSE),
    cross_lower = c(FALSE, FALSE)
  )

  result <- summarize_gs_sim(sim_df)
  expect_equal(result$power, 0.0)
  expect_equal(result$futility, 0.0)
})
