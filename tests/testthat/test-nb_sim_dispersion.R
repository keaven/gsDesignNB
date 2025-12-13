test_that("nb_sim generates Negative Binomial data with correct dispersion", {
  # This test verifies that the variance of the simulated counts matches
  # the theoretical Negative Binomial variance: Var = mean + k * mean^2

  set.seed(123)
  n_sim <- 5000 # Large sample for stable variance estimation

  enroll_rate <- data.frame(rate = 1000, duration = 1)
  # Use a simple setup: 1 year follow-up, mean rate 2, dispersion 0.5
  # Theoretical Mean = 2
  # Theoretical Var = 2 + 0.5 * 2^2 = 4
  target_mean <- 2
  target_k <- 0.5

  fail_rate <- data.frame(
    treatment = "A",
    rate = target_mean,
    dispersion = target_k
  )

  sim <- nb_sim(
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    max_followup = 1,
    n = n_sim,
    block = "A"
  )

  # Aggregate counts
  counts <- as.data.frame(table(factor(sim$id, levels = 1:n_sim)))
  # Note: table() counts rows. Each event is a row. Censoring is a row with event=0.
  # But nb_sim returns multiple rows: events + 1 censoring.
  # We want to sum(event) per id.

  library(data.table)
  dt <- as.data.table(sim)
  per_subj <- dt[, .(events = sum(event)), by = id]

  observed_mean <- mean(per_subj$events)
  observed_var <- var(per_subj$events)

  # Method of Moments estimator for k: k = (Var - Mean) / Mean^2
  estimated_k <- (observed_var - observed_mean) / (observed_mean^2)

  # Check Mean (Standard Error ~ sqrt(Var/n) = sqrt(4/5000) approx 0.028)
  expect_equal(observed_mean, target_mean, tolerance = 0.1)

  # Check Dispersion
  # The variance of the variance estimator is high, so we need a loose tolerance or larger N.
  # With N=5000, it should be reasonably close.
  expect_equal(estimated_k, target_k, tolerance = 0.1)
})
