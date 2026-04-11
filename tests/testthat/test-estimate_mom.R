test_that("estimate_nb_mom works for blinded data", {
  # Synthetic data: Poisson (k=0)
  # Rate = 0.5
  set.seed(123)
  n <- 500
  tte <- runif(n, 0.5, 2)
  # Poisson means variance = mean
  events <- rpois(n, lambda = 0.5 * tte)
  
  df <- data.frame(events = events, tte = tte)
  
  res <- estimate_nb_mom(df)
  
  expect_type(res, "list")
  expect_true(res$dispersion >= 0)
  # For Poisson, dispersion should be close to 0
  expect_lt(res$dispersion, 0.2) 
  # Rate should be close to 0.5
  expect_equal(res$lambda, sum(events)/sum(tte))
  
  # Synthetic data: NB (k=1)
  # Rate = 0.5, k = 1 => Variance = mu + mu^2
  # rnegbin(mu, theta) where theta = 1/k
  events_nb <- MASS::rnegbin(n, mu = 0.5 * tte, theta = 1) 
  df_nb <- data.frame(events = events_nb, tte = tte)
  
  res_nb <- estimate_nb_mom(df_nb)
  expect_type(res_nb, "list")
  # Expect k close to 1
  expect_true(abs(res_nb$dispersion - 1) < 0.3) 
})

test_that("estimate_nb_mom works for unblinded data", {
  set.seed(456)
  n <- 200
  group <- rep(c("A", "B"), each = 100)
  tte <- runif(n, 0.5, 2)
  
  # Rates: A=0.2, B=0.8. Common k=0.5
  lambda_vec <- ifelse(group == "A", 0.2, 0.8)
  events <- MASS::rnegbin(n, mu = lambda_vec * tte, theta = 1/0.5)
  
  df <- data.frame(events = events, tte = tte, arm = group)
  
  res <- estimate_nb_mom(df, group = "arm")
  
  expect_type(res, "list")
  expect_named(res$lambda, c("A", "B"))
  
  # Check rates
  expect_equal(unname(res$lambda["A"]), sum(events[group=="A"])/sum(tte[group=="A"]))
  expect_equal(unname(res$lambda["B"]), sum(events[group=="B"])/sum(tte[group=="B"]))
  
  # Check dispersion (should be close to 0.5)
  expect_true(abs(res$dispersion - 0.5) < 0.3)
})

test_that("estimate_nb_mom handles edge cases", {
  # Zero events
  df <- data.frame(events = c(0, 0), tte = c(1, 1))
  res <- estimate_nb_mom(df)
  expect_equal(res$lambda, 0)
  expect_equal(res$dispersion, 0)
  
  # Zero tte (should be filtered)
  df2 <- data.frame(events = c(1, 2), tte = c(0, 1))
  # Only row 2 is used
  res2 <- estimate_nb_mom(df2)
  expect_equal(res2$lambda, 2/1)
  
  # Filtered to empty
  df3 <- data.frame(events = c(1), tte = c(0))
  # Should warn
  expect_warning(res3 <- estimate_nb_mom(df3), "No data")
  expect_true(is.na(res3$lambda))
})

test_that("estimate_nb_mom errors on bad input", {
  expect_error(estimate_nb_mom("not a data frame"), "must be a data frame")
  expect_error(estimate_nb_mom(data.frame(x = 1)), "events.*tte")
})

test_that("estimate_nb_mom errors on missing group column", {
  df <- data.frame(events = c(1, 2), tte = c(1, 1))
  expect_error(estimate_nb_mom(df, group = "nonexistent"), "not found")
})
