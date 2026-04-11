test_that("nb_sim returns typed object", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = "A", rate = 0.5)
  # Use matching block since fail_rate uses "A"
  sim <- nb_sim(enroll_rate = enroll_rate, fail_rate = fail_rate, max_followup = 1, n = 5, block = "A")
  expect_s3_class(sim, "nb_sim_data")
})

test_that("cut_data_by_date truncates follow-up", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = "A", rate = 0.5)
  dropout_rate <- data.frame(treatment = "A", rate = 0.1, duration = 10)
  # Use matching block
  sim <- nb_sim(enroll_rate, fail_rate, dropout_rate, max_followup = 2, n = 5, block = "A")
  cut <- cut_data_by_date(sim, cut_date = 0.5)
  eligible_ids <- unique(sim$id[sim$enroll_time < 0.5])
  expect_equal(nrow(cut), length(eligible_ids))
  expect_true(all(cut$tte <= 0.5))
  expect_true(all(cut$events >= 0))
})

test_that("cut_data_by_date.default errors on unsupported class", {
  dummy <- data.frame(x = 1)
  class(dummy) <- "unsupported_class"
  expect_error(cut_data_by_date(dummy, cut_date = 1), "No cut_data_by_date")
})

test_that("cut_data_by_date.nb_sim_data errors on invalid cut_date", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = "A", rate = 0.5)
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 1, n = 5, block = "A")

  expect_error(cut_data_by_date(sim, cut_date = NULL), "single finite numeric")
  expect_error(cut_data_by_date(sim, cut_date = Inf), "single finite numeric")
  expect_error(cut_data_by_date(sim, cut_date = c(1, 2)), "single finite numeric")
})

test_that("cut_data_by_date.nb_sim_data returns empty for early cut", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = "A", rate = 0.5)
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 1, n = 5, block = "A")

  # Cut before anyone enrolls
  cut <- cut_data_by_date(sim, cut_date = 0)
  expect_equal(nrow(cut), 0)
  expect_true(all(c("id", "treatment", "enroll_time", "tte", "tte_total", "events") %in% names(cut)))
})

test_that("cut_data_by_date.nb_sim_data works with function gap", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = "A", rate = 0.5)
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 2, n = 10, block = "A")

  # Function-based gap
  cut_func <- cut_data_by_date(sim, cut_date = 1, event_gap = function() 0.1)
  cut_fixed <- cut_data_by_date(sim, cut_date = 1, event_gap = 0.1)

  # Both should produce valid output
  expect_true(nrow(cut_func) > 0)
  expect_true(all(cut_func$events >= 0))
})

test_that("cut_data_by_date.nb_sim_seasonal works", {
  enroll_rate <- data.frame(rate = 20 / (5 / 12), duration = 5 / 12)
  fail_rate <- data.frame(
    treatment = rep(c("Control", "Experimental"), each = 4),
    season = rep(c("Winter", "Spring", "Summer", "Fall"), times = 2),
    rate = c(0.6, 0.5, 0.4, 0.5, 0.4, 0.3, 0.2, 0.3)
  )
  sim <- nb_sim_seasonal(
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    max_followup = 1,
    randomization_start_date = as.Date("2020-01-01"),
    n = 20
  )

  cut <- cut_data_by_date(sim, cut_date = 0.5)
  expect_true(is.data.frame(cut))
  expect_true(all(c("id", "treatment", "season", "tte", "events") %in% names(cut)))
})

test_that("cut_data_by_date.nb_sim_seasonal errors on invalid cut_date", {
  enroll_rate <- data.frame(rate = 20 / (5 / 12), duration = 5 / 12)
  fail_rate <- data.frame(
    treatment = rep(c("Control", "Experimental"), each = 4),
    season = rep(c("Winter", "Spring", "Summer", "Fall"), times = 2),
    rate = rep(0.5, 8)
  )
  sim <- nb_sim_seasonal(
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    max_followup = 1,
    randomization_start_date = as.Date("2020-01-01"),
    n = 10
  )

  expect_error(cut_data_by_date(sim, cut_date = NULL), "single finite numeric")
})

test_that("cut_data_by_date.nb_sim_seasonal returns empty for early cut", {
  enroll_rate <- data.frame(rate = 20 / (5 / 12), duration = 5 / 12)
  fail_rate <- data.frame(
    treatment = rep(c("Control", "Experimental"), each = 4),
    season = rep(c("Winter", "Spring", "Summer", "Fall"), times = 2),
    rate = rep(0.5, 8)
  )
  sim <- nb_sim_seasonal(
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    max_followup = 1,
    randomization_start_date = as.Date("2020-01-01"),
    n = 10
  )

  cut <- cut_data_by_date(sim, cut_date = 0)
  expect_equal(nrow(cut), 0)
  expect_true("season" %in% names(cut))
})

test_that("cut_data_by_date.nb_sim_seasonal handles event_gap", {
  enroll_rate <- data.frame(rate = 20 / (5 / 12), duration = 5 / 12)
  fail_rate <- data.frame(
    treatment = rep(c("Control", "Experimental"), each = 4),
    season = rep(c("Winter", "Spring", "Summer", "Fall"), times = 2),
    rate = c(0.6, 0.5, 0.4, 0.5, 0.4, 0.3, 0.2, 0.3)
  )
  sim <- nb_sim_seasonal(
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    max_followup = 1,
    randomization_start_date = as.Date("2020-01-01"),
    n = 20
  )

  cut_gap <- cut_data_by_date(sim, cut_date = 0.8, event_gap = 0.05)
  cut_no_gap <- cut_data_by_date(sim, cut_date = 0.8, event_gap = 0)

  expect_true(is.data.frame(cut_gap))
  expect_true(nrow(cut_gap) > 0)
})
