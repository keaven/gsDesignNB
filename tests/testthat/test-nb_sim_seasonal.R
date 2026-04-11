test_that("nb_sim_seasonal returns correct structure", {
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

  expect_s3_class(sim, "nb_sim_seasonal")
  expect_s3_class(sim, "data.frame")

  # Check required columns
  expected_cols <- c("id", "treatment", "season", "enroll_time",
                     "start", "end", "event")
  expect_true(all(expected_cols %in% names(sim)))

  # All events should be 0 or 1

  expect_true(all(sim$event %in% c(0, 1)))

  # All treatments should be from assigned set
  expect_true(all(sim$treatment %in% c("Control", "Experimental")))

  # All seasons should be valid
  expect_true(all(sim$season %in% c("Winter", "Spring", "Summer", "Fall")))

  # Start should be <= end for each row
  expect_true(all(sim$start <= sim$end + 1e-10))
})

test_that("nb_sim_seasonal errors without start date", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(
    treatment = rep(c("Control", "Experimental"), each = 4),
    season = rep(c("Winter", "Spring", "Summer", "Fall"), times = 2),
    rate = rep(0.5, 8)
  )
  expect_error(
    nb_sim_seasonal(enroll_rate, fail_rate, max_followup = 1),
    "randomization_start_date is required"
  )
})

test_that("nb_sim_seasonal handles dispersion", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(
    treatment = rep(c("Control", "Experimental"), each = 4),
    season = rep(c("Winter", "Spring", "Summer", "Fall"), times = 2),
    rate = c(0.6, 0.5, 0.4, 0.5, 0.4, 0.3, 0.2, 0.3),
    dispersion = 0.5
  )
  sim <- nb_sim_seasonal(
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    max_followup = 1,
    randomization_start_date = as.Date("2020-06-01"),
    n = 30
  )

  expect_s3_class(sim, "nb_sim_seasonal")
  # Should have data for all subjects
  expect_true(length(unique(sim$id)) == 30)
})

test_that("nb_sim_seasonal validates inputs", {
  fail_rate <- data.frame(
    treatment = rep(c("Control", "Experimental"), each = 4),
    season = rep(c("Winter", "Spring", "Summer", "Fall"), times = 2),
    rate = rep(0.5, 8)
  )

  expect_error(
    nb_sim_seasonal(
      enroll_rate = NULL, fail_rate = fail_rate,
      max_followup = 1, randomization_start_date = as.Date("2020-01-01")
    ),
    "enroll_rate must be a data frame"
  )

  expect_error(
    nb_sim_seasonal(
      enroll_rate = data.frame(rate = 10, duration = 1),
      fail_rate = NULL,
      max_followup = 1, randomization_start_date = as.Date("2020-01-01")
    ),
    "fail_rate must be a data frame"
  )

  expect_error(
    nb_sim_seasonal(
      enroll_rate = data.frame(rate = 10, duration = 1),
      fail_rate = fail_rate,
      max_followup = NULL, randomization_start_date = as.Date("2020-01-01")
    ),
    "max_followup must be provided"
  )
})

test_that("nb_sim_seasonal handles dropout", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(
    treatment = rep(c("Control", "Experimental"), each = 4),
    season = rep(c("Winter", "Spring", "Summer", "Fall"), times = 2),
    rate = rep(0.5, 8)
  )
  dropout_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.1, 0.1),
    duration = c(1, 1)
  )
  sim <- nb_sim_seasonal(
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    dropout_rate = dropout_rate,
    max_followup = 1,
    randomization_start_date = as.Date("2020-01-01"),
    n = 30
  )

  expect_s3_class(sim, "nb_sim_seasonal")
  expect_true(nrow(sim) > 0)
})
