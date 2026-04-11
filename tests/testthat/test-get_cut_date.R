test_that("get_cut_date returns planned_calendar when no other criteria", {
  set.seed(456)
  enroll_rate <- data.frame(rate = 15, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.6, 0.4)
  )
  sim_data <- nb_sim(enroll_rate, fail_rate, max_followup = 2, n = 30)

  result <- get_cut_date(sim_data, planned_calendar = 0.5, event_gap = 0)
  expect_equal(result, 0.5)
})

test_that("get_cut_date returns max of multiple criteria", {
  set.seed(789)
  enroll_rate <- data.frame(rate = 30, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(1.0, 0.8)
  )
  sim_data <- nb_sim(enroll_rate, fail_rate, max_followup = 3, n = 60)

  # Request a calendar time AND a target number of events

  # The result should be the max of the two dates
  result <- get_cut_date(
    sim_data, planned_calendar = 0.5,
    target_events = 100,
    event_gap = 0
  )

  # Should be >= the calendar time criterion
  expect_true(result >= 0.5)
})

test_that("get_cut_date respects min_date and max_date", {
  set.seed(101)
  enroll_rate <- data.frame(rate = 20, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.6, 0.4)
  )
  sim_data <- nb_sim(enroll_rate, fail_rate, max_followup = 2, n = 40)

  # With min_date > planned_calendar
  result <- get_cut_date(sim_data, planned_calendar = 0.3, min_date = 0.5, event_gap = 0)
  expect_true(result >= 0.5)

  # With max_date limiting the result
  result2 <- get_cut_date(sim_data, planned_calendar = 5.0, max_date = 2.0, event_gap = 0)
  expect_true(result2 <= 2.0)
})

test_that("get_cut_date returns max_date when no criteria", {
  set.seed(202)
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.5, 0.3)
  )
  sim_data <- nb_sim(enroll_rate, fail_rate, max_followup = 2, n = 20)

  # No criteria provided
  result <- get_cut_date(sim_data, max_date = 3.0)
  expect_equal(result, 3.0)
})

test_that("get_cut_date errors when lambda missing for target_info", {
  set.seed(303)
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.5, 0.3)
  )
  sim_data <- nb_sim(enroll_rate, fail_rate, max_followup = 2, n = 20)

  expect_error(
    get_cut_date(sim_data, target_info = 50),
    "lambda1 and lambda2 must be provided"
  )
})

test_that("get_cut_date handles target_info criterion", {
  set.seed(404)
  enroll_rate <- data.frame(rate = 30, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.8, 0.5)
  )
  sim_data <- nb_sim(enroll_rate, fail_rate, max_followup = 3, n = 60)

  # Use a modest target info so uniroot can find it
  result <- get_cut_date(
    sim_data,
    target_info = 5,
    lambda1 = 0.8,
    lambda2 = 0.5,
    event_gap = 0,
    min_date = 0.1,
    max_date = 3.0
  )

  expect_true(is.numeric(result))
  expect_true(result >= 0.1)
  expect_true(result <= 3.0)
})

test_that("get_cut_date handles target_events", {
  set.seed(505)
  enroll_rate <- data.frame(rate = 20, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.6, 0.4)
  )
  sim_data <- nb_sim(enroll_rate, fail_rate, max_followup = 3, n = 40)

  result <- get_cut_date(sim_data, target_events = 10, event_gap = 0)
  expect_true(is.numeric(result))
  expect_true(result > 0)
})

test_that("get_cut_date handles target_completers", {
  set.seed(606)
  enroll_rate <- data.frame(rate = 20, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.6, 0.4)
  )
  sim_data <- nb_sim(enroll_rate, fail_rate, max_followup = 2, n = 40)

  result <- get_cut_date(sim_data, target_completers = 5)
  expect_true(is.numeric(result))
  expect_true(result > 0)
})

test_that("get_cut_date target_info with already-met info at min_date", {
  set.seed(707)
  enroll_rate <- data.frame(rate = 30, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.8, 0.5)
  )
  sim_data <- nb_sim(enroll_rate, fail_rate, max_followup = 3, n = 60)

  # Very low target info, should be met immediately
  result <- get_cut_date(
    sim_data,
    target_info = 0.001,
    lambda1 = 0.8,
    lambda2 = 0.5,
    event_gap = 0,
    min_date = 0.5,
    max_date = 3.0
  )

  expect_true(result <= 3.0)
})

test_that("get_cut_date target_info unreachable returns max_date", {
  set.seed(808)
  enroll_rate <- data.frame(rate = 5, duration = 0.5)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.5, 0.3)
  )
  sim_data <- nb_sim(enroll_rate, fail_rate, max_followup = 1, n = 5)

  # Extremely high target info, unreachable
  result <- get_cut_date(
    sim_data,
    target_info = 100000,
    lambda1 = 0.5,
    lambda2 = 0.3,
    event_gap = 0,
    min_date = 0.1,
    max_date = 1.5
  )

  expect_equal(result, 1.5)
})
