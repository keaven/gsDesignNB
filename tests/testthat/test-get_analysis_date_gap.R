test_that("get_analysis_date respects event gaps", {
  # ID 1: Events at 0.1, 0.15, 0.3. Gap 0.1.
  # Valid: 0.1 (gap ends 0.2). 0.15 skipped. 0.3 valid.
  # Total valid: 2.
  # ID 2: Events at 0.2. Gap 0.1.
  # Valid: 0.2.
  # Total valid: 1.
  # Total in study: 3.

  # Calendar times:
  # ID 1 enroll 0. Events cal: 0.1, 0.15, 0.3.
  # ID 2 enroll 0. Events cal: 0.2.

  # Valid calendar times: 0.1, 0.2, 0.3.

  # Target: 3 events. Date should be 0.3.

  d <- data.frame(
    id = c(1, 1, 1, 2),
    treatment = "A",
    enroll_time = 0,
    tte = c(0.1, 0.15, 0.3, 0.2),
    event = c(1, 1, 1, 1),
    calendar_time = c(0.1, 0.15, 0.3, 0.2)
  )
  class(d) <- c("nb_sim_data", "data.frame")

  res <- get_analysis_date(d, planned_events = 3, event_gap = 0.1)
  expect_equal(res, 0.3)

  # Target: 2 events. Sorted valid: 0.1, 0.2. Date should be 0.2.
  res2 <- get_analysis_date(d, planned_events = 2, event_gap = 0.1)
  expect_equal(res2, 0.2)

  # With gap = 0, all 4 events valid.
  # Sorted: 0.1, 0.15, 0.2, 0.3.
  # Target 3 -> 0.2.
  res3 <- get_analysis_date(d, planned_events = 3, event_gap = 0)
  expect_equal(res3, 0.2)
})
