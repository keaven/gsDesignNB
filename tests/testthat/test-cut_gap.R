test_that("cut_data_by_date applies fixed gap correctly", {
  # ID 1: Enroll 0. Events: 0.1, 0.15, 0.3, 0.4. Cut 1.0. Gap 0.1.
  # Events expected:
  # - 0.1 (valid). Gap [0.1, 0.2].
  # - 0.15 (in gap). Ignored.
  # - 0.3 (valid, > 0.2). Gap [0.3, 0.4].
  # - 0.4 (valid, >= 0.4). Gap [0.4, 0.5].
  # Total events: 3.
  # Total gap deducted:
  # - Gap 1: 0.1 (ends 0.2 <= 1.0)
  # - Gap 2: 0.1 (ends 0.4 <= 1.0)
  # - Gap 3: 0.1 (ends 0.5 <= 1.0)
  # Total gap: 0.3.
  # Time at risk: 1.0 - 0.3 = 0.7.

  # ID 2: Enroll 0. Events: 0.95. Cut 1.0. Gap 0.1.
  # Event 0.95 (valid). Gap [0.95, 1.05].
  # Gap used: min(1.05, 1.0) - 0.95 = 0.05.
  # Time at risk: 1.0 - 0.05 = 0.95.

  d <- data.frame(
    id = c(1, 1, 1, 1, 1, 2, 2),
    treatment = "A",
    enroll_time = 0,
    tte = c(0.1, 0.15, 0.3, 0.4, 1.0, 0.95, 1.0),
    event = c(1, 1, 1, 1, 0, 1, 0)
  )
  d$calendar_time <- d$enroll_time + d$tte
  class(d) <- c("nb_sim_data", "data.frame")

  res <- cut_data_by_date(d, cut_date = 1.0, event_gap = 0.1)

  # Check ID 1
  r1 <- res[res$id == 1, ]
  expect_equal(r1$events, 3)
  expect_equal(r1$tte, 0.7)

  # Check ID 2
  r2 <- res[res$id == 2, ]
  expect_equal(r2$events, 1)
  expect_equal(r2$tte, 0.95)
})

test_that("cut_data_by_date applies random gap correctly", {
  # ID 1: Enroll 0. Events: 0.1, 0.15. Cut 1.0.
  # Gap function: returns 0.2 always (mocked).
  # Event 0.1. Gap [0.1, 0.3].
  # Event 0.15 ignored.
  # Events: 1.
  # Gap deducted: 0.2.
  # TTE: 0.8.

  d <- data.frame(
    id = c(1, 1, 1),
    treatment = "A",
    enroll_time = 0,
    tte = c(0.1, 0.15, 1.0),
    event = c(1, 1, 0)
  )
  d$calendar_time <- d$enroll_time + d$tte
  class(d) <- c("nb_sim_data", "data.frame")

  gap_fun <- function() 0.2

  res <- cut_data_by_date(d, cut_date = 1.0, event_gap = gap_fun)

  expect_equal(res$events, 1)
  expect_equal(res$tte, 0.8)
})

test_that("cut_data_by_date uses default gap", {
  # Default gap is 5/365.25 ~ 0.01369
  # Event at 0.1. Next event at 0.11 (diff 0.01). Should be ignored.

  gap <- 5 / 365.25

  d <- data.frame(
    id = c(1, 1, 1),
    treatment = "A",
    enroll_time = 0,
    tte = c(0.1, 0.1 + gap / 2, 1.0),
    event = c(1, 1, 0)
  )
  d$calendar_time <- d$enroll_time + d$tte
  class(d) <- c("nb_sim_data", "data.frame")

  res <- cut_data_by_date(d, cut_date = 1.0)

  expect_equal(res$events, 1)
  expect_equal(res$tte, 1.0 - gap)
})
