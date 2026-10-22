assert("nb_sim returns typed object", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = "A", rate = 0.5)
  # Use matching block since fail_rate uses "A"
  sim <- nb_sim(enroll_rate = enroll_rate, fail_rate = fail_rate, max_followup = 1, n = 5, block = "A")
  (inherits(sim, "nb_sim_data"))
})

assert("cut_data_by_date truncates follow-up", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = "A", rate = 0.5)
  dropout_rate <- data.frame(treatment = "A", rate = 0.1, duration = 10)
  # Use matching block
  sim <- nb_sim(enroll_rate, fail_rate, dropout_rate, max_followup = 2, n = 5, block = "A")
  cut <- cut_data_by_date(sim, cut_date = 0.5)
  eligible_ids <- unique(sim$id[sim$enroll_time < 0.5])
  (nrow(cut) %==% length(eligible_ids))
  (cut$tte <= 0.5)
  (cut$events >= 0)
})

assert("cut_data_by_date.default errors on unsupported class", {
  dummy <- data.frame(x = 1)
  class(dummy) <- "unsupported_class"
  (has_error(cut_data_by_date(dummy, cut_date = 1), "No cut_data_by_date"))
})

assert("cut_data_by_date.nb_sim_data errors on invalid cut_date", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = "A", rate = 0.5)
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 1, n = 5, block = "A")

  (has_error(cut_data_by_date(sim, cut_date = NULL), "single finite numeric"))
  (has_error(cut_data_by_date(sim, cut_date = Inf), "single finite numeric"))
  (has_error(cut_data_by_date(sim, cut_date = c(1, 2)), "single finite numeric"))
})

assert("cut_data_by_date.nb_sim_data returns empty for early cut", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = "A", rate = 0.5)
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 1, n = 5, block = "A")

  # Cut before anyone enrolls
  cut <- cut_data_by_date(sim, cut_date = 0)
  (nrow(cut) %==% 0L)
  (c("id", "treatment", "enroll_time", "tte", "tte_total", "events") %in% names(cut))
})

assert("cut_data_by_date.nb_sim_data works with function gap", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = "A", rate = 0.5)
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 2, n = 10, block = "A")

  # Function-based gap
  cut_func <- cut_data_by_date(sim, cut_date = 1, event_gap = function() 0.1)
  cut_fixed <- cut_data_by_date(sim, cut_date = 1, event_gap = 0.1)

  # Both should produce valid output
  (nrow(cut_func) > 0)
  (cut_func$events >= 0)
})

assert("cut_data_by_date.nb_sim_seasonal works", {
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
  (is.data.frame(cut))
  (c("id", "treatment", "season", "tte", "events") %in% names(cut))
})

assert("cut_data_by_date.nb_sim_seasonal errors on invalid cut_date", {
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

  (has_error(cut_data_by_date(sim, cut_date = NULL), "single finite numeric"))
})

assert("cut_data_by_date.nb_sim_seasonal returns empty for early cut", {
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
  (nrow(cut) %==% 0L)
  ("season" %in% names(cut))
})

assert("cut_data_by_date.nb_sim_seasonal handles event_gap", {
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

  (is.data.frame(cut_gap))
  (nrow(cut_gap) > 0)
})
