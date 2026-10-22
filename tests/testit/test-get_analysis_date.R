# Simple simulation
enroll_rate <- data.frame(rate = 10, duration = 1)
fail_rate <- data.frame(treatment = "A", rate = 5)
sim <- nb_sim(enroll_rate, fail_rate, max_followup = 5, n = 10, block = "A")

# Skip if not enough events were generated
if (sum(sim$event == 1) >= 5) assert("get_analysis_date returns correct time", {
  # Calculate manually
  events_df <- sim[sim$event == 1, ]
  events_df <- events_df[order(events_df$calendar_time), ]
  target_time <- events_df$calendar_time[5]

  res <- get_analysis_date(sim, planned_events = 5, event_gap = 0)
  (all.equal(res, target_time))
})

assert("get_analysis_date handles insufficient events", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = "A", rate = 0.1) # Low rate
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 0.1, n = 10, block = "A")

  # Likely very few events
  total_events <- sum(sim$event == 1)
  target <- total_events + 10

  (has_message(
    res <- get_analysis_date(sim, planned_events = target),
    sprintf("Only %d events in trial", total_events)
  ))

  (all.equal(res, max(sim$calendar_time)))
})
