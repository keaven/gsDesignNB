assert("nb_sim generates correct structure", {
  enroll_rate <- data.frame(rate = c(10, 10), duration = c(1, 1))
  fail_rate <- data.frame(treatment = c("A", "B"), rate = c(0.5, 0.8))
  dropout_rate <- data.frame(treatment = c("A", "B"), rate = c(0.1, 0.1), duration = c(100, 100))
  max_followup <- 5

  # Use block matching treatments "A", "B"
  sim <- nb_sim(enroll_rate, fail_rate, dropout_rate, max_followup, n = 100, block = c("A", "B"))

  (inherits(sim, "data.frame"))
  (c("id", "treatment", "enroll_time", "tte", "event", "calendar_time") %in% names(sim))
  (length(unique(sim$id)) %==% 100L)

  # Check max followup
  # tte should never exceed max_followup (approx, allow for floating point)
  (sim$tte <= max_followup + 1e-5)

  # Check censoring rows exist
  censored <- sim[sim$event == 0, ]
  (nrow(censored) %==% 100L) # Each subject should have exactly one censoring row

  # Check calendar_time calculation
  (all.equal(sim$calendar_time, sim$enroll_time + sim$tte))
})

assert("nb_sim handles dropout", {
  # High dropout rate -> shorter follow-up
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = "A", rate = 0.1)
  dropout_rate <- data.frame(treatment = "A", rate = 100, duration = 10) # Instant dropout
  max_followup <- 10

  sim <- nb_sim(enroll_rate, fail_rate, dropout_rate, max_followup, n = 50, block = "A")

  # Follow-up times should be very short
  followup_times <- sim$tte[sim$event == 0]
  (mean(followup_times) < 1)
})

assert("nb_sim respects n argument", {
  enroll_rate <- data.frame(rate = 100, duration = 1) # implies 100
  fail_rate <- data.frame(treatment = "A", rate = 0.1)
  max_followup <- 1

  # Force n=50
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = max_followup, n = 50, block = "A")
  (length(unique(sim$id)) %==% 50L)
})
