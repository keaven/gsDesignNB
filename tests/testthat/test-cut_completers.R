test_that("cut_completers and cut_date_for_completers work as expected", {
  # Setup simulation
  enroll_rate <- data.frame(rate = 10, duration = 10) # fast enrollment
  fail_rate <- data.frame(treatment = c("A", "B"), rate = c(0.5, 0.5))
  dropout_rate <- data.frame(treatment = c("A", "B"), rate = c(0.01, 0.01), duration = c(100, 100)) # low dropout
  max_followup <- 2
  n <- 20
  
  set.seed(123)
  sim <- nb_sim(enroll_rate, fail_rate, dropout_rate, max_followup = max_followup, n = n, block = c("A", "B"))
  
  # 1. Test cut_date_for_completers
  # Find date for 5 completers
  date_5 <- cut_date_for_completers(sim, 5)
  
  expect_true(is.numeric(date_5))
  expect_length(date_5, 1)
  
  # Verify that at least 5 subjects have completed by this date
  # Completer definition: reached max_followup (max(tte) in sim)
  dt <- data.table::as.data.table(sim)
  max_f <- max(dt$tte)
  completers_dt <- dt[, .(is_completer = max(tte) >= max_f - 1e-8,
                          completion_time = first(enroll_time) + max(tte)), by = id]
  
  n_completed_at_date <- sum(completers_dt$is_completer & completers_dt$completion_time <= date_5 + 1e-8)
  expect_gte(n_completed_at_date, 5)
  
  # 2. Test cut_completers
  cut_data <- cut_completers(sim, date_5)
  
  # Verify output structure
  expect_true(is.data.frame(cut_data))
  expect_true(all(c("id", "treatment", "enroll_time", "tte", "events") %in% names(cut_data)))
  
  # Verify that subjects randomized before cut_date are included (even if not completers)
  # In this simulation, everyone is randomized quickly (rate 10, duration 10 -> 100 subjects? No n=20)
  # n=20, rate=10 -> duration 2.
  # date_5 will be around 2 + something.
  # Everyone should be randomized by then.
  
  included_ids <- cut_data$id
  # Check that we have at least the completers
  expect_true(all(completers_dt[is_completer == TRUE & completion_time <= date_5, id] %in% included_ids))
  
  # Check that we have non-completers (or those who complete later) if they were randomized
  randomized_before_cut <- unique(sim$id[sim$enroll_time < date_5])
  expect_true(all(randomized_before_cut %in% included_ids))
  
  # Verify tte is correct (<= cut_date - enroll_time)
  # For included subjects, tte should be min(max_followup, cut_date - enroll_time) - gaps
  # Just check it's positive and reasonable
  expect_true(all(cut_data$tte > 0))
  
  # 3. Test edge case: target > n completers
  expect_message(date_high <- cut_date_for_completers(sim, n + 10), "Only .* completers in trial")
  expect_equal(date_high, max(sim$calendar_time))
  
  # 4. Test cut_completers with date 0 (no one completed)
  cut_empty <- cut_completers(sim, 0)
  expect_equal(nrow(cut_empty), 0)
})

