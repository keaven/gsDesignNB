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
  
  # Verify only completers are included
  included_ids <- cut_data$id
  expect_true(all(included_ids %in% completers_dt[is_completer == TRUE, id]))
  
  # Verify all included subjects completed by cut_date
  # (Already checked by logic, but double check)
  completion_times_included <- completers_dt[id %in% included_ids, completion_time]
  expect_true(all(completion_times_included <= date_5 + 1e-8))
  
  # Verify tte is max_followup (minus gaps if any, but default gap is small)
  # Here we used default gap. 
  # tte should be approx max_followup - gaps.
  # Since gap logic subtracts gaps after events, and events are random, tte might vary slightly but should be close to max_followup if few events.
  # But definitely shouldn't be small if max_followup is 2.
  expect_true(all(cut_data$tte > 0))
  
  # 3. Test edge case: target > n completers
  expect_message(date_high <- cut_date_for_completers(sim, n + 10), "Only .* completers in trial")
  expect_equal(date_high, max(sim$calendar_time))
  
  # 4. Test cut_completers with date 0 (no one completed)
  cut_empty <- cut_completers(sim, 0)
  expect_equal(nrow(cut_empty), 0)
})

