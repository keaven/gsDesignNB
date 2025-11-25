test_that("nb_sim returns typed object", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = "A", rate = 0.5)
  sim <- nb_sim(enroll_rate = enroll_rate, fail_rate = fail_rate, max_followup = 1, n = 5)
  expect_s3_class(sim, "nb_sim_data")
})

test_that("cut_data_by_date truncates follow-up", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = "A", rate = 0.5)
  dropout_rate <- data.frame(treatment = "A", rate = 0.1, duration = 10)
  sim <- nb_sim(enroll_rate, fail_rate, dropout_rate, max_followup = 2, n = 5)
  cut <- cut_data_by_date(sim, cut_date = 0.5)
  eligible_ids <- unique(sim$id[sim$enroll_time < 0.5])
  expect_equal(nrow(cut), length(eligible_ids))
  expect_true(all(cut$tte <= 0.5))
  expect_true(all(cut$events >= 0))
})
