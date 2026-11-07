assert("sim_gs_nbinom uses planned information dates without moving the final", {
  fixed <- sample_size_nbinom(
    lambda1 = 0.8, lambda2 = 0.5, dispersion = 0.2, power = NULL,
    accrual_rate = 60, accrual_duration = 1, trial_duration = 4,
    max_followup = 3
  )
  design <- toInteger(gsNB(fixed, k = 3, timing = c(0.4, 0.7)))
  n_target <- utils::tail(design$n_total, 1)
  expected_info <- vapply(design$T, function(time) {
    compute_info_at_time(
      analysis_time = time, accrual_rate = n_target, accrual_duration = 1,
      lambda1 = 0.8, lambda2 = 0.5, dispersion = 0.2, max_followup = 3
    )
  }, numeric(1))
  (all.equal(expected_info, design$n.I, tolerance = 1e-6))
  sim_args <- list(
    n_sims = 2, enroll_rate = data.frame(rate = n_target, duration = 1),
    fail_rate = data.frame(
      treatment = c("Control", "Experimental"),
      rate = c(0.8, 0.5), dispersion = 0.2
    ),
    max_followup = 3, n_target = n_target, design = design, seed = 55
  )
  inherited <- do.call(sim_gs_nbinom, sim_args)
  explicit_args <- sim_args
  explicit_args$analysis_times <- design$T
  explicit <- do.call(sim_gs_nbinom, explicit_args)
  (inherited %==% explicit)
  (all.equal(inherited$analysis_time, rep(design$T, 2)))
  (all.equal(inherited$planned_info, rep(design$n.I, 2)))
  (all.equal(inherited$planned_info_fraction, rep(design$timing, 2)))
  (all.equal(inherited$planned_final_info, rep(design$n.I[3], 6)))
  for (metric in c("unblinded_ml", "blinded_ml", "unblinded_mom", "blinded_mom")) {
    (all.equal(inherited[[paste0("info_fraction_", metric)]],
               inherited[[paste0("info_", metric)]] / design$n.I[3]))
  }

  low_args <- sim_args
  low_args$fail_rate$rate <- c(0.08, 0.05)
  low <- do.call(sim_gs_nbinom, low_args)
  high_args <- sim_args
  high_args$fail_rate$rate <- c(4, 2.5)
  high <- do.call(sim_gs_nbinom, high_args)
  (all.equal(low$analysis_time, inherited$analysis_time))
  (all.equal(high$analysis_time, inherited$analysis_time))
  (low$analysis_time[low$analysis == 3] == 4)
  (high$analysis_time[high$analysis == 3] == 4)
  (low$info_fraction_blinded_ml[low$analysis == 3] < 1)
  (high$info_fraction_blinded_ml[high$analysis == 3] > 1)
})

assert("sim_gs_nbinom runs basic simulation", {
  set.seed(123)
  enroll_rate <- data.frame(rate = 10, duration = 3)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.6, 0.4),
    dispersion = 0.2
  )
  dropout_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.05, 0.05),
    duration = c(6, 6)
  )
  design <- sample_size_nbinom(
    lambda1 = 0.6, lambda2 = 0.4, dispersion = 0.2, power = 0.8,
    accrual_rate = enroll_rate$rate, accrual_duration = enroll_rate$duration,
    trial_duration = 6
  )
  cuts <- list(
    list(planned_calendar = 2),
    list(planned_calendar = 4)
  )
  sim_results <- sim_gs_nbinom(
    n_sims = 2,
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    dropout_rate = dropout_rate,
    max_followup = 4,
    n_target = 30,
    design = design,
    cuts = cuts,
    seed = TRUE
  )

  (is.data.frame(sim_results))

  # Expected columns
  expected_cols <- c("sim", "analysis", "analysis_time", "n_enrolled",
                     "n_ctrl", "n_exp", "events_total", "z_stat",
                     "blinded_info", "unblinded_info",
                     "info_unblinded_ml", "info_blinded_ml",
                     "info_unblinded_mom", "info_blinded_mom")
  (expected_cols %in% names(sim_results))

  # 2 sims x 2 analyses = 4 rows
  (nrow(sim_results) %==% 4L)

  # Sim IDs should be 1 and 2
  (all.equal(sort(unique(sim_results$sim)), c(1, 2)))

  # Analysis IDs should be 1 and 2
  (all.equal(sort(unique(sim_results$analysis)), c(1, 2)))

  # Analysis times should be non-decreasing within each sim
  (tapply(sim_results$analysis_time, sim_results$sim, function(x) all(diff(x) >= 0)))
})

assert("sim_gs_nbinom errors without design", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.5, 0.3),
    dispersion = 0.1
  )

  (has_error(
    sim_gs_nbinom(
      n_sims = 1,
      enroll_rate = enroll_rate,
      fail_rate = fail_rate,
      max_followup = 2,
      cuts = list(list(planned_calendar = 1))
    ),
    "design object must be provided"
  ))
})

assert("sim_gs_nbinom errors without cuts or analysis_times", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.5, 0.3),
    dispersion = 0.1
  )
  design <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
    accrual_rate = 10, accrual_duration = 1, trial_duration = 3
  )

  (has_error(
    sim_gs_nbinom(
      n_sims = 1,
      enroll_rate = enroll_rate,
      fail_rate = fail_rate,
      max_followup = 2,
      design = design
    ),
    "Either analysis_times or cuts must be provided"
  ))
})

assert("sim_gs_nbinom accepts analysis_times", {
  set.seed(42)
  enroll_rate <- data.frame(rate = 10, duration = 3)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.6, 0.4),
    dispersion = 0.2
  )
  design <- sample_size_nbinom(
    lambda1 = 0.6, lambda2 = 0.4, dispersion = 0.2, power = 0.8,
    accrual_rate = 10, accrual_duration = 3, trial_duration = 6
  )

  # Use analysis_times instead of cuts
  sim_results <- sim_gs_nbinom(
    n_sims = 2,
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    max_followup = 4,
    n_target = 30,
    design = design,
    analysis_times = c(2, 4),
    seed = TRUE
  )

  (is.data.frame(sim_results))
  (nrow(sim_results) %==% 4L)
})

assert("sim_gs_nbinom is reproducible with seed", {
  enroll_rate <- data.frame(rate = 10, duration = 2)
  fail_rate <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.6, 0.4),
    dispersion = 0.2
  )
  design <- sample_size_nbinom(
    lambda1 = 0.6, lambda2 = 0.4, dispersion = 0.2, power = 0.8,
    accrual_rate = 10, accrual_duration = 2, trial_duration = 4
  )
  cuts <- list(list(planned_calendar = 2), list(planned_calendar = 4))

  set.seed(42)
  res1 <- sim_gs_nbinom(
    n_sims = 3, enroll_rate = enroll_rate, fail_rate = fail_rate,
    max_followup = 3, n_target = 20, design = design, cuts = cuts, seed = TRUE
  )

  set.seed(42)
  res2 <- sim_gs_nbinom(
    n_sims = 3, enroll_rate = enroll_rate, fail_rate = fail_rate,
    max_followup = 3, n_target = 20, design = design, cuts = cuts, seed = TRUE
  )

  (res1 %==% res2)
})
