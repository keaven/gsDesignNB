assert("one-look score and Wald power match fixed-size formulas", {
  for (test_type in c("score", "wald")) {
    args <- list(lambda1 = 0.5, lambda2 = 0.35, dispersion = 0.5,
                 accrual_rate = 20, accrual_duration = 12, trial_duration = 24,
                 max_followup = 12, test_type = test_type)
    fixed <- do.call(sample_size_nbinom, c(args, list(power = NULL)))
    result <- do.call(gsNBPower, args)
    (all.equal(result$power, fixed$power, tolerance = 1e-8))
    (result$expected_n %==% fixed$n_total)
    (result$expected_analysis_time %==% 24)
  }
})

assert("power uses signed effects and preserves the reference design", {
  fixed <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.35, dispersion = 0.5, power = 0.9,
    accrual_rate = 1, accrual_duration = 12, trial_duration = 24, max_followup = 12
  )
  design <- toInteger(gsNB(fixed, test.type = 1))
  original <- design
  result <- gsNBPower(design, boundary_update = "fixed")
  reference <- gsDesign::gsProbability(
    k = 3, theta = log(0.5 / 0.35), n.I = design$n.I,
    a = rep(-Inf, 3), b = design$upper$bound
  )
  (all.equal(result$power, sum(reference$upper$prob), tolerance = 1e-6))
  (all.equal(result$look_results$information, design$n.I, tolerance = 1e-6))
  (design %==% original)
  (gsNBPower(design, lambda2 = 0.7)$power < gsNBPower(design, lambda2 = 0.5)$power)
})

assert("recomputed boundaries retain fixed information and full final spending", {
  fixed <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.35, dispersion = 0.5, power = 0.9,
    accrual_rate = 1, accrual_duration = 12, trial_duration = 24, max_followup = 12
  )
  design <- toInteger(gsNB(fixed))
  result <- gsNBPower(design)
  (all.equal(result$look_results$upper_bound, design$upper$bound, tolerance = 1e-6))
  (all.equal(result$look_results$lower_bound, design$lower$bound, tolerance = 1e-6))
  (all.equal(result$look_results$information, design$n.I, tolerance = 1e-6))
  (result$look_results$spending_time[3] %==% 1)
  (result$expected_n <= design$n_total[3])
  (result$power + result$futility_probability + result$harm_probability <= 1 + 1e-7)
})

assert("expected cutoff requirements combine and unreachable interims are skipped", {
  fixed <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.35, dispersion = 0.5, power = 0.9,
    accrual_rate = 1, accrual_duration = 12, trial_duration = 24, max_followup = 12
  )
  design <- toInteger(gsNB(fixed, test.type = 1))
  expected <- gsNBPower(
    design, planned_calendar = c(NA, NA, 24),
    target_info = c(design$n.I[1:2], NA)
  )
  (all.equal(expected$look_results$analysis_time, design$T, tolerance = 1e-6))
  delayed <- gsNBPower(design, min_n = c(design$n_total[3], NA, NA),
                        min_followup = c(1, NA, NA))
  (delayed$look_results$analysis_time[1] >= 13)
  skipped <- gsNBPower(design, target_info = c(1e9, NA, NA))
  (skipped$look_results$status %==% c("skipped", "analyzed", "analyzed"))
  (skipped$look_results$analysis_time[3] %==% 24)
  (skipped$look_results$efficacy_prob[1] %==% 0)
  final_only <- gsNBPower(design, target_info = c(1e9, 1e9, NA))
  (final_only$look_results$status %==% c("skipped", "skipped", "analyzed"))
  (is.finite(final_only$power))
  (has_error(gsNBPower(design, target_info = Inf), "finite"))
  (has_error(gsNBPower(design, min_followup = 1), "requires min_n"))
})

assert("multi-look score probabilities match canonical simulation", {
  set.seed(40)
  n_sims <- 60000L
  info <- c(20, 40, 60)
  info_null <- c(25, 55, 70)
  theta <- 0.2
  upper <- c(2.8, 2.4, 2)
  expected <- .gs_nb_probabilities(info, info_null, theta, upper,
                                    rep(-Inf, 3), rep(-Inf, 3))
  increments <- sweep(matrix(rnorm(n_sims * 3), nrow = n_sims), 2,
                       sqrt(c(info[1], diff(info))), "*")
  cumulative <- increments
  cumulative[, 2] <- increments[, 1] + increments[, 2]
  cumulative[, 3] <- rowSums(increments)
  estimates <- sweep(cumulative, 2, info, "/") + theta
  score <- sweep(estimates, 2, sqrt(info_null), "*")
  rejected <- rowSums(sweep(score, 2, upper, ">")) > 0
  empirical <- mean(rejected)
  mcse <- sqrt(empirical * (1 - empirical) / n_sims)
  (abs(empirical - sum(expected$upper)) < 5 * mcse)
})

assert("all spending rules complete nominal null efficacy spending", {
  fixed <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.35, dispersion = 0.5, power = 0.9,
    accrual_rate = 1, accrual_duration = 12, trial_duration = 24, max_followup = 12
  )
  design <- toInteger(gsNB(fixed, test.type = 1))
  for (rate in c(0.15, 1)) {
    for (strategy in c("min_planned_actual", "planned_final", "expected_final", "calendar")) {
      result <- gsNBPower(design, lambda1 = rate, lambda2 = rate, spending = strategy)
      (all.equal(result$power, design$alpha, tolerance = 1e-5))
      (result$look_results$spending_time[3] %==% 1)
    }
  }
})

assert("legacy serialized designs retain usable spending functions", {
  cache <- readRDS(system.file("extdata", "gs_simulation_results.rds", package = "gsDesignNB"))
  result <- gsNBPower(cache$gs_nb)
  (is.finite(result$power))
  (all.equal(result$look_results$information, cache$gs_nb$n.I, tolerance = 1e-6))
})