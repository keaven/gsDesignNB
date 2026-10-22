assert("gsNBCalendar creates valid gsNB object", {
  # Create sample size object
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )

  # Create group sequential design
  gs_design <- gsNBCalendar(nb_ss, k = 3, test.type = 4, analysis_times = c(12, 18, 24))

  # Check class inheritance
  (inherits(gs_design, "gsNB"))
  (inherits(gs_design, "gsDesign"))
  (inherits(gs_design, "sample_size_nbinom_result"))

  # Check that nb_design is preserved
  (gs_design$nb_design %==% nb_ss)

  # Check sample size vectors
  (length(gs_design$n1) == 3)
  (length(gs_design$n2) == 3)
  (length(gs_design$n_total) == 3)

  # Cumulative sample sizes should increase
  (diff(gs_design$n1) > 0)
  (diff(gs_design$n2) > 0)

  # Final sample sizes should match original (approximately, due to inflation)
  (all.equal(
    gs_design$n1[3] + gs_design$n2[3],
    gs_design$n_total[3]
  ))
})

assert("gsNBCalendar rejects invalid input", {
  # Not a sample_size_nbinom_result object
  (has_error(
    gsNBCalendar(list(n_total = 100)),
    "must be an object of class 'sample_size_nbinom_result'"
  ))
})

assert("gsNBCalendar respects allocation ratio", {
  # Create sample size with 2:1 allocation
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    ratio = 2,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )

  gs_design <- gsNBCalendar(nb_ss, k = 2, analysis_times = c(12, 24))

  # Check ratio is preserved
  (all.equal(gs_design$n2[2] / gs_design$n1[2], 2))
})

assert("gsNBCalendar works with different test types", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )

  # One-sided
  gs1 <- gsNBCalendar(nb_ss, k = 2, test.type = 1, analysis_times = c(12, 24))
  (inherits(gs1, "gsNB"))

  # Two-sided symmetric
  gs2 <- gsNBCalendar(nb_ss, k = 2, test.type = 2, analysis_times = c(12, 24))
  (inherits(gs2, "gsNB"))

  # Two-sided asymmetric non-binding
  gs4 <- gsNBCalendar(nb_ss, k = 2, test.type = 4, analysis_times = c(12, 24))
  (inherits(gs4, "gsNB"))
})

assert("gsNBCalendar supports harm-bound designs from gsDesign", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )

  gs_design <- gsNBCalendar(
    nb_ss,
    k = 3,
    test.type = 8,
    astar = 0.025,
    analysis_times = c(10, 18, 24),
    testHarm = c(FALSE, TRUE, TRUE)
  )

  (inherits(gs_design, "gsNB"))
  (!is.null(gs_design$harm))
  (gs_design$testHarm %==% c(FALSE, TRUE, TRUE))
  (gs_design$harm$bound[1] < 0)
  (all.equal(gs_design$harm$spend[1], 0))
  (all.equal(as.numeric(gs_design$harm$prob[1, ]), rep(0, length(gs_design$theta))))
  (is.finite(gs_design$harm$bound[-1]))
  ("Harm" %in% names(gsDesign::gsBoundSummary(gs_design)))
})

assert("gsNBCalendar works with custom spending functions", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )

  # O'Brien-Fleming-like spending
  gs_design <- gsNBCalendar(
    nb_ss,
    k = 3,
    sfu = gsDesign::sfHSD,
    sfupar = -4,
    sfl = gsDesign::sfHSD,
    sflpar = -2,
    analysis_times = c(12, 18, 24)
  )

  (inherits(gs_design, "gsNB"))
})

assert("summary.gsNB returns gsNBsummary object", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  gs_design <- gsNBCalendar(nb_ss, k = 3, analysis_times = c(12, 18, 24))

  s <- summary(gs_design)
  (inherits(s, "gsNBsummary"))
  (is.character(s))

  s_collapsed <- paste(s, collapse = " ")
  (nchar(s_collapsed) > 0)

  # Should mention key terms
  (grepl("negative binomial", s_collapsed, ignore.case = TRUE))
  (grepl("0.5", s_collapsed))
  (grepl("0.3", s_collapsed))
})

assert("print.gsNBsummary prints without error", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  gs_design <- gsNBCalendar(nb_ss, k = 2, analysis_times = c(12, 24))
  s <- summary(gs_design)

  (length(capture.output(print(s))) > 0)
})

assert("toInteger.gsNB rounds final sample size", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  gs_design <- gsNBCalendar(nb_ss, k = 2, analysis_times = c(12, 24))

  gs_int <- toInteger(gs_design)

  (inherits(gs_int, "gsNB"))

  # Final sample size should be integer (within rounding)
  k <- gs_int$k
  (all.equal(gs_int$n_total[k], round(gs_int$n_total[k])))

  # n1 + n2 should equal n_total
  (all.equal(gs_int$n1[k] + gs_int$n2[k], gs_int$n_total[k]))
})

assert("toInteger.gsNB respects ratio", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    ratio = 2,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )
  gs_design <- gsNBCalendar(nb_ss, k = 2, analysis_times = c(12, 24))

  gs_int <- toInteger(gs_design)

  # Final total should be divisible by (ratio + 1) = 3
  (all.equal(gs_int$n_total[2] %% 3, 0))
})

assert("toInteger.gsNB preserves calendar enrollment at interim analyses", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 12, trial_duration = 24
  )
  analysis_times <- c(10, 18, 24)
  gs_design <- gsNBCalendar(nb_ss, k = 3, analysis_times = analysis_times)

  gs_int <- toInteger(gs_design)

  expected_n <- gs_int$n_total[gs_int$k] *
    pmin(analysis_times, nb_ss$accrual_duration) / nb_ss$accrual_duration

  (all.equal(gs_int$n_total, expected_n, tolerance = 1e-8))
  (diff(gs_int$n.I) > 0)
})

assert("toInteger.gsNB preserves piecewise accrual shape", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.35, dispersion = 0.5, power = 0.9,
    accrual_rate = c(12, 24), accrual_duration = c(6, 6),
    trial_duration = 24, dropout_rate = 0.10 / 12,
    max_followup = 12, event_gap = 20 / 30.4375,
    test_type = "score"
  )
  analysis_times <- c(9, 14, 24)
  gs_design <- gsNBCalendar(nb_ss, k = 3, analysis_times = analysis_times)

  gs_int <- toInteger(gs_design)

  information_scale <- gs_design$n.I[gs_design$k] / (1 / nb_ss$variance)
  rounding_scale <- gs_int$n_total[gs_int$k] / gs_design$n_total[gs_design$k]
  expected_n <- vapply(
    analysis_times,
    function(t) {
      sample_size_nbinom(
        lambda1 = nb_ss$inputs$lambda1,
        lambda2 = nb_ss$inputs$lambda2,
        rr0 = nb_ss$inputs$rr0,
        dispersion = nb_ss$inputs$dispersion,
        power = NULL,
        alpha = nb_ss$inputs$alpha,
        sided = nb_ss$inputs$sided,
        ratio = nb_ss$inputs$ratio,
        accrual_rate = nb_ss$accrual_rate * information_scale * rounding_scale,
        accrual_duration = nb_ss$accrual_duration,
        trial_duration = t,
        dropout_rate = nb_ss$inputs$dropout_rate,
        max_followup = nb_ss$inputs$max_followup,
        event_gap = nb_ss$inputs$event_gap
      )$n_total
    },
    numeric(1)
  )

  (all.equal(gs_int$n_total, expected_n, tolerance = 1e-8))
  (gs_int$n_total[gs_int$k] < 1.05 * gs_design$n_total[gs_design$k])
})

assert("toInteger.gsDesign dispatches correctly", {
  gs <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 2)
  gs_int <- toInteger(gs)
  (inherits(gs_int, "gsDesign"))
})

assert("update_gsNB preserves selective harm-bound settings", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )

  gs_design <- gsNBCalendar(
    nb_ss,
    k = 3,
    test.type = 8,
    astar = 0.025,
    analysis_times = c(10, 18, 24),
    testUpper = c(FALSE, TRUE, TRUE),
    testLower = c(TRUE, TRUE, FALSE),
    testHarm = c(FALSE, TRUE, TRUE)
  )

  original <- gs_design
  updated <- update_gsNB(gs_design, observed_info = gs_design$n.I[1])

  (gs_design %==% original)
  (updated$design$testUpper %==% gs_design$testUpper)
  (updated$design$testLower %==% gs_design$testLower)
  (updated$design$testHarm %==% gs_design$testHarm)
  (updated$design$upper$bound[1] %==% gs_design$upper$bound[1])
  (updated$design$harm$bound[1] %==% gs_design$harm$bound[1])
  (updated$design$lower$bound[3] %==% gs_design$lower$bound[3])
  (all.equal(updated$design$upper$spend[1], 0))
  (all.equal(updated$design$harm$spend[1], 0))
  (inherits(updated$bounds, "gsBoundSummary"))
  (all.equal(updated$bounds, gsDesign::gsBoundSummary(updated$design)))
  (names(updated$bounds) %==% c("Analysis", "Value", "Harm", "Futility", "Efficacy"))
  z_rows <- which(updated$bounds$Value == "Z")
  (length(z_rows) == gs_design$k)
  (is.na(updated$bounds$Harm[z_rows[1]:(z_rows[2] - 1)]))
  (is.na(updated$bounds$Efficacy[z_rows[1]:(z_rows[2] - 1)]))
  (all.equal(
    updated$bounds$Harm[z_rows[-1]],
    round(updated$design$harm$bound[-1], 4)
  ))
})

assert("update_gsNB retains the usual summary for active bounds", {
  for (test_type in c(2, 4, 7, 8)) {
    design <- gsDesign::gsDesign(
      k = 3, test.type = test_type, astar = 0.025,
      n.fix = 100, timing = c(0.35, 0.7, 1)
    )
    updated <- update_gsNB(design, observed_info = 0.35 * design$n.fix)

    (all.equal(updated$bounds, gsDesign::gsBoundSummary(design)))
  }
})

assert("update_gsNB summaries preserve harm designs with absent bounds", {
  for (test_type in c(7, 8)) {
    for (zero_spend in c(FALSE, TRUE)) {
      args <- list(
        k = 3, test.type = test_type, astar = 0.025,
        n.fix = 100, timing = c(0.35, 0.7, 1)
      )
      if (zero_spend) {
        # Zero spending also creates an absent bound when testHarm is TRUE.
        args$sfharm <- gsDesign::sfPoints
        args$sfharmparam <- c(0, 0.5, 1)
        inactive <- 1L
      } else {
        args$testHarm <- c(FALSE, FALSE, TRUE)
        inactive <- 1:2
      }
      design <- do.call(gsDesign::gsDesign, args)
      original <- design
      updated <- update_gsNB(design, observed_info = 0.35 * design$n.fix)

      (design %==% original)
      components <- c("n.I", "theta", "upper", "lower", "harm", "testHarm")
      (all.equal(updated$design[components], design[components]))

      bounds <- updated$bounds
      (inherits(bounds, "gsBoundSummary"))
      (!any(bounds$Value %in% c("CP", "CP H1", "PP")))
      analysis <- cumsum(bounds$Value == "Z")
      (is.na(bounds$Harm[analysis %in% inactive]))
      (is.finite(bounds$Harm[!analysis %in% inactive]))
      active <- setdiff(seq_len(design$k), inactive)
      (all.equal(
        bounds$Harm[bounds$Value == "Z"][-inactive],
        round(design$harm$bound[active], 4)
      ))
      (all.equal(
        bounds$Harm[bounds$Value == "p (1-sided)"][-inactive],
        round(pnorm(design$harm$bound[active], lower.tail = FALSE), 4)
      ))
      (all.equal(
        bounds$Harm[bounds$Value == "~delta at bound"][-inactive],
        round(gsDesign::gsDelta(design, z = design$harm$bound, i = 1:3)[active], 4)
      ))
      for (j in seq_along(design$theta)) {
        cross_rows <- bounds$Value == paste0("P(Cross) if delta=", design$theta[j] / design$delta)
        for (bound in c("upper", "lower", "harm")) {
          column <- c(upper = "Efficacy", lower = "Futility", harm = "Harm")[[bound]]
          expected <- round(cumsum(design[[bound]]$prob[, j]), 4)
          if (bound == "harm") expected[inactive] <- NA_real_
          (all.equal(bounds[[column]][cross_rows], expected))
        }
      }
    }
  }
})

assert("compute_info_at_time returns positive value", {
  info <- compute_info_at_time(
    analysis_time = 12,
    accrual_rate = 10,
    accrual_duration = 10,
    lambda1 = 0.5,
    lambda2 = 0.3,
    dispersion = 0.1
  )

  (is.numeric(info))
  (info > 0)
})

assert("compute_info_at_time increases with analysis time", {
  base_args <- list(
    accrual_rate = 10,
    accrual_duration = 10,
    lambda1 = 0.5,
    lambda2 = 0.3,
    dispersion = 0.1
  )

  info_early <- do.call(compute_info_at_time, c(list(analysis_time = 6), base_args))
  info_late <- do.call(compute_info_at_time, c(list(analysis_time = 12), base_args))

  (info_late > info_early)
})

assert("gsNBCalendar errors without analysis_times", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )

  (has_error(
    gsNBCalendar(nb_ss, k = 3),
    "analysis_times must be provided"
  ))
})

assert("gsNBCalendar errors when analysis_times length != k", {
  nb_ss <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
    accrual_rate = 10, accrual_duration = 20, trial_duration = 24
  )

  (has_error(
    gsNBCalendar(nb_ss, k = 3, analysis_times = c(12, 24)),
    "analysis_times must have length k"
  ))
})
