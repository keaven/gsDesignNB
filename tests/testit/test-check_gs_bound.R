assert("new monitoring strategies share expected boundaries with gsNBPower", {
  fixed <- sample_size_nbinom(
    lambda1 = 0.5, lambda2 = 0.35, dispersion = 0.5, power = 0.9,
    accrual_rate = 1, accrual_duration = 12, trial_duration = 24, max_followup = 12
  )
  design <- toInteger(gsNB(fixed, test.type = 1))
  for (strategy in c("min_planned_actual", "planned_final", "expected_final", "calendar")) {
    expected <- gsNBPower(design, spending = strategy)
    rows <- data.frame(sim = 1, analysis = 1:3, z_stat = 0,
                       analysis_time = design$T, info_blinded_ml = design$n.I,
                       info_blinded_final_projected = design$n.I[3])
    checked <- check_gs_bound(rows, design, info_col = "info_blinded_ml", spending = strategy)
    (checked$bound_status %==% rep("ok", 3))
    (all.equal(checked$upper_bound, expected$look_results$upper_bound, tolerance = 1e-6))
    (checked$spending_time[3] %==% 1)
  }
  rows$z_stat <- c(-5, 0, 0)
  checked <- check_gs_bound(rows, design, info_col = "info_blinded_ml", boundary_update = "fixed")
  (checked$cross_upper[1])
  (all.equal(checked$upper_bound[1], design$upper$bound[1]))
  rows$info_blinded_ml <- NA_real_
  checked <- check_gs_bound(rows, design, info_col = "info_blinded_ml", boundary_update = "fixed")
  (checked$cross_upper[1])
  rows$info_blinded_ml <- design$n.I
  rows$z_stat <- 0
  rows$status <- c("skipped", "skipped", "analyzed")
  checked <- check_gs_bound(rows, design, info_col = "info_blinded_ml", spending = "expected_final")
  (checked$bound_status %==% c("skipped", "skipped", "ok"))
  rows$status <- "analyzed"
  rows$info_blinded_ml[2] <- rows$info_blinded_ml[1] / 2
  checked <- check_gs_bound(rows, design, info_col = "info_blinded_ml", spending = "planned_final")
  (checked$bound_status[2] %==% "failed")
  rows$info_blinded_ml <- design$n.I
  rows$status <- "analyzed"
  baseline <- check_gs_bound(rows, design, info_col = "info_blinded_ml", spending = "min_planned_actual")
  rows$status[2] <- "skipped"
  altered <- check_gs_bound(rows, design, info_col = "info_blinded_ml", spending = "min_planned_actual")
  (all.equal(altered$upper_bound[1], baseline$upper_bound[1]))
})

assert("check_gs_bound returns correct structure", {
  design <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 2, timing = c(0.5, 1))
  sim_df <- data.frame(
    sim = c(1, 1, 2, 2),
    analysis = c(1, 2, 1, 2),
    z_stat = c(2.5, NA, -0.2, 2.2),
    blinded_info = c(50, 100, 50, 100),
    unblinded_info = c(50, 100, 50, 100)
  )

  result <- check_gs_bound(sim_df, design)

  (is.data.frame(result))
  (c("cross_upper", "cross_lower", "cross_harm") %in% names(result))

  # Same number of rows
 (nrow(result) %==% nrow(sim_df))
})

assert("check_gs_bound detects upper crossings", {
  design <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 2, timing = c(0.5, 1))
  # Large z should cross upper bound
  sim_df <- data.frame(
    sim = c(1, 1),
    analysis = c(1, 2),
    z_stat = c(5.0, NA),
    blinded_info = c(50, 100),
    unblinded_info = c(50, 100)
  )

  result <- check_gs_bound(sim_df, design)
  (result$cross_upper[1])
})

assert("check_gs_bound supports efficacy-only designs", {
  design <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 1,
                               timing = c(0.5, 1))
  sim_df <- data.frame(
    sim = c(1, 1, 2, 2),
    analysis = c(1, 2, 1, 2),
    z_stat = c(5, NA, 0.5, 2.5),
    blinded_info = c(50, 100, 50, 100),
    unblinded_info = c(50, 100, 50, 100)
  )

  result <- check_gs_bound(sim_df, design)
  (result$cross_upper[result$sim == 1 & result$analysis == 1])
  (result$cross_upper[result$sim == 2 & result$analysis == 2])
  (!any(result$cross_lower))
})

assert("check_gs_bound detects lower crossings for test.type >= 2", {
  design <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 4, timing = c(0.5, 1))
  # Very negative z should cross lower bound
  sim_df <- data.frame(
    sim = c(1, 1),
    analysis = c(1, 2),
    z_stat = c(-5.0, NA),
    blinded_info = c(50, 100),
    unblinded_info = c(50, 100)
  )

  result <- check_gs_bound(sim_df, design)
  (result$cross_lower[1])
})

assert("check_gs_bound detects harm crossings for harm-bound designs", {
  design <- gsDesign::gsDesign(
    k = 2,
    n.fix = 100,
    test.type = 8,
    astar = 0.025,
    timing = c(0.5, 1)
  )
  sim_df <- data.frame(
    sim = c(1, 1),
    analysis = c(1, 2),
    z_stat = c(-5.0, NA),
    blinded_info = c(50, 100),
    unblinded_info = c(50, 100)
  )

  result <- check_gs_bound(sim_df, design)
  (result$cross_harm[1])
  (!result$cross_lower[1])
})

assert("check_gs_bound errors on bad design", {
  sim_df <- data.frame(
    sim = 1, analysis = 1, z_stat = 1.5,
    blinded_info = 50, unblinded_info = 50
  )
  (has_error(check_gs_bound(sim_df, list(n.fix = 100)), "gsDesign or gsNB"))
})

assert("check_gs_bound errors on missing z_stat", {
  design <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 1)
  sim_df <- data.frame(sim = c(1, 1), analysis = c(1, 2), blinded_info = c(50, 100), unblinded_info = c(50, 100))
  (has_error(check_gs_bound(sim_df, design), "z_stat"))
})

assert("check_gs_bound errors on missing info column", {
  design <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 1)
  sim_df <- data.frame(sim = c(1, 1), analysis = c(1, 2), z_stat = c(2.0, 1.0))
  (has_error(check_gs_bound(sim_df, design), "blinded_info"))
})

assert("check_gs_bound accepts custom info_col", {
  design <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 2, timing = c(0.5, 1))
  sim_df <- data.frame(
    sim = c(1, 1),
    analysis = c(1, 2),
    z_stat = c(1.0, 2.0),
    blinded_info = c(50, 100),
    unblinded_info = c(50, 100)
  )

  result <- check_gs_bound(sim_df, design, info_col = "unblinded_info")
  (is.data.frame(result))
  ("cross_upper" %in% names(result))
})

assert("check_gs_bound handles NA z_stat gracefully", {
  design <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 2, timing = c(0.5, 1))
  sim_df <- data.frame(
    sim = c(1, 1),
    analysis = c(1, 2),
    z_stat = c(NA_real_, NA_real_),
    blinded_info = c(50, 100),
    unblinded_info = c(50, 100)
  )

  result <- check_gs_bound(sim_df, design)
  (!result$cross_upper)
  (!result$cross_lower)
})

assert("check_gs_bound works with multiple simulations", {
  design <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 2, timing = c(0.5, 1))
  sim_df <- data.frame(
    sim = c(1, 1, 2, 2, 3, 3),
    analysis = c(1, 2, 1, 2, 1, 2),
    z_stat = c(5.0, NA, 0.5, 1.0, -5.0, NA),
    blinded_info = c(50, 100, 50, 100, 50, 100),
    unblinded_info = c(50, 100, 50, 100, 50, 100)
  )

  result <- check_gs_bound(sim_df, design)
  (nrow(result) %==% 6L)

  # Sim 1: crosses upper at analysis 1
  (result$cross_upper[result$sim == 1 & result$analysis == 1])
  # Sim 3: crosses lower at analysis 1
  (result$cross_lower[result$sim == 3 & result$analysis == 1])
})
