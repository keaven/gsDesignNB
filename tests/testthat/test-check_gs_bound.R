test_that("check_gs_bound returns correct structure", {
  design <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 2, timing = c(0.5, 1))
  sim_df <- data.frame(
    sim = c(1, 1, 2, 2),
    analysis = c(1, 2, 1, 2),
    z_stat = c(2.5, NA, -0.2, 2.2),
    blinded_info = c(50, 100, 50, 100),
    unblinded_info = c(50, 100, 50, 100)
  )

  result <- check_gs_bound(sim_df, design)

  expect_true(is.data.frame(result))
  expect_true("cross_upper" %in% names(result))
  expect_true("cross_lower" %in% names(result))
  expect_true("cross_harm" %in% names(result))

  # Same number of rows
 expect_equal(nrow(result), nrow(sim_df))
})

test_that("check_gs_bound detects upper crossings", {
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
  expect_true(result$cross_upper[1])
})

test_that("check_gs_bound detects lower crossings for test.type >= 2", {
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
  expect_true(result$cross_lower[1])
})

test_that("check_gs_bound detects harm crossings for harm-bound designs", {
  skip_if_not("testHarm" %in% names(formals(gsDesign::gsDesign)))

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
  expect_true(result$cross_harm[1])
  expect_false(result$cross_lower[1])
})

test_that("check_gs_bound errors on bad design", {
  sim_df <- data.frame(
    sim = 1, analysis = 1, z_stat = 1.5,
    blinded_info = 50, unblinded_info = 50
  )
  expect_error(check_gs_bound(sim_df, list(n.fix = 100)), "gsDesign or gsNB")
})

test_that("check_gs_bound errors on missing z_stat", {
  design <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 1)
  sim_df <- data.frame(sim = c(1, 1), analysis = c(1, 2), blinded_info = c(50, 100), unblinded_info = c(50, 100))
  expect_error(check_gs_bound(sim_df, design), "z_stat")
})

test_that("check_gs_bound errors on missing info column", {
  design <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 1)
  sim_df <- data.frame(sim = c(1, 1), analysis = c(1, 2), z_stat = c(2.0, 1.0))
  expect_error(check_gs_bound(sim_df, design), "blinded_info")
})

test_that("check_gs_bound accepts custom info_col", {
  design <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 2, timing = c(0.5, 1))
  sim_df <- data.frame(
    sim = c(1, 1),
    analysis = c(1, 2),
    z_stat = c(1.0, 2.0),
    blinded_info = c(50, 100),
    unblinded_info = c(50, 100)
  )

  result <- check_gs_bound(sim_df, design, info_col = "unblinded_info")
  expect_true(is.data.frame(result))
  expect_true("cross_upper" %in% names(result))
})

test_that("check_gs_bound handles NA z_stat gracefully", {
  design <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 2, timing = c(0.5, 1))
  sim_df <- data.frame(
    sim = c(1, 1),
    analysis = c(1, 2),
    z_stat = c(NA_real_, NA_real_),
    blinded_info = c(50, 100),
    unblinded_info = c(50, 100)
  )

  result <- check_gs_bound(sim_df, design)
  expect_false(any(result$cross_upper))
  expect_false(any(result$cross_lower))
})

test_that("check_gs_bound works with multiple simulations", {
  design <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 2, timing = c(0.5, 1))
  sim_df <- data.frame(
    sim = c(1, 1, 2, 2, 3, 3),
    analysis = c(1, 2, 1, 2, 1, 2),
    z_stat = c(5.0, NA, 0.5, 1.0, -5.0, NA),
    blinded_info = c(50, 100, 50, 100, 50, 100),
    unblinded_info = c(50, 100, 50, 100, 50, 100)
  )

  result <- check_gs_bound(sim_df, design)
  expect_equal(nrow(result), 6)

  # Sim 1: crosses upper at analysis 1
  expect_true(result$cross_upper[result$sim == 1 & result$analysis == 1])
  # Sim 3: crosses lower at analysis 1
  expect_true(result$cross_lower[result$sim == 3 & result$analysis == 1])
})
