# Common setup for sim_gs_nbinom tests
setup_sim_args <- function() {
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
  list(
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    dropout_rate = dropout_rate,
    design = design,
    cuts = cuts
  )
}

test_that("sim_gs_nbinom is reproducible with seed = TRUE (sequential)", {
  skip_if_not_installed("future.apply")

  args <- setup_sim_args()

  set.seed(42)
  res1 <- sim_gs_nbinom(
    n_sims = 5,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    n_target = 30,
    design = args$design,
    cuts = args$cuts,
    seed = TRUE
  )

  set.seed(42)
  res2 <- sim_gs_nbinom(
    n_sims = 5,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    n_target = 30,
    design = args$design,
    cuts = args$cuts,
    seed = TRUE
  )

  expect_identical(res1, res2)
})

test_that("sim_gs_nbinom is reproducible with integer seed", {
  skip_if_not_installed("future.apply")

  args <- setup_sim_args()

  res1 <- sim_gs_nbinom(
    n_sims = 5,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    n_target = 30,
    design = args$design,
    cuts = args$cuts,
    seed = 12345
  )

  res2 <- sim_gs_nbinom(
    n_sims = 5,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    n_target = 30,
    design = args$design,
    cuts = args$cuts,
    seed = 12345
  )

  expect_identical(res1, res2)
})

test_that("sim_gs_nbinom returns correct structure with seed", {
  skip_if_not_installed("future.apply")

  args <- setup_sim_args()
  n_sims <- 3

  set.seed(99)
  res <- sim_gs_nbinom(
    n_sims = n_sims,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    n_target = 30,
    design = args$design,
    cuts = args$cuts,
    seed = TRUE
  )

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), n_sims * length(args$cuts))
  expect_true(all(c("sim", "analysis", "z_stat", "estimate") %in% names(res)))
  expect_equal(sort(unique(res$sim)), seq_len(n_sims))
})

test_that("sim_gs_nbinom parallel results match sequential results", {
  skip_if_not_installed("future.apply")
  skip_if_not_installed("future")
  skip_on_cran()

  args <- setup_sim_args()

  # Run sequentially
  future::plan(future::sequential)
  set.seed(42)
  res_seq <- sim_gs_nbinom(
    n_sims = 10,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    n_target = 30,
    design = args$design,
    cuts = args$cuts,
    seed = TRUE
  )

  # Run with multisession (2 workers)
  future::plan(future::multisession, workers = 2)
  set.seed(42)
  res_par <- sim_gs_nbinom(
    n_sims = 10,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    n_target = 30,
    design = args$design,
    cuts = args$cuts,
    seed = TRUE
  )
  future::plan(future::sequential)

  # Integer/count columns should be identical (same RNG streams)
  int_cols <- c("sim", "analysis", "n_enrolled", "n_ctrl", "n_exp",
                "events_total", "events_ctrl", "events_exp")
  expect_identical(res_seq[, int_cols], res_par[, int_cols])

  # Floating-point columns match within tolerance (iterative optimization

  # in MASS::glm.nb can produce tiny differences across processes)
  num_cols <- c("z_stat", "estimate", "se", "blinded_info", "unblinded_info",
                "info_unblinded_ml", "info_blinded_ml",
                "info_unblinded_mom", "info_blinded_mom")
  for (col in num_cols) {
    expect_equal(res_seq[[col]], res_par[[col]], tolerance = 1e-3,
                 label = paste("Column", col))
  }
})

test_that("sim_gs_nbinom different seeds give different results", {
  skip_if_not_installed("future.apply")

  args <- setup_sim_args()

  res1 <- sim_gs_nbinom(
    n_sims = 5,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    n_target = 30,
    design = args$design,
    cuts = args$cuts,
    seed = 111
  )

  res2 <- sim_gs_nbinom(
    n_sims = 5,
    enroll_rate = args$enroll_rate,
    fail_rate = args$fail_rate,
    dropout_rate = args$dropout_rate,
    max_followup = 4,
    n_target = 30,
    design = args$design,
    cuts = args$cuts,
    seed = 222
  )

  expect_false(identical(res1$z_stat, res2$z_stat))
})
