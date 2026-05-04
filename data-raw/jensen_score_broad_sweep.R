# Broad Jensen correction sweep with Wald and score analyses
#
# This extends data-raw/jensen_correction_broad_sweep.R by analyzing each
# corrected-vs-naive paired data set with both the Wald and score tests and by
# adding matched null simulations for empirical Type I error.
#
# Output:
#   inst/extdata/jensen_score_broad_sweep_results.rds
#
# Usage:
#   Rscript data-raw/jensen_score_broad_sweep.R
#   Rscript data-raw/jensen_score_broad_sweep.R --debug

library(parallel)
devtools::load_all(".")

args <- commandArgs(trailingOnly = TRUE)
debug_mode <- "--debug" %in% args

n_cores <- as.integer(Sys.getenv("GSDESIGNNB_N_CORES", unset = "0"))
if (!is.finite(n_cores) || n_cores < 1L) {
  detected_cores <- detectCores()
  if (!is.finite(detected_cores) || detected_cores < 2L) {
    detected_cores <- 2L
  }
  n_cores <- max(1L, min(6L, detected_cores - 1L))
}
cat(sprintf("Using %d cores for parallel simulation\n", n_cores))

set.seed(20260502)

# ---- Fixed parameters ----
rr <- 0.7
alpha <- 0.025
power_target <- 0.9
dropout <- 0.1 / 12
max_fu <- 12
trial_dur <- 24
accrual_rel <- c(1, 2)
accrual_dur <- c(6, 6)

if (debug_mode) {
  nsim_power <- 40L
  nsim_null <- 80L
  lambda1_grid <- c(0.5, 1.0)
  k_grid <- c(0.5, 1.0)
  gap_days_grid <- c(20, 45)
} else {
  nsim_power <- as.integer(Sys.getenv("GSDESIGNNB_JENSEN_POWER_NSIM", unset = "3500"))
  nsim_null <- as.integer(Sys.getenv("GSDESIGNNB_JENSEN_NULL_NSIM", unset = "3500"))
  lambda1_grid <- c(0.3, 0.5, 1.0, 2.0)
  k_grid <- c(0.1, 0.3, 0.5, 1.0)
  gap_days_grid <- c(10, 20, 30, 45)
}

# ---- Scenario grid ----
scenarios <- expand.grid(
  lambda1 = lambda1_grid,
  k = k_grid,
  gap_days = gap_days_grid,
  stringsAsFactors = FALSE
)
scenarios$lambda2 <- scenarios$lambda1 * rr
scenarios$gap <- scenarios$gap_days / 30.42

cat(sprintf("Total scenarios: %d\n", nrow(scenarios)))
cat(sprintf("Power simulations per scenario: %d\n", nsim_power))
cat(sprintf("Null simulations per scenario: %d\n", nsim_null))

# ---- Helper: naive sample size without Jensen correction ----
compute_naive_n <- function(design, lambda1, lambda2, k, gap,
                            test_type = c("wald", "score")) {
  test_type <- match.arg(test_type)
  if (gap == 0 || k == 0) return(design$n_total)

  inputs <- design$inputs
  ratio <- inputs$ratio
  t_bar <- design$exposure

  naive1 <- lambda1 / (1 + lambda1 * gap)
  naive2 <- lambda2 / (1 + lambda2 * gap)

  corr1 <- naive1 * (1 - k * lambda1 * gap / (1 + lambda1 * gap)^2)
  corr2 <- naive2 * (1 - k * lambda2 * gap / (1 + lambda2 * gap)^2)

  n1 <- design$n_total / (1 + ratio)
  mu1_corr <- corr1 * t_bar[1]
  mu2_corr <- corr2 * t_bar[2]

  adjusted_dispersion <- (
    design$variance * n1 - 1 / mu1_corr - 1 / (ratio * mu2_corr)
  ) / (1 + 1 / ratio)

  mu1_naive <- naive1 * t_bar[1]
  mu2_naive <- naive2 * t_bar[2]
  v1 <- (1 / mu1_naive + adjusted_dispersion) +
    (1 / ratio) * (1 / mu2_naive + adjusted_dispersion)

  theta <- log(lambda1 * inputs$rr0 / lambda2)
  z_alpha <- qnorm(1 - inputs$alpha / inputs$sided)
  z_beta <- qnorm(power_target)

  if (test_type == "score") {
    lambda0_naive <- (naive1 + ratio * naive2 * inputs$rr0) / (1 + ratio)
    mu0 <- lambda0_naive * t_bar[1]
    v0 <- (1 / mu0 + adjusted_dispersion) * (1 + 1 / ratio)
    n1_naive <- (z_alpha * sqrt(v0) + z_beta * sqrt(v1))^2 / theta^2
  } else {
    n1_naive <- (z_alpha + z_beta)^2 * v1 / theta^2
  }

  ceiling(n1_naive) + ceiling(n1_naive * ratio)
}

analyze_cut <- function(cut_dt) {
  wald <- tryCatch(
    mutze_test(cut_dt, method = "nb", sided = 1, test_type = "wald"),
    error = function(e) NULL
  )
  score <- tryCatch(
    mutze_test(cut_dt, method = "nb", sided = 1, test_type = "score"),
    error = function(e) NULL
  )

  list(
    p_wald = if (!is.null(wald)) wald$p_value else NA_real_,
    p_score = if (!is.null(score)) score$p_value else NA_real_,
    fallback_wald = if (!is.null(wald)) wald$fallback else NA_character_,
    fallback_score = if (!is.null(score)) score$fallback else NA_character_
  )
}

run_trial_pair <- function(enroll_df, fail_df, drop_df, gap, n_naive) {
  sim_data <- nb_sim(enroll_df, fail_df, drop_df,
                     max_followup = max_fu, event_gap = gap)

  cut_full <- cut_data_by_date(sim_data, cut_date = trial_dur, event_gap = gap)
  res_full <- analyze_cut(cut_full)

  if (n_naive < length(unique(sim_data$id))) {
    first_ids <- sort(unique(sim_data$id))[seq_len(n_naive)]
    sim_naive <- sim_data[sim_data$id %in% first_ids, ]
    class(sim_naive) <- class(sim_data)
    cut_naive <- cut_data_by_date(sim_naive, cut_date = trial_dur, event_gap = gap)
    res_naive <- analyze_cut(cut_naive)
  } else {
    res_naive <- res_full
  }

  data.frame(
    p_wald_corrected = res_full$p_wald,
    p_score_corrected = res_full$p_score,
    p_wald_naive = res_naive$p_wald,
    p_score_naive = res_naive$p_score,
    fallback_wald_corrected = res_full$fallback_wald,
    fallback_score_corrected = res_full$fallback_score,
    fallback_wald_naive = res_naive$fallback_wald,
    fallback_score_naive = res_naive$fallback_score,
    stringsAsFactors = FALSE
  )
}

summarize_rejections <- function(raw, prefix, alpha) {
  p_wald_corrected <- raw$p_wald_corrected
  p_score_corrected <- raw$p_score_corrected
  p_wald_naive <- raw$p_wald_naive
  p_score_naive <- raw$p_score_naive

  out <- data.frame(
    metric = prefix,
    n_valid_wald_corrected = sum(!is.na(p_wald_corrected)),
    n_valid_score_corrected = sum(!is.na(p_score_corrected)),
    n_valid_wald_naive = sum(!is.na(p_wald_naive)),
    n_valid_score_naive = sum(!is.na(p_score_naive)),
    rate_wald_corrected = mean(p_wald_corrected < alpha, na.rm = TRUE),
    rate_score_corrected = mean(p_score_corrected < alpha, na.rm = TRUE),
    rate_wald_naive = mean(p_wald_naive < alpha, na.rm = TRUE),
    rate_score_naive = mean(p_score_naive < alpha, na.rm = TRUE)
  )

  for (nm in c(
    "wald_corrected", "score_corrected",
    "wald_naive", "score_naive"
  )) {
    rate_nm <- paste0("rate_", nm)
    n_nm <- paste0("n_valid_", nm)
    se_nm <- paste0("se_", nm)
    out[[se_nm]] <- sqrt(out[[rate_nm]] * (1 - out[[rate_nm]]) / out[[n_nm]])
  }

  out$pct_fallback_mom_wald_corrected <-
    100 * mean(raw$fallback_wald_corrected == "mom", na.rm = TRUE)
  out$pct_fallback_mom_score_corrected <-
    100 * mean(raw$fallback_score_corrected == "mom", na.rm = TRUE)
  out
}

run_scenario <- function(sc_idx) {
  lambda1 <- scenarios$lambda1[sc_idx]
  lambda2 <- scenarios$lambda2[sc_idx]
  k <- scenarios$k[sc_idx]
  gap <- scenarios$gap[sc_idx]
  gap_days <- scenarios$gap_days[sc_idx]
  gap_arg <- if (gap > 0) gap else NULL

  cat(sprintf(
    "\n=== Scenario %d/%d: lambda1=%.1f, k=%.1f, gap=%d days ===\n",
    sc_idx, nrow(scenarios), lambda1, k, gap_days
  ))

  design <- tryCatch(
    sample_size_nbinom(
      lambda1 = lambda1, lambda2 = lambda2, dispersion = k,
      power = power_target, alpha = alpha,
      accrual_rate = accrual_rel, accrual_duration = accrual_dur,
      trial_duration = trial_dur, dropout_rate = dropout,
      max_followup = max_fu, event_gap = gap_arg,
      test_type = "wald"
    ),
    error = function(e) {
      cat(sprintf("  SKIPPED: %s\n", e$message))
      NULL
    }
  )
  if (is.null(design)) return(NULL)

  n_corr <- design$n_total
  n_naive <- compute_naive_n(design, lambda1, lambda2, k, gap, test_type = "wald")
  n_naive <- min(n_naive, n_corr)

  cat(sprintf("  n_corrected = %d, n_naive = %d (diff = %d)\n",
              n_corr, n_naive, n_corr - n_naive))

  enroll_df <- data.frame(rate = design$accrual_rate, duration = accrual_dur)
  drop_df <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(dropout, dropout),
    duration = c(100, 100)
  )

  make_fail_df <- function(lambda2_value) {
    data.frame(
      treatment = c("Control", "Experimental"),
      rate = c(lambda1, lambda2_value),
      dispersion = c(k, k)
    )
  }

  RNGkind("L'Ecuyer-CMRG")
  set.seed(20260502 + sc_idx)

  power_raw <- mclapply(seq_len(nsim_power), function(i) {
    run_trial_pair(enroll_df, make_fail_df(lambda2), drop_df, gap, n_naive)
  }, mc.cores = n_cores, mc.set.seed = TRUE)
  power_raw <- do.call(rbind, power_raw)

  set.seed(20270502 + sc_idx)
  null_raw <- mclapply(seq_len(nsim_null), function(i) {
    run_trial_pair(enroll_df, make_fail_df(lambda1), drop_df, gap, n_naive)
  }, mc.cores = n_cores, mc.set.seed = TRUE)
  null_raw <- do.call(rbind, null_raw)

  power_summary <- summarize_rejections(power_raw, "power", alpha)
  null_summary <- summarize_rejections(null_raw, "type1", alpha)

  scenario_row <- data.frame(
    lambda1 = lambda1,
    lambda2 = lambda2,
    k = k,
    gap_days = gap_days,
    n_corrected = n_corr,
    n_naive = n_naive,
    n_diff = n_corr - n_naive,
    correction_pct = if (gap > 0 && k > 0)
      100 * k * lambda1 * gap / (1 + lambda1 * gap)^2
    else 0,
    stringsAsFactors = FALSE
  )

  out <- cbind(
    scenario_row,
    data.frame(
      power_wald_corrected = power_summary$rate_wald_corrected,
      power_score_corrected = power_summary$rate_score_corrected,
      power_wald_naive = power_summary$rate_wald_naive,
      power_score_naive = power_summary$rate_score_naive,
      type1_wald_corrected = null_summary$rate_wald_corrected,
      type1_score_corrected = null_summary$rate_score_corrected,
      type1_wald_naive = null_summary$rate_wald_naive,
      type1_score_naive = null_summary$rate_score_naive,
      se_power_wald_corrected = power_summary$se_wald_corrected,
      se_power_score_corrected = power_summary$se_score_corrected,
      se_power_wald_naive = power_summary$se_wald_naive,
      se_power_score_naive = power_summary$se_score_naive,
      se_type1_wald_corrected = null_summary$se_wald_corrected,
      se_type1_score_corrected = null_summary$se_score_corrected,
      se_type1_wald_naive = null_summary$se_wald_naive,
      se_type1_score_naive = null_summary$se_score_naive,
      n_power = nsim_power,
      n_null = nsim_null
    )
  )

  cat(sprintf(
    "  Power Wald %.4f/%.4f; Score %.4f/%.4f (corr/naive)\n",
    out$power_wald_corrected, out$power_wald_naive,
    out$power_score_corrected, out$power_score_naive
  ))
  cat(sprintf(
    "  Type I Wald %.4f/%.4f; Score %.4f/%.4f (corr/naive)\n",
    out$type1_wald_corrected, out$type1_wald_naive,
    out$type1_score_corrected, out$type1_score_naive
  ))

  out
}

t_start <- Sys.time()
all_results <- lapply(seq_len(nrow(scenarios)), run_scenario)
t_end <- Sys.time()

all_results <- Filter(Negate(is.null), all_results)
sweep_df <- do.call(rbind, all_results)
sweep_df <- sweep_df[order(-sweep_df$correction_pct), ]

cat(sprintf("\nCompleted in %.1f minutes\n",
            as.numeric(difftime(t_end, t_start, units = "mins"))))

cat("\n=== SUMMARY TABLE ===\n")
print(sweep_df[, c(
  "lambda1", "k", "gap_days", "n_corrected", "n_naive", "n_diff",
  "power_wald_corrected", "power_wald_naive",
  "power_score_corrected", "power_score_naive",
  "type1_wald_corrected", "type1_score_corrected"
)], row.names = FALSE)

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
saveRDS(sweep_df, file = "inst/extdata/jensen_score_broad_sweep_results.rds")
cat("\nResults saved to inst/extdata/jensen_score_broad_sweep_results.rds\n")
