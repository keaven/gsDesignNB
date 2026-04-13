# Broad parameter sweep: Jensen correction impact across lambda, k, and gap
#
# Parallelized using parallel::mclapply (fork-based, macOS/Linux).
# Each scenario's 3,500 sims run in parallel across available cores.
# Total: 64 scenarios x 3,500 sims = 224,000 simulations.
#
# Output saved to inst/extdata/jensen_broad_sweep_results.rds

library(parallel)
devtools::load_all(".")

n_cores <- detectCores() - 1L
cat(sprintf("Using %d cores for parallel simulation\n", n_cores))

set.seed(20260411)

# ---- Fixed parameters ----
rr <- 0.7 # lambda2 / lambda1
alpha <- 0.025
power_target <- 0.9
dropout <- 0.1 / 12
max_fu <- 12
trial_dur <- 24
accrual_rel <- c(1, 2)
accrual_dur <- c(6, 6)
nsim <- 3500L

# ---- Scenario grid ----
scenarios <- expand.grid(
  lambda1 = c(0.3, 0.5, 1.0, 2.0),
  k = c(0.1, 0.3, 0.5, 1.0),
  gap_days = c(10, 20, 30, 45),
  stringsAsFactors = FALSE
)
scenarios$lambda2 <- scenarios$lambda1 * rr
scenarios$gap <- scenarios$gap_days / 30.42

cat(sprintf("Total scenarios: %d\n", nrow(scenarios)))

# ---- Helper: naive sample size (without Jensen correction) ----
compute_naive_n <- function(design, lambda1, lambda2, k, gap) {
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

  A <- design$variance * n1 - 1 / mu1_corr - 1 / (ratio * mu2_corr)
  k_adj <- A / (1 + 1 / ratio)

  mu1_naive <- naive1 * t_bar[1]
  mu2_naive <- naive2 * t_bar[2]
  V_naive <- (1 / mu1_naive + k_adj) + (1 / ratio) * (1 / mu2_naive + k_adj)

  theta <- log(lambda2 / lambda1)
  z_a <- qnorm(1 - alpha)
  z_b <- qnorm(power_target)
  n1_naive <- (z_a + z_b)^2 * V_naive / theta^2
  ceiling(n1_naive) + ceiling(n1_naive * ratio)
}

# ---- Helper: analyze one cut ----
analyze_cut <- function(cut_dt) {
  res <- tryCatch(mutze_test(cut_dt, method = "nb", sided = 1), error = function(e) NULL)
  if (is.null(res)) return(NULL)
  list(p_value = res$p_value)
}

# ---- Run one scenario (parallelized internally) ----
run_scenario <- function(sc_idx) {
  lambda1 <- scenarios$lambda1[sc_idx]
  lambda2 <- scenarios$lambda2[sc_idx]
  k <- scenarios$k[sc_idx]
  gap <- scenarios$gap[sc_idx]
  gap_days <- scenarios$gap_days[sc_idx]
  gap_arg <- if (gap > 0) gap else NULL

  cat(sprintf("\n=== Scenario %d/%d: lambda1=%.1f, k=%.1f, gap=%d days ===\n",
              sc_idx, nrow(scenarios), lambda1, k, gap_days))

  design <- tryCatch(
    sample_size_nbinom(
      lambda1 = lambda1, lambda2 = lambda2, dispersion = k,
      power = power_target, alpha = alpha,
      accrual_rate = accrual_rel, accrual_duration = accrual_dur,
      trial_duration = trial_dur, dropout_rate = dropout,
      max_followup = max_fu, event_gap = gap_arg
    ),
    error = function(e) {
      cat(sprintf("  SKIPPED: %s\n", e$message))
      NULL
    }
  )
  if (is.null(design)) return(NULL)

  n_corr <- design$n_total
  n_naive <- compute_naive_n(design, lambda1, lambda2, k, gap)

  cat(sprintf("  n_corrected = %d, n_naive = %d (diff = %d)\n",
              n_corr, n_naive, n_corr - n_naive))

  enroll_df <- data.frame(rate = design$accrual_rate, duration = accrual_dur)
  fail_df <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(lambda1, lambda2),
    dispersion = k
  )
  drop_df <- data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(dropout, dropout),
    duration = c(100, 100)
  )

  # Use L'Ecuyer-CMRG for reproducible parallel streams
  RNGkind("L'Ecuyer-CMRG")
  set.seed(20260411 + sc_idx)
  sim_results <- mclapply(seq_len(nsim), function(i) {
    sim_data <- nb_sim(enroll_df, fail_df, drop_df,
                       max_followup = max_fu, event_gap = gap)
    # Full analysis (corrected n)
    cut_full <- cut_data_by_date(sim_data, cut_date = trial_dur, event_gap = gap)
    res_full <- analyze_cut(cut_full)

    # Naive analysis (first n_naive subjects)
    if (n_naive < n_corr) {
      first_ids <- sort(unique(sim_data$id))[seq_len(min(n_naive, length(unique(sim_data$id))))]
      sim_naive <- sim_data[sim_data$id %in% first_ids, ]
      class(sim_naive) <- class(sim_data)
      cut_naive <- cut_data_by_date(sim_naive, cut_date = trial_dur, event_gap = gap)
      res_naive <- analyze_cut(cut_naive)
    } else {
      res_naive <- res_full
    }

    if (!is.null(res_full) && !is.null(res_naive)) {
      list(valid = TRUE,
           rej_corr = res_full$p_value < alpha,
           rej_naive = res_naive$p_value < alpha)
    } else {
      list(valid = FALSE, rej_corr = FALSE, rej_naive = FALSE)
    }
   }, mc.cores = n_cores, mc.set.seed = TRUE)

  valid <- vapply(sim_results, function(x) x$valid, logical(1))
  rej_corr <- vapply(sim_results, function(x) x$rej_corr, logical(1))
  rej_naive <- vapply(sim_results, function(x) x$rej_naive, logical(1))

  n_valid <- sum(valid)
  power_corr <- mean(rej_corr[valid])
  power_naive <- mean(rej_naive[valid])

  corr_only <- sum(rej_corr[valid] & !rej_naive[valid])
  naive_only <- sum(!rej_corr[valid] & rej_naive[valid])

  ci_corr <- binom.test(sum(rej_corr[valid]), n_valid)$conf.int
  ci_naive <- binom.test(sum(rej_naive[valid]), n_valid)$conf.int

  d <- as.integer(rej_corr[valid]) - as.integer(rej_naive[valid])
  diff_mean <- mean(d)
  diff_se <- sd(d) / sqrt(length(d))

  correction_pct <- 100 * k * lambda1 * gap / (1 + lambda1 * gap)^2

  cat(sprintf("  Power corrected: %.4f [%.4f, %.4f]\n", power_corr, ci_corr[1], ci_corr[2]))
  cat(sprintf("  Power naive:     %.4f [%.4f, %.4f]\n", power_naive, ci_naive[1], ci_naive[2]))
  cat(sprintf("  Paired diff:     %.4f pp (SE=%.4f)\n", 100 * diff_mean, 100 * diff_se))

  data.frame(
    lambda1 = lambda1,
    lambda2 = lambda2,
    k = k,
    gap_days = gap_days,
    n_corrected = n_corr,
    n_naive = n_naive,
    n_diff = n_corr - n_naive,
    n_valid = n_valid,
    power_corrected = power_corr,
    power_naive = power_naive,
    power_diff = diff_mean,
    power_diff_se = diff_se,
    ci_corr_lo = ci_corr[1],
    ci_corr_hi = ci_corr[2],
    ci_naive_lo = ci_naive[1],
    ci_naive_hi = ci_naive[2],
    corr_only = corr_only,
    naive_only = naive_only,
    correction_pct = correction_pct,
    stringsAsFactors = FALSE
  )
}

# ---- Run all scenarios sequentially (each parallelized internally) ----
t_start <- Sys.time()
all_results <- lapply(seq_len(nrow(scenarios)), run_scenario)
t_end <- Sys.time()

all_results <- Filter(Negate(is.null), all_results)
sweep_df <- do.call(rbind, all_results)

cat(sprintf("\n\nCompleted in %.1f minutes\n", as.numeric(difftime(t_end, t_start, units = "mins"))))

# ---- Summary ----
cat("\n=== SUMMARY TABLE ===\n")
sweep_df <- sweep_df[order(-sweep_df$correction_pct), ]
print(sweep_df[, c("lambda1", "k", "gap_days", "n_corrected", "n_naive", "n_diff",
                    "power_corrected", "power_naive", "power_diff", "correction_pct")],
      row.names = FALSE)

# ---- Save ----
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
saveRDS(sweep_df, file = "inst/extdata/jensen_broad_sweep_results.rds")
cat("\nResults saved to inst/extdata/jensen_broad_sweep_results.rds\n")
