# Scenario sweep: Jensen correction impact across k and gap values
#
# For each scenario, simulates at the corrected n and analyzes at both
# corrected and naive sample sizes (first n_naive subjects).
# Results are saved locally (not pushed to git due to size).

library(data.table)
data.table::setDTthreads(1)
devtools::load_all(".")

set.seed(54321)

# Fixed parameters
lambda1 <- 0.4; lambda2 <- 0.3
alpha <- 0.025; power_target <- 0.9
dropout <- 0.1 / 12; max_fu <- 12; trial_dur <- 24
accrual_rel <- c(1, 2); accrual_dur <- c(6, 6)

# Scenario grid: k x gap
# Only scenarios where both k > 0 AND gap > 0 are interesting;
# when either is zero the Jensen correction has no effect and
# naive = corrected by construction.
scenarios <- expand.grid(
  k = c(0.5, 1.0),
  gap_days = c(20, 30),
  stringsAsFactors = FALSE
)
scenarios$gap <- scenarios$gap_days / 30.42

nsim <- 10000 # per scenario

# Helper: compute naive sample size (without Jensen correction)
# by using the naive effective rate but same exposure/Q as the corrected design
compute_naive_n <- function(design, k, gap) {
  if (gap == 0 || k == 0) return(design$n_total)

  inputs <- design$inputs
  ratio <- inputs$ratio
  t_bar <- design$exposure

  # Naive effective rates
  naive1 <- lambda1 / (1 + lambda1 * gap)
  naive2 <- lambda2 / (1 + lambda2 * gap)

  # Corrected effective rates (what sample_size_nbinom used)
  corr1 <- naive1 * (1 - k * lambda1 * gap / (1 + lambda1 * gap)^2)
  corr2 <- naive2 * (1 - k * lambda2 * gap / (1 + lambda2 * gap)^2)

  # The Q-adjusted dispersion is the same (it doesn't depend on gap correction)
  # Back out from design: mu_corr = lambda_corr_eff * t_bar
  # V_corr = (1/mu1_corr + k_adj)/n1 + (1/mu2_corr + k_adj)/n2
  # We need k_adj = k * Q. Extract from: design$variance * n1 = (1/mu1 + k_adj) + ...

  # Simpler: recompute using the design's exposure and same k*Q
  # Since sample_size_nbinom reports exposure (t_bar) and we know k, Q = k_adj/k
  # But k_adj is internal. Let me extract it from the relationship:
  # variance = (1/mu1 + k_adj)/n1 + (1/mu2 + k_adj)/n2
  # with mu1 = corr1 * t_bar[1], mu2 = corr2 * t_bar[2]
  n1 <- design$n_total / (1 + ratio)
  n2 <- design$n_total * ratio / (1 + ratio)
  mu1_corr <- corr1 * t_bar[1]
  mu2_corr <- corr2 * t_bar[2]
  # k_adj from: V * n1 = (1/mu1 + k_adj) + (1/ratio)(1/mu2 + k_adj)
  # design$variance * n1 = (1/mu1_corr + k_adj) + (1/ratio)(1/mu2_corr + k_adj)
  # Solve for k_adj:
  # Let A = design$variance * n1 - 1/mu1_corr - 1/(ratio * mu2_corr)
  # A = k_adj + k_adj/ratio = k_adj * (1 + 1/ratio)
  A <- design$variance * n1 - 1 / mu1_corr - 1 / (ratio * mu2_corr)
  k_adj <- A / (1 + 1 / ratio)

  # Now compute naive variance factor
  mu1_naive <- naive1 * t_bar[1]
  mu2_naive <- naive2 * t_bar[2]
  V_naive <- (1 / mu1_naive + k_adj) + (1 / ratio) * (1 / mu2_naive + k_adj)

  theta <- log(lambda1 * 1 / lambda2) # rr0 = 1
  z_a <- qnorm(1 - alpha)
  z_b <- qnorm(power_target)
  n1_naive <- (z_a + z_b)^2 * V_naive / theta^2
  n_total_naive <- ceiling(n1_naive) + ceiling(n1_naive * ratio)
  n_total_naive
}

# Helper: analyze one cut dataset
analyze_cut <- function(cut_dt) {
  res <- tryCatch(mutze_test(cut_dt, method = "nb", sided = 1), error = function(e) NULL)
  if (is.null(res)) return(NULL)
  list(
    p_value = res$p_value,
    estimate = res$estimate,
    se = res$se
  )
}

# Run sweep
all_results <- list()

for (sc in seq_len(nrow(scenarios))) {
  k <- scenarios$k[sc]
  gap <- scenarios$gap[sc]
  gap_arg <- if (gap > 0) gap else NULL

  cat(sprintf("\n=== Scenario %d/%d: k=%.1f, gap=%d days ===\n",
              sc, nrow(scenarios), k, scenarios$gap_days[sc]))

  # Corrected design
  design <- sample_size_nbinom(
    lambda1 = lambda1, lambda2 = lambda2, dispersion = k,
    power = power_target, alpha = alpha,
    accrual_rate = accrual_rel, accrual_duration = accrual_dur,
    trial_duration = trial_dur, dropout_rate = dropout,
    max_followup = max_fu, event_gap = gap_arg
  )
  n_corr <- design$n_total
  n_naive <- compute_naive_n(design, k, gap)

  cat(sprintf("  n_corrected = %d, n_naive = %d (diff = %d)\n", n_corr, n_naive, n_corr - n_naive))

  # Setup simulation parameters
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

  rej_corr <- logical(nsim)
  rej_naive <- logical(nsim)
  valid <- logical(nsim)

  for (i in seq_len(nsim)) {
    if (i %% 2000 == 0) cat(sprintf("  %d / %d\n", i, nsim))

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
      valid[i] <- TRUE
      rej_corr[i] <- res_full$p_value < alpha
      rej_naive[i] <- res_naive$p_value < alpha
    }
  }

  n_valid <- sum(valid)
  power_corr <- mean(rej_corr[valid])
  power_naive <- mean(rej_naive[valid])

  # Paired comparison
  both <- sum(rej_corr[valid] & rej_naive[valid])
  corr_only <- sum(rej_corr[valid] & !rej_naive[valid])
  naive_only <- sum(!rej_corr[valid] & rej_naive[valid])

  ci_corr <- binom.test(sum(rej_corr[valid]), n_valid)$conf.int
  ci_naive <- binom.test(sum(rej_naive[valid]), n_valid)$conf.int

  # Paired CI for difference
  d <- as.integer(rej_corr[valid]) - as.integer(rej_naive[valid])
  diff_mean <- mean(d)
  diff_se <- sd(d) / sqrt(length(d))

  cat(sprintf("  Power corrected: %.4f [%.4f, %.4f]\n", power_corr, ci_corr[1], ci_corr[2]))
  cat(sprintf("  Power naive:     %.4f [%.4f, %.4f]\n", power_naive, ci_naive[1], ci_naive[2]))
  cat(sprintf("  Paired diff:     %.4f (SE=%.4f)\n", diff_mean, diff_se))
  cat(sprintf("  Discordant: corr_only=%d, naive_only=%d\n", corr_only, naive_only))

  all_results[[sc]] <- data.frame(
    k = k,
    gap_days = scenarios$gap_days[sc],
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
    correction_pct = if (gap > 0 && k > 0)
      round(100 * k * lambda1 * gap / (1 + lambda1 * gap)^2, 1)
    else 0,
    stringsAsFactors = FALSE
  )
}

sweep_df <- do.call(rbind, all_results)

cat("\n\n=== SUMMARY TABLE ===\n")
print(sweep_df[, c("k", "gap_days", "n_corrected", "n_naive", "n_diff",
                    "power_corrected", "power_naive", "power_diff", "correction_pct")])

# Save locally (not for git push)
saveRDS(sweep_df, file = "inst/extdata/scenario_sweep_results.rds")
cat("\nResults saved to inst/extdata/scenario_sweep_results.rds\n")
