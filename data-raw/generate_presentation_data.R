devtools::load_all(".", quiet = TRUE)
library(ggplot2)
library(rlang)
output_dir <- file.path("paper", "gsDesignNB-talk_files")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
colors <- c("Score" = "#007A73", "Wald" = "#C34633", "Naive" = "#C34633",
            "Corrected" = "#007A73", "No adaptation" = "#5B6269",
            "Blinded SSR" = "#007A73", "Unblinded SSR" = "#C34633")
theme_set(theme_minimal(base_size = 16) + theme(
  panel.grid.minor = element_blank(), legend.position = "bottom",
  plot.background = element_rect(fill = "white", color = NA),
  axis.title = element_text(color = "#232729"), plot.margin = margin(12, 16, 12, 12)
))
save_plot <- function(plot, name, width = 9, height = 4.7) {
  ggsave(file.path(output_dir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 180, bg = "white")
}
percentage <- scales::label_percent(accuracy = 1)
cache <- readRDS("inst/extdata/score_sweep_summary.rds")
jensen <- readRDS("inst/extdata/jensen_broad_sweep_results.rds")
jensen_score <- readRDS("inst/extdata/jensen_score_broad_sweep_results.rds")
strategy <- readRDS("inst/extdata/gs_strategy_summary.rds")
ssr <- readRDS("inst/extdata/ssr_sim_vignette_summary_score.rds")
fixed <- sample_size_nbinom(
  lambda1 = 0.5, lambda2 = 0.35, dispersion = 0.5, power = 0.90, alpha = 0.025,
  accrual_rate = c(12, 24), accrual_duration = c(6, 6), trial_duration = 24,
  dropout_rate = 0.10 / 12, max_followup = 12, event_gap = 20 / 30.4375,
  test_type = "score"
)
design <- fixed |>
  gsNBCalendar(k = 3, test.type = 4, alpha = 0.025, sfu = sfHSD, sfupar = -2,
               sfl = sfCauchy, sflpar = c(0.4, 0.75, 0.56, 0.63),
               analysis_times = c(9, 14, 24)) |>
  gsDesignNB::toInteger()
information <- expand.grid(mu = seq(0, 15, length.out = 150), k = c(0, 0.2, 0.5, 1))
information$value <- information$mu / (1 + information$k * information$mu)
information$dispersion <- factor(information$k, labels = c("Poisson: k = 0", "k = 0.2", "k = 0.5", "k = 1"))
save_plot(ggplot(information, aes(mu, value, color = dispersion)) + geom_line(linewidth = 1.3) +
            scale_color_manual(values = c("#C34633", "#007A73", "#A38224", "#5B6269")) +
            labs(x = "Expected events per subject", y = "Information per subject", color = NULL), "information")
gap <- expand.grid(lambda = seq(0.01, 3, length.out = 150), k = c(0.1, 0.5, 1))
gap$correction <- gap$k * gap$lambda / (1 + gap$lambda)^2
save_plot(ggplot(gap, aes(lambda, correction, color = factor(k))) + geom_line(linewidth = 1.3) +
            scale_color_manual(values = c("#5B6269", "#007A73", "#C34633")) +
            scale_y_continuous(labels = percentage) +
            labs(x = "Control rate per month (one-month gap)", y = "Effective-rate correction", color = "Dispersion"), "gap-correction")
save_plot(ggplot(jensen, aes(correction_pct, 100 * power_diff, color = factor(k))) +
            geom_hline(yintercept = 0, color = "#A9ADAE") + geom_point(size = 3, alpha = 0.8) +
            scale_color_manual(values = c("#5B6269", "#A38224", "#007A73", "#C34633")) +
            labs(x = "Theoretical effective-rate correction (%)", y = "Paired power gain (percentage points)", color = "Dispersion"), "jensen-wald")
save_plot(ggplot(jensen_score, aes(power_wald_corrected, power_score_corrected, color = factor(gap_days))) +
            geom_abline(slope = 1, intercept = 0, color = "#A9ADAE") +
            geom_hline(yintercept = 0.9, linetype = "dashed", color = "#007A73") +
            geom_point(size = 3, alpha = 0.8) +
            scale_color_manual(values = c("#5B6269", "#A38224", "#007A73", "#C34633")) +
            scale_x_continuous(labels = percentage) + scale_y_continuous(labels = percentage) +
            labs(x = "Wald-test power", y = "Score-test power", color = "Gap (days)"), "jensen-score")
long_rates <- function(data) {
  rbind(data.frame(sizing = data$sizing, test = "Wald", rate = data$rate_wald),
        data.frame(sizing = data$sizing, test = "Score", rate = data$rate_score))
}
null_rates <- long_rates(cache$null_summary)
power_rates <- long_rates(cache$power_summary)
rate_plot <- function(data, reference, ylabel) {
  ggplot(data, aes(.data$sizing, .data$rate, color = .data$test)) +
    geom_hline(yintercept = reference, linetype = "dashed", color = "#5B6269") +
    geom_boxplot(position = position_dodge(width = 0.6), width = 0.5, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.6, seed = 1),
               alpha = 0.6, size = 2) +
    scale_color_manual(values = colors) +
    scale_y_continuous(labels = scales::label_percent(accuracy = if (reference == 0.025) 0.01 else 0.1)) +
    scale_x_discrete(labels = c(score = "Score-sized", wald = "Wald-sized")) +
    labs(x = NULL, y = ylabel, color = "Final test")
}
save_plot(rate_plot(null_rates, 0.025, "Null rejection probability"), "score-null")
save_plot(rate_plot(power_rates, 0.90, "Power at RR = 0.70"), "score-power")
looks <- data.frame(month = design$T, information = design$timing,
                    enrolled = design$n_total / utils::tail(design$n_total, 1))
look_plot <- rbind(data.frame(month = looks$month, fraction = looks$information, metric = "Information"),
                   data.frame(month = looks$month, fraction = looks$enrolled, metric = "Enrollment"))
save_plot(ggplot(look_plot, aes(month, fraction, color = metric)) +
            geom_line(linewidth = 1.3) + geom_point(size = 4) +
            scale_x_continuous(breaks = design$T) + scale_y_continuous(labels = percentage, limits = c(0, 1)) +
            scale_color_manual(values = c(Information = "#007A73", Enrollment = "#C34633")) +
            labs(x = "Planned analysis month", y = "Fraction of final target", color = NULL), "calendar")
sensitivity <- expand.grid(rr = seq(0.45, 1.1, by = 0.05), k = c(0.5, 1))
sensitivity$power <- vapply(seq_len(nrow(sensitivity)), function(index) {
  gsNBPower(design, lambda2 = 0.5 * sensitivity$rr[index], dispersion = sensitivity$k[index],
            test_type = "score", boundary_update = "fixed")$power
}, numeric(1))
save_plot(ggplot(sensitivity, aes(rr, power, color = factor(k))) + geom_line(linewidth = 1.3) +
            geom_hline(yintercept = 0.9, linetype = "dashed", color = "#A9ADAE") +
            geom_vline(xintercept = 0.7, linetype = "dotted", color = "#A9ADAE") +
            scale_color_manual(values = c("#007A73", "#C34633")) +
            scale_y_continuous(labels = percentage, breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
            labs(x = "True experimental/control rate ratio", y = "Asymptotic power; no resizing", color = "Dispersion"), "sensitivity")
coverage <- strategy$coverage
coverage$label <- sprintf("%d / %d", coverage$completed, coverage$planned)
coverage_plot <- rbind(data.frame(tier = coverage$tier, value = coverage$completed, status = "Completed"),
                       data.frame(tier = coverage$tier, value = coverage$planned - coverage$completed,
                                  status = "Remaining planned"))
coverage_plot$status <- factor(coverage_plot$status, levels = c("Remaining planned", "Completed"))
save_plot(ggplot(coverage_plot, aes(tier, value, fill = status)) + geom_col(width = 0.6) +
            geom_text(data = coverage, aes(x = tier, y = planned + 3, label = label),
                      inherit.aes = FALSE, size = 6) +
            scale_fill_manual(values = c("Completed" = "#007A73", "Remaining planned" = "#E3E7E7")) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
            scale_x_discrete(labels = c(A = "A: Planned nuisance", B = "B: Misspecification", C = "C: Timing / SSR")) +
            labs(x = NULL, y = "Simulation cells", fill = NULL), "coverage")
null_strategy <- subset(strategy$summary, tier == "A" & truth == "null" &
                          policy != "legacy" & information_estimator == "info_unblinded_ml")
null_strategy$failure_rate <- null_strategy$n_failed / null_strategy$n_sims
worst <- null_strategy[which.max(null_strategy$failure_rate), ]
failure <- data.frame(metric = c("Observed rejection", "Failed paths", "Rejection or failure"),
                      probability = c(worst$rejection_rate, worst$failure_rate, worst$rejection_or_failure))
failure$metric <- factor(failure$metric, levels = failure$metric)
save_plot(ggplot(failure, aes(metric, probability, fill = metric)) + geom_col(width = 0.6) +
            geom_hline(yintercept = 0.025, linetype = "dashed", color = "#5B6269") +
            geom_text(aes(label = scales::percent(probability, accuracy = 0.1)), vjust = -0.6, size = 6) +
            scale_fill_manual(values = c("#007A73", "#C34633", "#A38224")) +
            scale_y_continuous(labels = percentage, expand = expansion(mult = c(0, 0.15))) +
            labs(x = NULL, y = "Fraction of 20,000 simulated paths") + theme(legend.position = "none"), "failure")
ssr_null <- do.call(rbind, ssr$null_nonbinding_by_test)
ssr_null$Test <- tools::toTitleCase(ssr_null$analysis_test)
ssr_null$lower <- vapply(seq_len(nrow(ssr_null)), function(index) {
  binom.test(round(ssr_null$type1_error[index] * ssr_null$n_sims[index]), ssr_null$n_sims[index])$conf.int[1]
}, numeric(1))
ssr_null$upper <- vapply(seq_len(nrow(ssr_null)), function(index) {
  binom.test(round(ssr_null$type1_error[index] * ssr_null$n_sims[index]), ssr_null$n_sims[index])$conf.int[2]
}, numeric(1))
save_plot(ggplot(ssr_null, aes(strategy, type1_error, color = Test)) +
            geom_hline(yintercept = 0.025, linetype = "dashed", color = "#5B6269") +
            geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12, position = position_dodge(width = 0.5)) +
            geom_point(size = 3, position = position_dodge(width = 0.5)) +
            facet_wrap(~ k_true, labeller = label_both) + scale_color_manual(values = colors) +
            scale_y_continuous(labels = scales::label_percent(accuracy = 0.1)) +
            scale_x_discrete(labels = c("No adaptation" = "None", "Blinded SSR" = "Blinded", "Unblinded SSR" = "Unblinded")) +
            labs(x = "Adaptation strategy", y = "Archived null rejection", color = "Test"), "ssr-null")
ssr_power <- subset(ssr$summary_dt, abs(rr_true - 0.7) < 0.0001)
save_plot(ggplot(ssr_power, aes(strategy, rejection_rate, color = strategy)) +
            geom_hline(yintercept = 0.9, linetype = "dashed", color = "#5B6269") +
            geom_boxplot(width = 0.5, outlier.shape = NA) +
            geom_point(position = position_jitter(width = 0.12, seed = 1), size = 2, alpha = 0.7) +
            scale_color_manual(values = colors) + scale_y_continuous(labels = percentage) +
            labs(x = NULL, y = "Archived score-test power at RR = 0.70") + theme(legend.position = "none"), "ssr-power")
payload <- list(
  generated_at = strategy$generated_at, coverage = coverage,
  completed_cells = nrow(strategy$cells), planned_cells = strategy$expected_cells,
  excluded = strategy$excluded_checkpoints,
  fixed_n = fixed$n_total, gs_n = utils::tail(design$n_total, 1),
  final_information = utils::tail(design$n.I, 1),
  final_events = utils::tail(design$events, 1),
  score_null_range = range(null_rates$rate[null_rates$test == "Score"]),
  wald_null_range = range(null_rates$rate[null_rates$test == "Wald"]),
  jensen_gain = max(jensen$power_diff),
  jensen_score_range = range(jensen_score$power_score_corrected),
  failure_overlap = worst$rejection_rate + worst$failure_rate - worst$rejection_or_failure,
  worst_failure = worst[c("scenario", "timing_method", "spending", "n_sims",
                          "rejection_rate", "failure_rate", "rejection_or_failure")]
)
jsonlite::write_json(payload, file.path(output_dir, "data.json"), auto_unbox = TRUE,
                     pretty = TRUE, digits = 8)
stopifnot(length(list.files(output_dir, pattern = "\\.png$")) == 12L,
          nrow(ssr_power) > 0, nrow(strategy$cells) == sum(coverage$completed))
cat("Generated 12 presentation figures and data.json without new simulations.\n")