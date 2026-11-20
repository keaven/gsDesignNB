summarize_spending_thresholds <- function(result, thresholds = c(0.95, 1)) {
  if (result$cell$cut_rule == "ssr" || is.null(result$raw)) return(NULL)
  raw <- as.data.frame(result$raw)
  metrics <- c("info_blinded_ml", "info_unblinded_ml")
  n_paths <- length(unique(raw$sim))
  stopifnot(n_paths == result$cell$n_sims)
  final_look <- result$design$k
  active <- !is.na(raw$status) & raw$status == "analyzed"
  projected <- raw$info_blinded_final_projected
  do.call(rbind, lapply(metrics, function(metric) {
    current <- raw[[metric]]
    valid <- active & is.finite(current) & current > 0 &
      is.finite(projected) & projected > 0
    time <- rep(NA_real_, nrow(raw))
    time[valid] <- current[valid] / projected[valid]
    do.call(rbind, lapply(thresholds, function(threshold) {
      reached <- valid & !is.na(time) & time >= threshold
      early <- unique(raw$sim[reached & raw$analysis < final_look])
      final <- unique(raw$sim[reached & raw$analysis == final_look])
      rate <- length(early) / n_paths
      data.frame(
        result$cell, spending = "expected_final", information_estimator = metric,
        final_information_estimator = "info_blinded_final_projected",
        threshold = threshold, n_paths = n_paths,
        n_early = length(early), early_rate = rate,
        early_mcse = sqrt(rate * (1 - rate) / n_paths),
        n_final_only = length(setdiff(final, early)),
        n_without_observed_hit = n_paths - length(union(early, final)),
        n_invalid_interim = length(unique(raw$sim[active & !valid & raw$analysis < final_look])),
        n_invalid_final = length(unique(raw$sim[active & !valid & raw$analysis == final_look])),
        n_skipped_interim = length(unique(raw$sim[!active & raw$analysis < final_look]))
      )
    }))
  }))
}

checkpoint_dir <- file.path("paper", "gs-strategy-results", "cells")
paths <- list.files(checkpoint_dir, pattern = "\\.rds$", full.names = TRUE)
source_signature <- tools::md5sum(c(
  "R/gsNBPower.R", "R/check_gs_bound.R", "R/sim_gs_nbinom.R",
  "R/blinded_ssr.R", "R/get_cut_date.R", "data-raw/gs_strategy_study.R"
))
checkpoints <- lapply(paths, function(path) {
  result <- readRDS(path)
  expected_n <- if (result$cell$truth == "null") 20000L else 5000L
  if (!identical(result$source_signature, source_signature) ||
      result$cell$n_sims != expected_n) {
    return(NULL)
  }
  result$spending_threshold_summary <- summarize_spending_thresholds(result)
  result[c("cell", "summary", "look_summary", "ssr_summary", "spending_threshold_summary")]
})
included <- !vapply(checkpoints, is.null, logical(1))
checkpoints <- checkpoints[included]
if (!length(checkpoints)) stop("No production checkpoints match the current source.")
bind_component <- function(component) {
  do.call(rbind, lapply(checkpoints, `[[`, component))
}
cells <- bind_component("cell")
stopifnot(!anyDuplicated(cells[c("tier", "scenario", "timing_method", "truth",
                                "nuisance", "test_type", "cut_rule")]))
expected_by_tier <- c(A = 32L, B = 48L, C = 12L)
coverage <- data.frame(
  tier = names(expected_by_tier),
  completed = as.integer(table(factor(cells$tier, levels = names(expected_by_tier)))),
  planned = unname(expected_by_tier)
)
snapshot <- list(
  summary = bind_component("summary"),
  look_summary = bind_component("look_summary"),
  ssr_summary = bind_component("ssr_summary"),
  spending_threshold_summary = bind_component("spending_threshold_summary"),
  cells = cells,
  coverage = coverage,
  expected_cells = sum(expected_by_tier),
  complete = all(coverage$completed == coverage$planned),
  excluded_checkpoints = sum(!included),
  mode = "production",
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  source_signature = source_signature,
  session_info = utils::sessionInfo()
)
if ("--spending-thresholds-only" %in% commandArgs(trailingOnly = TRUE)) {
  saveRDS(snapshot, file.path("paper", "gs-strategy-results", "spending_threshold_summary.rds"))
} else {
  saveRDS(snapshot, file.path("inst", "extdata", "gs_strategy_summary.rds"))
}
thresholds <- snapshot$spending_threshold_summary
report_columns <- c("tier", "scenario", "timing_method", "truth", "nuisance", "test_type",
                    "cut_rule", "information_estimator", "threshold", "n_paths", "n_early",
                    "early_rate", "early_mcse", "n_final_only", "n_without_observed_hit",
                    "n_invalid_interim", "n_invalid_final", "n_skipped_interim")
report <- thresholds[report_columns]
report$early_rate <- sprintf("%.2f%%", 100 * report$early_rate)
report$early_mcse <- sprintf("%.3f pp", 100 * report$early_mcse)
report_lines <- apply(report, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
writeLines(c(
  "# Near-final spending-time diagnostics", "",
  paste("Generated:", snapshot$generated_at), "",
  "These are potential skipped-interim counts, not operating characteristics of the revised skip-interim/retain-final rule. Earlier efficacy and futility stops are not applied, so a candidate may occur after a trial would already have stopped. Use data-raw/evaluate_interim_skip.R for sequential reanalysis with actual skip counts. No new trials were simulated.", "",
  "Spending time is the uncapped current information divided by blinded projected-final information (the expected_final strategy). Both blinded and unblinded current-information estimates are shown separately. Thresholds 0.95 and 1.00 are inclusive. These counts do not describe the calendar or planned_final spending strategies.", "",
  "n_early counts distinct trials reaching the threshold before the scheduled final look; n_final_only counts trials reaching it at the final look but at no interim. n_without_observed_hit includes trials with invalid or skipped looks and must not be interpreted as confirmed non-crossing. Invalid and skipped counts can overlap threshold hits at other looks. All rates use n_paths as the denominator; early_mcse is the Monte Carlo standard error of early_rate.", "",
  paste0("Compatible production cells: ", nrow(cells), "/", snapshot$expected_cells,
         ". Cells with a non-SSR threshold diagnostic: ",
         length(unique(do.call(paste, thresholds[c("tier", "scenario", "timing_method", "truth",
                                                  "nuisance", "test_type", "cut_rule")]))),
         ". Adaptive SSR cells are not screened by this diagnostic."), "",
  paste0("| ", paste(report_columns, collapse = " | "), " |"),
  paste0("| ", paste(rep("---", length(report_columns)), collapse = " | "), " |"),
  report_lines
), file.path("paper", "gs-strategy-results", "spending-threshold-summary.md"))
print(coverage, row.names = FALSE)
cat("Excluded incompatible or non-production checkpoints:", sum(!included), "\n")
cat("Full study complete:", snapshot$complete, "\n")
cat("Spending-threshold diagnostics: paper/gs-strategy-results/spending-threshold-summary.md\n")