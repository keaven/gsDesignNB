evaluate_interim_skip <- function(raw, design, info_col, threshold = 0.95,
                                  ignore_futility = FALSE,
                                  final_info_col = "info_blinded_final_projected") {
  if (length(threshold) != 1L || !is.finite(threshold) ||
      threshold <= 0 || threshold > 1) {
    stop("threshold must be a single number in (0, 1].")
  }
  raw <- as.data.frame(raw)
  required <- c("sim", "analysis", "status", info_col, final_info_col)
  if (!all(required %in% names(raw))) stop("Missing monitoring columns.")
  current <- raw[[info_col]]
  projected <- raw[[final_info_col]]
  active <- !is.na(raw$status) & raw$status == "analyzed"
  valid <- is.finite(current) & current > 0 & is.finite(projected) & projected > 0
  raw$uncapped_spending_time <- NA_real_
  raw$uncapped_spending_time[valid] <- current[valid] / projected[valid]
  raw$near_final_candidate <- active & raw$analysis < design$k & valid &
    !is.na(raw$uncapped_spending_time) & raw$uncapped_spending_time >= threshold
  raw$status[raw$near_final_candidate] <- "skipped"
  final <- active & raw$analysis == design$k
  raw[[final_info_col]][final] <- current[final]
  checked <- gsDesignNB::check_gs_bound(
    raw, design, info_col = info_col, boundary_update = "recompute",
    spending = "expected_final", final_info_col = final_info_col,
    ignore_futility = ignore_futility
  )
  groups <- split(seq_len(nrow(checked)), checked$sim)
  for (rows in groups) {
    stopped <- rows[checked$cross_upper[rows] | checked$cross_lower[rows] |
                      checked$cross_harm[rows]]
    if (length(stopped)) {
      after_stop <- rows[checked$analysis[rows] > min(checked$analysis[stopped])]
      checked$bound_status[after_stop] <- "not_evaluated"
    }
  }
  checked$near_final_skip <- checked$near_final_candidate & checked$bound_status == "skipped"
  checked
}

check_interim_skip <- function() {
  design <- gsDesign::gsDesign(
    k = 3, test.type = 1, alpha = 0.025, beta = 0.1, n.fix = 100,
    timing = c(0.4, 0.75, 1), sfu = gsDesign::sfHSD, sfupar = -2
  )
  raw <- data.frame(
    sim = 1L, analysis = 1:3, analysis_time = c(10, 18, 24),
    status = "analyzed", z_stat = 0, information = c(40, 95, 100),
    info_blinded_final_projected = 100
  )
  checked <- evaluate_interim_skip(raw, design, "information")
  stopifnot(identical(checked$bound_status, c("ok", "skipped", "ok")),
            identical(checked$near_final_skip, c(FALSE, TRUE, FALSE)),
            is.na(checked$spending_time[2]), checked$spending_time[3] == 1,
            checked$analysis_time[3] == 24)
  baseline <- gsDesignNB::check_gs_bound(
    raw, design, info_col = "information", spending = "expected_final"
  )
  stopifnot(isTRUE(all.equal(checked$upper_bound[1], baseline$upper_bound[1])))
  below <- raw
  below$information[2] <- 94.99
  stopifnot(!any(evaluate_interim_skip(below, design, "information")$near_final_skip))
  above <- raw
  above$information[2] <- 101
  stopifnot(evaluate_interim_skip(above, design, "information")$near_final_skip[2])
  stopped <- raw
  stopped$z_stat[1] <- 10
  stopped <- evaluate_interim_skip(stopped, design, "information")
  stopifnot(stopped$cross_upper[1], !any(stopped$near_final_skip),
            identical(stopped$bound_status, c("ok", "not_evaluated", "not_evaluated")))
  final_only <- raw
  final_only$information <- c(95, 100, 110)
  final_only$info_blinded_final_projected[3] <- NA_real_
  final_only <- evaluate_interim_skip(final_only, design, "information")
  stopifnot(identical(final_only$bound_status, c("skipped", "skipped", "ok")),
            isTRUE(all.equal(final_only$upper_bound[3], stats::qnorm(0.975))))
  invalid <- raw
  invalid$info_blinded_final_projected[2] <- NA_real_
  invalid <- evaluate_interim_skip(invalid, design, "information")
  stopifnot(!invalid$near_final_skip[2], invalid$bound_status[2] == "failed")
  pre_skipped <- raw
  pre_skipped$status[2] <- "skipped"
  stopifnot(!any(evaluate_interim_skip(pre_skipped, design, "information")$near_final_skip))
  multiple <- rbind(raw, transform(below, sim = 2L))
  multiple <- evaluate_interim_skip(multiple, design, "information")
  stopifnot(sum(multiple$near_final_skip) == 1L)
  for (threshold in c(0, -0.1, 1.1, NA_real_)) {
    stopifnot(inherits(try(evaluate_interim_skip(raw, design, "information", threshold),
                          silent = TRUE), "try-error"))
  }
  cat("PASS: threshold inclusion, no spending at skipped looks, prior boundary preservation, final date, full final spending, prior stops, missing projections, existing skips, multiple trials, and input validation.\n")
}

summarize_interim_skip <- function(checked, cell, metric, threshold, ignore_futility, final_look) {
  bound_status <- near_final_skip <- NULL
  paths <- data.table::as.data.table(checked)[, list(
    reject = any(cross_upper), futility = any(cross_lower),
    failed = any(bound_status %in% c("failed", "missing_information_or_statistic")),
    skipped = any(near_final_skip),
    final_evaluated = any(analysis == final_look & bound_status == "ok")
  ), by = "sim"]
  n_paths <- nrow(paths)
  skip_rate <- mean(paths$skipped)
  interval <- stats::binom.test(sum(paths$reject), n_paths)$conf.int
  data.frame(
    cell, information_estimator = metric, threshold = threshold,
    ignore_futility = ignore_futility, n_paths = n_paths,
    n_skipped_paths = sum(paths$skipped), skip_rate = skip_rate,
    skip_mcse = sqrt(skip_rate * (1 - skip_rate) / n_paths),
    n_skipped_looks = sum(checked$near_final_skip),
    n_final_evaluated = sum(paths$final_evaluated),
    n_failed = sum(paths$failed), failure_rate = mean(paths$failed),
    rejection_rate = mean(paths$reject), ci_lower = interval[1], ci_upper = interval[2],
    rejection_or_failure = mean(paths$reject | paths$failed),
    futility_rate = mean(paths$futility)
  )
}

run_interim_skip_study <- function(args) {
  bound_status <- near_final_skip <- NULL
  option <- function(name, default) {
    values <- sub(paste0("^--", name, "="), "", args[startsWith(args, paste0("--", name, "="))])
    if (length(values)) values[1] else default
  }
  thresholds <- as.numeric(strsplit(option("thresholds", "0.95"), ",", fixed = TRUE)[[1]])
  metrics <- strsplit(option("metrics", "info_blinded_ml,info_unblinded_ml"), ",", fixed = TRUE)[[1]]
  max_sims <- as.numeric(option("max-sims", "Inf"))
  stopifnot(length(thresholds) > 0L, all(is.finite(thresholds)),
            all(thresholds > 0 & thresholds <= 1),
            all(metrics %in% c("info_blinded_ml", "info_unblinded_ml")),
            length(max_sims) == 1L, !is.na(max_sims), max_sims >= 1,
            is.infinite(max_sims) || max_sims == floor(max_sims))
  source_files <- c("R/gsNBPower.R", "R/check_gs_bound.R", "R/sim_gs_nbinom.R",
                    "R/blinded_ssr.R", "R/get_cut_date.R", "data-raw/gs_strategy_study.R")
  source_signature <- tools::md5sum(source_files)
  evaluation_signature <- tools::md5sum("data-raw/evaluate_interim_skip.R")
  mode <- if (is.infinite(max_sims)) "production" else paste0("subset-", max_sims)
  output_dir <- file.path("paper", "gs-strategy-results", "interim-skip", mode)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- list.files(file.path("paper", "gs-strategy-results", "cells"),
                      pattern = "\\.rds$", full.names = TRUE)
  results <- list()
  for (path in paths) {
    saved <- readRDS(path)
    expected_n <- if (saved$cell$truth == "null") 20000L else 5000L
    if (!identical(saved$source_signature, source_signature) || saved$cell$n_sims != expected_n ||
        saved$cell$cut_rule == "ssr") next
    filters <- c("scenario", "truth", "timing_method", "tier")
    selected <- vapply(filters, function(name) {
      value <- option(gsub("_", "-", name), "all")
      value == "all" || saved$cell[[name]] == value
    }, logical(1))
    if (!all(selected)) next
    cache_path <- file.path(output_dir, basename(path))
    settings <- list(thresholds = thresholds, metrics = metrics, max_sims = max_sims)
    if (file.exists(cache_path)) {
      cached <- readRDS(cache_path)
      if (identical(cached$source_signature, source_signature) &&
          identical(cached$evaluation_signature, evaluation_signature) &&
          identical(cached$settings, settings)) {
        results[[length(results) + 1L]] <- cached$summary
        cat("Reused:", basename(path), "\n")
        next
      }
    }
    raw <- as.data.frame(saved$raw)
    ids <- unique(raw$sim)
    if (is.finite(max_sims)) ids <- head(ids, max_sims)
    raw <- raw[raw$sim %in% ids, , drop = FALSE]
    rows <- looks <- list()
    for (metric in metrics) {
      for (threshold in thresholds) {
        modes <- if (saved$design$test.type == 1L) FALSE else c(FALSE, TRUE)
        for (ignore_futility in modes) {
          checked <- evaluate_interim_skip(raw, saved$design, metric, threshold, ignore_futility)
          row <- summarize_interim_skip(checked, saved$cell, metric, threshold,
                                         ignore_futility, saved$design$k)
          rows[[length(rows) + 1L]] <- row
          per_look <- data.table::as.data.table(checked)[, list(
            n_paths = .N, n_skipped = sum(near_final_skip),
            n_evaluated = sum(bound_status == "ok"),
            n_failed = sum(bound_status %in% c("failed", "missing_information_or_statistic"))
          ), by = "analysis"]
          looks[[length(looks) + 1L]] <- cbind(
            saved$cell[rep(1, nrow(per_look)), ], information_estimator = metric,
            threshold = threshold, ignore_futility = ignore_futility, as.data.frame(per_look)
          )
        }
      }
    }
    summary <- do.call(rbind, rows)
    payload <- list(summary = summary, look_summary = do.call(rbind, looks),
                    settings = settings, mode = mode,
                    source_signature = source_signature, evaluation_signature = evaluation_signature,
                    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
    saveRDS(payload, paste0(cache_path, ".tmp"))
    if (!file.rename(paste0(cache_path, ".tmp"), cache_path)) stop("Could not save reanalysis checkpoint.")
    results[[length(results) + 1L]] <- summary
    cat("Evaluated:", basename(path), "using", length(ids), "saved trials\n")
  }
  if (!length(results)) stop("No compatible checkpoints matched the requested settings.")
  summary <- do.call(rbind, results)
  payload <- list(summary = summary, mode = mode, settings = settings,
                  source_signature = source_signature, evaluation_signature = evaluation_signature,
                  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  saveRDS(payload, file.path(output_dir, "summary.rds"))
  columns <- c("tier", "scenario", "timing_method", "truth", "nuisance", "test_type",
               "cut_rule", "information_estimator", "threshold", "ignore_futility", "n_paths",
               "n_skipped_paths", "skip_rate", "skip_mcse", "n_final_evaluated", "n_failed",
               "rejection_rate", "ci_lower", "ci_upper")
  report <- summary[columns]
  for (column in c("skip_rate", "rejection_rate", "ci_lower", "ci_upper")) {
    report[[column]] <- sprintf("%.3f%%", 100 * report[[column]])
  }
  report$skip_mcse <- sprintf("%.3f pp", 100 * report$skip_mcse)
  writeLines(c(
    "# Skip near-final interims, retain the planned final analysis", "",
    paste("Generated:", payload$generated_at), paste("Mode:", mode), "",
    "At an active interim, skip formal testing when uncapped current information / blinded projected-final information is at least the threshold. Spend no alpha at that look. Preserve prior spending, retain the scheduled final analysis, and complete spending there using observed information history. This does not promote an interim to the final analysis.", "",
    "Counts apply the revised rule sequentially to saved trial histories, including prior efficacy/futility stops; candidates after a valid stop are not counted as skips. Existing skipped looks are not counted as near-final skips. Adaptive SSR cells are excluded. Both current-information estimators are labeled; the projected-final estimate is blinded. No new trials were generated. Unresolved analysis failures remain explicitly reported, so rejection rates alone do not establish calibration.", "",
    "n_paths is the number of saved trials actually evaluated, not necessarily the checkpoint replication count n_sims. A subset-* mode is only a smoke test, not a production calibration result. skip_mcse is the Monte Carlo standard error of the trial-level skip rate. ci_lower and ci_upper give the exact 95% binomial interval for rejection_rate. Per-look counts and provenance are stored in the individual RDS checkpoints.", "",
    paste0("| ", paste(columns, collapse = " | "), " |"),
    paste0("| ", paste(rep("---", length(columns)), collapse = " | "), " |"),
    apply(report, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  ), file.path(output_dir, "summary.md"))
  cat("Report:", file.path(output_dir, "summary.md"), "\n")
  invisible(payload)
}

if (sys.nframe() == 0L) {
  devtools::load_all(".", quiet = TRUE)
  args <- commandArgs(trailingOnly = TRUE)
  if ("--self-test" %in% args) check_interim_skip() else run_interim_skip_study(args)
}