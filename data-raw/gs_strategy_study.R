#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
production <- "--production" %in% args
if (production == ("--debug" %in% args)) stop("Specify --debug or --production.")
option <- function(name, default) {
  values <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(values) > 1L) stop("Duplicate option: ", name)
  if (length(values)) sub(paste0("^--", name, "="), "", values) else default
}
output_dir <- option("output-dir", file.path("paper", "gs-strategy-results"))
workers <- as.integer(option("workers", "1"))
if (!is.finite(workers) || workers < 1L) stop("workers must be a positive integer.")
tiers <- strsplit(option("tiers", "A,B,C"), ",", fixed = TRUE)[[1]]
if (!all(tiers %in% c("A", "B", "C"))) stop("tiers must be A, B, or C.")
devtools::load_all(".", quiet = TRUE)
data.table::setDTthreads(1L)
source_signature <- tools::md5sum(c(
  "R/gsNBPower.R", "R/check_gs_bound.R", "R/sim_gs_nbinom.R",
  "R/blinded_ssr.R", "R/get_cut_date.R", "data-raw/gs_strategy_study.R"
))

scenarios <- data.frame(
  scenario = c("reference", "low-rate-ramp", "high-dispersion", "low-rate-gap",
               "moderate-gap", "low-rate-front-loaded", "high-rate-front-loaded", "mixed-stress"),
  lambda1 = c(0.5, 0.15, 0.5, 0.15, 0.5, 0.15, 1, 0.3),
  dispersion = c(0.2, 0.2, 1, 1, 0.2, 1, 0.2, 0.5),
  gap_days = c(0, 0, 0, 30, 30, 0, 15, 15),
  dropout_annual = c(0, 0.1, 0.1, 0, 0.1, 0.1, 0, 0.05),
  accrual_pattern = c("constant", "ramp", "ramp", "constant", "constant",
                      "front_loaded", "front_loaded", "ramp")
)
shape <- function(pattern) switch(pattern,
  constant = list(rate = 1, duration = 12),
  ramp = list(rate = c(1, 2), duration = c(6, 6)),
  front_loaded = list(rate = c(2, 1), duration = c(6, 6))
)
make_design <- function(scenario, timing) {
  enrollment <- shape(scenario$accrual_pattern)
  fixed <- sample_size_nbinom(
    lambda1 = scenario$lambda1, lambda2 = scenario$lambda1 * 0.7,
    dispersion = scenario$dispersion, power = 0.9, alpha = 0.025,
    accrual_rate = enrollment$rate, accrual_duration = enrollment$duration,
    trial_duration = 24, max_followup = 12,
    dropout_rate = -log1p(-scenario$dropout_annual) / 12,
    event_gap = scenario$gap_days / 30.4375, test_type = "score"
  )
  design <- if (timing == "information_fraction") {
    gsNB(fixed, test.type = 1)
  } else gsNBCalendar(fixed, test.type = 1, analysis_times = c(10, 18, 24))
  toInteger(design)
}
make_cell <- function(tier, scenario_index, timing, truth, nuisance = "planned",
                      test_type = "score", cut_rule = "calendar") {
  data.frame(tier, scenario_index, timing_method = timing, truth, nuisance,
             test_type, cut_rule, stringsAsFactors = FALSE)
}
cells <- list()
if ("A" %in% tiers) {
  for (scenario_index in seq_len(nrow(scenarios))) {
    for (timing in c("information_fraction", "fixed_calendar")) {
      for (truth in c("null", "alternative")) {
        cells[[length(cells) + 1L]] <- make_cell("A", scenario_index, timing, truth)
      }
    }
  }
}
if ("B" %in% tiers) {
  for (scenario_index in c(1L, 3L, 4L, 8L)) {
    for (truth in c("null", "alternative")) {
      for (nuisance in c("shortfall", "excess", "slow_accrual")) {
        for (test_type in c("score", "wald")) {
          cells[[length(cells) + 1L]] <- make_cell(
            "B", scenario_index, "information_fraction", truth, nuisance, test_type
          )
        }
      }
    }
  }
}
if ("C" %in% tiers) {
  for (truth in c("null", "alternative")) {
    for (cut_rule in c("minimum_age", "blinded_information", "ssr")) {
      for (timing in c("information_fraction", "futility")) {
        cells[[length(cells) + 1L]] <- make_cell("C", 1L, timing, truth,
                                              cut_rule = cut_rule)
      }
    }
  }
}
cells <- do.call(rbind, cells)
cells$scenario <- scenarios$scenario[cells$scenario_index]
selected_scenario <- option("scenario", "")
if (nzchar(selected_scenario)) cells <- cells[cells$scenario == selected_scenario, , drop = FALSE]
if (!nrow(cells)) stop("No cells match the selected scenario.")
cells$n_sims <- if (production) ifelse(cells$truth == "null", 20000L, 5000L) else 30L
cells$seed <- 20261106L + cells$scenario_index * 1000L +
  ifelse(cells$timing_method == "information_fraction", 100L, 200L) +
  ifelse(cells$truth == "null", 1L, 2L) +
  ifelse(cells$tier == "A", 0L,
         100000L * match(cells$nuisance, c("planned", "shortfall", "excess", "slow_accrual")) +
           10000L * match(cells$cut_rule, c("calendar", "minimum_age", "blinded_information", "ssr")))
if (!production && !nzchar(selected_scenario)) {
  keep <- c(which(cells$tier == "A" & cells$scenario_index == 1L),
            which(cells$tier == "B" & cells$scenario_index == 1L & cells$nuisance == "shortfall"),
            which(cells$tier == "C" & cells$timing_method == "futility"))
  cells <- cells[keep, , drop = FALSE]
}
cell_index <- as.integer(option("cell-index", "0"))
if (!is.finite(cell_index) || cell_index < 0L || cell_index > nrow(cells)) {
  stop("cell-index must identify an available cell, or be zero for all cells.")
}
if (cell_index > 0L) cells <- cells[cell_index, , drop = FALSE]
cell_dir <- file.path(output_dir, "cells")
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)

summarize_paths <- function(checked, prediction, cell, policy, spending, metric, ignore_futility) {
  groups <- split(seq_len(nrow(checked)), checked$sim)
  outcomes <- lapply(groups, function(rows) {
    rows <- rows[order(checked$analysis[rows])]
    crosses <- checked$cross_upper[rows] | checked$cross_lower[rows] | checked$cross_harm[rows]
    stop_row <- if (any(crosses)) rows[which(crosses)[1]] else utils::tail(rows, 1)
    failed <- if ("bound_status" %in% names(checked)) {
      any(checked$bound_status[rows] %in% c("failed", "missing_information_or_statistic"))
    } else any(!is.finite(checked$z_stat[utils::tail(rows, 1)]))
    data.frame(reject = any(checked$cross_upper[rows]),
               futility = any(checked$cross_lower[rows]), failed = failed,
               n = checked$n_enrolled[stop_row], events = checked$events_total[stop_row],
               time = checked$analysis_time[stop_row])
  })
  outcomes <- do.call(rbind, outcomes)
  rejection <- mean(outcomes$reject)
  ci <- stats::binom.test(sum(outcomes$reject), nrow(outcomes))$conf.int
  predicted <- if (is.null(prediction)) NA_real_ else prediction$power
  data.frame(
    cell, policy = policy, spending = spending, information_estimator = metric,
    ignore_futility = ignore_futility, rejection_rate = rejection,
    mcse = sqrt(rejection * (1 - rejection) / nrow(outcomes)),
    ci_lower = ci[1], ci_upper = ci[2], predicted_power = predicted,
    power_difference = rejection - predicted,
    n_failed = sum(outcomes$failed),
    rejection_or_failure = mean(outcomes$reject | outcomes$failed),
    futility_rate = mean(outcomes$futility),
    mean_n = mean(outcomes$n), mean_events = mean(outcomes$events),
    mean_time = mean(outcomes$time),
    predicted_n = if (is.null(prediction)) NA_real_ else prediction$expected_n,
    predicted_events = if (is.null(prediction)) NA_real_ else prediction$expected_events,
    predicted_time = if (is.null(prediction)) NA_real_ else prediction$expected_analysis_time
  )
}

run_cell <- function(cell) {
  key <- paste(c(unlist(cell[c("tier", "scenario", "timing_method", "truth", "nuisance",
                              "test_type", "cut_rule")]), cell$n_sims), collapse = "__")
  cell_path <- file.path(cell_dir, paste0("strategy-v2__", key, ".rds"))
  if (file.exists(cell_path)) {
    cached <- readRDS(cell_path)
    if (identical(cached$source_signature, source_signature)) {
      message(key, " | cached")
      return(cached)
    }
  }
  message(key, " | simulate")
  scenario <- scenarios[cell$scenario_index, ]
  design <- if (cell$timing_method == "futility") {
    cache <- readRDS("inst/extdata/gs_simulation_results.rds")
    cache$gs_nb
  } else make_design(scenario, cell$timing_method)
  inputs <- design$nb_design$inputs
  true_lambda1 <- inputs$lambda1 * switch(cell$nuisance, shortfall = 0.5, excess = 2, 1)
  true_dispersion <- inputs$dispersion * if (cell$nuisance == "shortfall") 2 else 1
  rr <- if (cell$truth == "null") 1 else inputs$lambda2 / inputs$lambda1
  enrollment <- .gs_nb_accrual(design)
  if (cell$nuisance == "slow_accrual") {
    enrollment$rate <- enrollment$rate * 0.75
    enrollment$duration <- enrollment$duration / 0.75
  }
  n_target <- as.integer(ceiling(utils::tail(design$n_total, 1)))
  if (cell$tier == "A") {
    historical_shape <- shape(scenario$accrual_pattern)
    enrollment$rate <- historical_shape$rate * n_target /
      sum(historical_shape$rate * historical_shape$duration)
    enrollment$duration <- historical_shape$duration
  }
  fail_rate <- data.frame(treatment = c("Control", "Experimental"),
                          rate = c(true_lambda1, true_lambda1 * rr),
                          dispersion = rep(true_dispersion, length.out = 2))
  dropout <- inputs$dropout_rate
  dropout_scalar <- if (is.numeric(dropout)) rep(dropout, length.out = 2) else dropout$rate[1]
  dropout_rate <- data.frame(treatment = c("Control", "Experimental"),
                             rate = dropout_scalar, duration = if (cell$tier == "A") 100 else Inf)
  prediction_args <- list(
    design = design, lambda1 = true_lambda1, lambda2 = true_lambda1 * rr,
    dispersion = true_dispersion, accrual_rate = enrollment$rate,
    accrual_duration = enrollment$duration, test_type = cell$test_type
  )
  cuts <- lapply(design$T, function(time) list(planned_calendar = time))
  if (cell$cut_rule == "minimum_age") {
    cuts[[1]] <- list(min_n = floor(n_target / 2), min_followup = 3)
    prediction_args$planned_calendar <- c(NA, design$T[-1])
    prediction_args$min_n <- c(floor(n_target / 2), NA, NA)
    prediction_args$min_followup <- c(3, NA, NA)
  }
  if (cell$cut_rule == "blinded_information") {
    cuts <- c(lapply(design$n.I[-design$k], function(info) list(target_info = info)),
              list(list(planned_calendar = design$T[design$k])))
    prediction_args$planned_calendar <- c(rep(NA_real_, design$k - 1L), design$T[design$k])
    prediction_args$target_info <- c(design$n.I[-design$k], NA)
  }
  if (cell$cut_rule == "ssr") {
    raw <- sim_ssr_nbinom(
      n_sims = cell$n_sims, enroll_rate = data.frame(rate = n_target / 12, duration = 24),
      fail_rate = fail_rate, dropout_rate = dropout_rate, max_followup = 12,
      design = design, n_max = ceiling(1.5 * n_target),
      strategies = c("No adaptation", "Blinded SSR"),
      bound_info = "blinded_ml", ignore_futility = cell$truth == "null",
      test_type = cell$test_type, seed = cell$seed
    )
    summary <- summarize_ssr_sim(raw)$trial_summary
    summary <- cbind(cell[rep(1, nrow(summary)), ], summary)
    result <- list(cell = cell, design = design, raw = raw, ssr_summary = summary)
  } else {
    raw <- sim_gs_nbinom(
      n_sims = cell$n_sims, enroll_rate = data.frame(rate = enrollment$rate, duration = enrollment$duration),
      fail_rate = fail_rate, dropout_rate = dropout_rate, max_followup = 12,
      event_gap = inputs$event_gap, cuts = cuts, n_target = n_target,
      design = design, test_type = cell$test_type, seed = cell$seed
    )
    configurations <- data.frame(policy = "fixed", spending = "min_planned_actual",
                                  metric = "info_unblinded_ml")
    dynamic <- expand.grid(
      policy = "recompute",
      spending = c("min_planned_actual", "planned_final", "expected_final", "calendar"),
      metric = c("info_blinded_ml", "info_unblinded_ml"), stringsAsFactors = FALSE
    )
    configurations <- rbind(configurations, dynamic)
    rows <- list()
    looks <- list()
    for (config_index in seq_len(nrow(configurations))) {
      config <- configurations[config_index, ]
      modes <- if (design$test.type == 1L) FALSE else c(FALSE, TRUE)
      for (ignore_futility in modes) {
        prediction <- tryCatch(do.call(gsNBPower, c(prediction_args, list(
          boundary_update = config$policy, spending = config$spending
        ))), error = function(error) NULL)
        if (!is.null(prediction) && ignore_futility) {
          active <- prediction$look_results$status == "analyzed"
          stage <- prediction$look_results[active, ]
          predicted <- .gs_nb_probabilities(
            stage$information,
            if (cell$test_type == "score") stage$information_null else stage$information,
            log(1 / rr),
            stage$upper_bound, rep(-Inf, nrow(stage)), stage$harm_bound
          )
          prediction$power <- sum(predicted$upper)
          mass <- predicted$upper + predicted$harm
          mass[length(mass)] <- mass[length(mass)] + utils::tail(predicted$continuation, 1)
          prediction$expected_n <- sum(mass * stage$n_total)
          prediction$expected_events <- sum(mass * stage$expected_events)
          prediction$expected_analysis_time <- sum(mass * stage$analysis_time)
        }
        checked <- check_gs_bound(
          raw, design, info_col = config$metric,
          boundary_update = config$policy, spending = config$spending,
          ignore_futility = ignore_futility
        )
        rows[[length(rows) + 1L]] <- summarize_paths(
          checked, prediction, cell, config$policy, config$spending, config$metric, ignore_futility
        )
        per_look <- lapply(split(checked, checked$analysis), function(look) {
          fractions <- look$info_blinded_ml / look$info_blinded_final_projected
          valid <- fractions[is.finite(fractions)]
          expected_fraction <- if (is.null(prediction)) NA_real_ else {
            prediction$look_results$information[look$analysis[1]] /
              utils::tail(prediction$look_results$information, 1)
          }
          data.frame(analysis = look$analysis[1], mean_time = mean(look$analysis_time),
                     prob_upper = mean(look$cross_upper), prob_lower = mean(look$cross_lower),
                     n_failed = sum(look$bound_status == "failed"),
                     n_skipped = sum(look$status == "skipped"),
                     mean_projected_final = mean(look$info_blinded_final_projected, na.rm = TRUE),
                     mean_information = mean(look[[config$metric]], na.rm = TRUE),
                     mean_fraction = if (length(valid)) mean(valid) else NA_real_,
                     sd_fraction = if (length(valid) > 1L) stats::sd(valid) else NA_real_,
                     expected_fraction = expected_fraction,
                     fraction_bias = if (length(valid)) mean(valid) - expected_fraction else NA_real_,
                     n_projection_failed = sum(!is.finite(fractions)),
                     n_mom_fallback = sum(look$blinded_fallback == "mom", na.rm = TRUE),
                     n_poisson_fallback = sum(look$blinded_fallback == "poisson", na.rm = TRUE))
        })
        looks[[length(looks) + 1L]] <- cbind(
          cell[rep(1, design$k), ], policy = config$policy, spending = config$spending,
          information_estimator = config$metric, ignore_futility = ignore_futility,
          do.call(rbind, per_look)
        )
      }
    }
    legacy <- check_gs_bound(raw, design, info_col = "info_unblinded_ml")
    rows[[length(rows) + 1L]] <- summarize_paths(
      legacy, NULL, cell, "legacy", "planned_final", "info_unblinded_ml", FALSE
    )
    result <- list(cell = cell, design = design, raw = raw,
                   summary = do.call(rbind, rows), look_summary = do.call(rbind, looks))
  }
  result$source_signature <- source_signature
  saveRDS(result, paste0(cell_path, ".tmp"), compress = "gzip")
  if (!file.rename(paste0(cell_path, ".tmp"), cell_path)) stop("Could not save cell: ", key)
  result
}

run_study <- function() {
  previous <- future::plan()
  on.exit(future::plan(previous), add = TRUE)
  if (workers > 1L) {
    if (.Platform$OS.type == "windows") future::plan(future::multisession, workers = workers) else
      future::plan(future::multicore, workers = workers)
  } else future::plan(future::sequential)
  result <- future.apply::future_lapply(
    seq_len(nrow(cells)), function(index) run_cell(cells[index, ]), future.seed = TRUE
  )
  summary <- do.call(rbind, lapply(result, function(cell) cell$summary))
  look_summary <- do.call(rbind, lapply(result, function(cell) cell$look_summary))
  ssr_summary <- do.call(rbind, lapply(result, function(cell) cell$ssr_summary))
  legacy_path <- "paper/gs-validation-results/operating-characteristics.csv"
  legacy_comparison <- NULL
  if (file.exists(legacy_path) && !is.null(summary)) {
    historical <- utils::read.csv(legacy_path)
    current <- summary[summary$tier == "A" & summary$policy == "legacy", ]
    legacy_comparison <- merge(historical, current,
                               by = c("scenario", "timing_method", "truth"),
                               suffixes = c("_historical", "_current"))
  }
  payload <- list(summary = summary, look_summary = look_summary,
                  ssr_summary = ssr_summary, scenarios = scenarios, cells = cells,
                  legacy_comparison = legacy_comparison,
                  mode = if (production) "production" else "debug",
                  generated_at = as.character(Sys.time()), session_info = utils::sessionInfo())
  payload$source_signature <- source_signature
  saveRDS(payload, file.path(output_dir, "gs_strategy_summary.rds"), compress = "xz")
  if (!is.null(summary)) utils::write.csv(summary, file.path(output_dir, "strategy-summary.csv"), row.names = FALSE)
  if ("--publish" %in% args) saveRDS(payload, "inst/extdata/gs_strategy_summary.rds", compress = "xz")
  print(summary[c("tier", "scenario", "truth", "spending", "rejection_rate", "n_failed")])
  message("Saved ", nrow(cells), " cells under ", normalizePath(output_dir))
}
run_study()