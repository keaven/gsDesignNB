# Generate production SSR simulation outputs for the score-test workflow.
#
# This script mirrors the SSR simulation-study vignette, but runs outside
# pkgdown so future workers can use the installed package namespace directly.
# It checkpoints scenario and chunk output under /tmp by default, then writes
# the compact cache consumed by vignettes/ssr-simulation-study.Rmd.
#
# Example:
#   Rscript data-raw/generate_ssr_score_outputs.R
#   GSDESIGNNB_PRODUCTION_SSR=true GSDESIGNNB_FORCE_SSR_SIM=true \
#     Rscript data-raw/generate_ssr_score_outputs.R

suppressPackageStartupMessages({
  library(gsDesignNB)
  library(gsDesign)
  library(data.table)
  library(future)
})

is_true <- function(x) identical(tolower(Sys.getenv(x)), "true")

project_root <- if (file.exists("DESCRIPTION")) "." else
  if (file.exists("../DESCRIPTION")) ".." else "."

production <- is_true("GSDESIGNNB_PRODUCTION_SSR")
force_main <- is_true("GSDESIGNNB_FORCE_SSR_SIM")
force_type1 <- is_true("GSDESIGNNB_FORCE_TYPE1_SIM")

analysis_test_type <- "score"
analysis_test_label <- tools::toTitleCase(analysis_test_type)

env_int <- function(name, default) {
  value <- as.integer(Sys.getenv(name, as.character(default)))
  if (!is.finite(value) || value < 1L) default else value
}

n_sims_initial <- env_int("GSDESIGNNB_SSR_N_INITIAL", 50L)
n_sims_production_power <- env_int("GSDESIGNNB_SSR_N_POWER", 3600L)
n_sims_production_type1 <- env_int("GSDESIGNNB_SSR_N_TYPE1", 20000L)
n_sims_production_rr_gt1 <- env_int("GSDESIGNNB_SSR_N_RR_GT1", 1000L)

chunk_size <- as.integer(Sys.getenv("GSDESIGNNB_SSR_CHUNK_SIZE", "1000"))
if (!is.finite(chunk_size) || chunk_size < 1L) chunk_size <- 1000L

default_workers <- max(1L, future::availableCores() - 1L)
workers <- as.integer(Sys.getenv("GSDESIGNNB_WORKERS", as.character(default_workers)))
if (!is.finite(workers) || workers < 1L) workers <- default_workers

cache_dir <- Sys.getenv(
  "GSDESIGNNB_SSR_CACHE_DIR",
  file.path(tempdir(), "gsDesignNB_ssr_score_cache")
)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
cached_scenarios_at_start <- sum(file.exists(file.path(
  cache_dir,
  sprintf("scenario_%03d.rds", seq_len(90L))
)))

extdata_dir <- file.path(project_root, "inst", "extdata")
dir.create(extdata_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(extdata_dir, "ssr_sim_vignette_outputs_score.rds")

lambda1_plan <- 0.5
rr_plan <- 0.7
lambda2_plan <- lambda1_plan * rr_plan
k_plan <- 0.5
power_plan <- 0.9
alpha_plan <- 0.025
accrual_rate_plan <- 30
accrual_scenario_plan <- 18
accrual_dur_plan <- 12
max_followup <- 12
trial_dur_plan <- accrual_dur_plan + max_followup
event_gap_val <- 20 / 30.4375
analysis_times_plan <- c(9, 14, 24)

design_plan <- sample_size_nbinom(
  lambda1 = lambda1_plan, lambda2 = lambda2_plan,
  dispersion = k_plan, power = power_plan, alpha = alpha_plan,
  accrual_rate = accrual_rate_plan,
  accrual_duration = accrual_dur_plan,
  trial_duration = trial_dur_plan,
  max_followup = max_followup,
  event_gap = event_gap_val,
  test_type = analysis_test_type
)

gs_plan <- design_plan |>
  gsNBCalendar(
    k = 3, test.type = 4, alpha = alpha_plan,
    sfu = sfHSD, sfupar = -2,
    sfl = sfCauchy, sflpar = c(0.4, 0.75, 0.56, 0.63),
    analysis_times = analysis_times_plan
  ) |>
  gsDesignNB::toInteger()

n_planned <- gs_plan$n_total[gs_plan$k]
planned_timing <- gs_plan$timing
n_max <- 2 * n_planned
min_if_futility <- 0.3
max_enrollment_frac_for_ia2 <- 1.00
min_months_to_close_for_adapt <- 2
analysis_lag_months <- 2
gs_inflation <- n_planned / design_plan$n_total
accrual_rate_plan_eff <- n_planned / accrual_dur_plan
design_note <- paste0(
  "Design: lambda1=", lambda1_plan,
  ", RR=", rr_plan,
  ", k=", k_plan,
  ", planned accrual=", round(accrual_rate_plan_eff, 1), "/mo",
  ", planned N=", n_planned,
  ", max follow-up=", max_followup, " mo"
)

scenarios <- expand.grid(
  lambda1_true = c(0.3, 0.5, 0.8),
  k_true = c(0.5, 1.0),
  accrual_true = c(12, 18, 24),
  rr_true = c(0.6, 0.7, 0.85, 1.0, 1.1),
  stringsAsFactors = FALSE
)
scenarios$n_sims <- if (production) {
  as.integer(ifelse(
    scenarios$rr_true == 1,
    n_sims_production_type1,
    ifelse(scenarios$rr_true > 1, n_sims_production_rr_gt1, n_sims_production_power)
  ))
} else {
  rep(as.integer(n_sims_initial), nrow(scenarios))
}
scenarios$accrual_eff <- scenarios$accrual_true

max_scenarios <- as.integer(Sys.getenv("GSDESIGNNB_SSR_MAX_SCENARIOS", "0"))
if (is.finite(max_scenarios) && max_scenarios > 0L) {
  scenarios <- scenarios[seq_len(min(max_scenarios, nrow(scenarios))), , drop = FALSE]
}

make_enroll_rate <- function(accrual_eff) {
  data.frame(rate = accrual_eff, duration = n_max / accrual_eff)
}

make_fail_rate <- function(lambda1_true, rr_true, k_true) {
  data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(lambda1_true, lambda1_true * rr_true),
    dispersion = k_true
  )
}

dropout_rate_sim <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0, 0),
  duration = c(100, 100)
)

split_chunks <- function(n, size) {
  starts <- seq.int(1L, n, by = size)
  ends <- pmin(starts + size - 1L, n)
  data.table(chunk = seq_along(starts), start = starts, end = ends, n = ends - starts + 1L)
}

cache_key_num <- function(x) {
  gsub("[^0-9A-Za-z]+", "p", format(x, trim = TRUE, scientific = FALSE))
}

run_ssr_chunk <- function(sc_idx, sc, n_chunk, sim_offset, seed) {
  sim_res <- sim_ssr_nbinom(
    n_sims = n_chunk,
    enroll_rate = make_enroll_rate(sc$accrual_eff),
    fail_rate = make_fail_rate(sc$lambda1_true, sc$rr_true, sc$k_true),
    dropout_rate = dropout_rate_sim,
    max_followup = max_followup,
    design = gs_plan,
    n_max = n_max,
    min_if_futility = min_if_futility,
    max_enrollment_frac_for_adapt = max_enrollment_frac_for_ia2,
    min_months_to_close_for_adapt = min_months_to_close_for_adapt,
    analysis_lag_months = analysis_lag_months,
    event_gap = event_gap_val,
    ignore_futility = FALSE,
    test_type = analysis_test_type,
    metadata = list(
      lambda1_true = sc$lambda1_true,
      k_true = sc$k_true,
      accrual_true = sc$accrual_true,
      accrual_eff = sc$accrual_eff,
      analysis_test = analysis_test_type,
      rr_true = sc$rr_true
    ),
    seed = seed
  )

  out <- as.data.table(sim_res$trial_results)
  out[, sim := sim + sim_offset]
  out[, scenario_index := sc_idx]
  out[]
}

run_scenario <- function(sc_idx) {
  sc <- scenarios[sc_idx, ]
  scenario_path <- file.path(cache_dir, sprintf("scenario_%03d.rds", sc_idx))
  if (!force_main && file.exists(scenario_path)) {
    message(sprintf("[%s] Loading scenario %03d cache", Sys.time(), sc_idx))
    return(readRDS(scenario_path))
  }

  message(sprintf(
    "[%s] Scenario %03d/%03d: RR=%.2f lambda1=%.2f k=%.1f accrual=%.0f n=%d",
    Sys.time(), sc_idx, nrow(scenarios), sc$rr_true, sc$lambda1_true,
    sc$k_true, sc$accrual_true, sc$n_sims
  ))

  chunks <- split_chunks(sc$n_sims, chunk_size)
  chunk_results <- vector("list", nrow(chunks))
  for (j in seq_len(nrow(chunks))) {
    chunk_path <- file.path(cache_dir, sprintf("scenario_%03d_chunk_%03d.rds", sc_idx, j))
    if (!force_main && file.exists(chunk_path)) {
      chunk_results[[j]] <- readRDS(chunk_path)
      next
    }
    message(sprintf(
      "  [%s] chunk %03d/%03d: sims %d-%d",
      Sys.time(), j, nrow(chunks), chunks$start[j], chunks$end[j]
    ))
    chunk_results[[j]] <- run_ssr_chunk(
      sc_idx = sc_idx,
      sc = sc,
      n_chunk = chunks$n[j],
      sim_offset = chunks$start[j] - 1L,
      seed = 1000000L + sc_idx * 1000L + j
    )
    saveRDS(chunk_results[[j]], chunk_path)
  }

  scenario_dt <- rbindlist(chunk_results, use.names = TRUE, fill = TRUE)
  saveRDS(scenario_dt, scenario_path)
  scenario_dt
}

run_null_type1 <- function(test_label, test_type_x, null_n) {
  key <- tolower(test_label)
  cache_path <- file.path(extdata_dir, paste0("ssr_type1_null_alpha025_", key, ".rds"))
  if (!force_type1 && file.exists(cache_path)) {
    cr <- readRDS(cache_path)
    if (isTRUE(as.integer(cr$null_nonbinding_n) == as.integer(null_n))) {
      return(list(summary = as.data.table(cr$summary), runtime_seconds = cr$runtime_seconds))
    }
  }

  t0 <- Sys.time()
  null_k_scenarios <- c(k_plan, 1.0)
  all_null <- vector("list", length(null_k_scenarios))
  for (i in seq_along(null_k_scenarios)) {
    k_null <- null_k_scenarios[i]
    message(sprintf(
      "[%s] Type I %s test: k=%.1f n=%d",
      Sys.time(), test_label, k_null, null_n
    ))
    chunks <- split_chunks(null_n, chunk_size)
    k_chunks <- vector("list", nrow(chunks))
    for (j in seq_len(nrow(chunks))) {
      chunk_path <- file.path(
        cache_dir,
        sprintf("type1_%s_k%s_chunk_%03d.rds", key, cache_key_num(k_null), j)
      )
      if (!force_type1 && file.exists(chunk_path)) {
        message(sprintf(
          "  [%s] Loading Type I %s k=%.1f chunk %03d/%03d cache",
          Sys.time(), test_label, k_null, j, nrow(chunks)
        ))
        k_chunks[[j]] <- readRDS(chunk_path)
        next
      }
      message(sprintf(
        "  [%s] Type I %s k=%.1f chunk %03d/%03d: sims %d-%d",
        Sys.time(), test_label, k_null, j, nrow(chunks), chunks$start[j], chunks$end[j]
      ))
      sim_res <- sim_ssr_nbinom(
        n_sims = chunks$n[j],
        enroll_rate = make_enroll_rate(accrual_scenario_plan),
        fail_rate = make_fail_rate(lambda1_plan, 1.0, k_null),
        dropout_rate = dropout_rate_sim,
        max_followup = max_followup,
        design = gs_plan,
        n_max = n_max,
        min_if_futility = min_if_futility,
        max_enrollment_frac_for_adapt = max_enrollment_frac_for_ia2,
        min_months_to_close_for_adapt = min_months_to_close_for_adapt,
        analysis_lag_months = analysis_lag_months,
        event_gap = event_gap_val,
        ignore_futility = TRUE,
        test_type = test_type_x,
        metadata = list(k_true = k_null, analysis_test = test_type_x),
        seed = 2000000L + i * 1000L + j
      )
      k_chunks[[j]] <- as.data.table(sim_res$trial_results)
      k_chunks[[j]][, sim := sim + chunks$start[j] - 1L]
      saveRDS(k_chunks[[j]], chunk_path)
    }
    all_null[[i]] <- rbindlist(k_chunks, use.names = TRUE, fill = TRUE)
  }

  null_dt <- rbindlist(all_null, use.names = TRUE, fill = TRUE)
  sm <- as.data.table(
    summarize_ssr_sim(null_dt, by = c("k_true", "strategy"))$trial_summary
  )
  sm[, `:=`(
    type1_error = rejection_rate,
    adapted_rate = pct_adapted / 100,
    analysis_test = test_type_x
  )]
  sm <- sm[, .(
    k_true, analysis_test, strategy, n_sims, type1_error,
    cross_ia1, cross_ia2, cross_final, mean_n, adapted_rate
  )]
  rt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  saveRDS(
    list(
      summary = as.data.frame(sm),
      runtime_seconds = rt,
      null_nonbinding_n = null_n,
      alpha_design = alpha_plan,
      test_type = test_type_x,
      generated_at = as.character(Sys.time())
    ),
    cache_path
  )
  list(summary = sm, runtime_seconds = rt)
}

message("=== SSR score production generation ===")
message("production: ", production)
message("workers: ", workers)
message("chunk_size: ", chunk_size)
message("cache_dir: ", cache_dir)
message("output_path: ", output_path)

old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)
future::plan(future::multisession, workers = workers)

sim_start <- Sys.time()
scenario_results <- vector("list", nrow(scenarios))
for (i in seq_len(nrow(scenarios))) {
  scenario_results[[i]] <- run_scenario(i)
}
all_results <- rbindlist(scenario_results, use.names = TRUE, fill = TRUE)
sim_runtime_seconds <- as.numeric(difftime(Sys.time(), sim_start, units = "secs"))

required_cols <- c(
  "ia2_adapt_cut_time",
  "ia2_enroll_frac", "ia2_months_to_close_pred",
  "ia2_adapt_allowed", "ia2_adapt_applied"
)
missing_cols <- setdiff(required_cols, names(all_results))
if (length(missing_cols) > 0) {
  for (nm in missing_cols) all_results[[nm]] <- NA
}

null_n <- if (production) n_sims_production_type1 else n_sims_initial
wald_type1 <- run_null_type1("Wald", "wald", null_n)
score_type1 <- run_null_type1("Score", "score", null_n)
null_nonbinding_by_test <- list(
  Wald = as.data.frame(wald_type1$summary),
  Score = as.data.frame(score_type1$summary)
)
null_nonbinding_runtime_by_test <- list(
  Wald = wald_type1$runtime_seconds,
  Score = score_type1$runtime_seconds
)

plot_cols <- c(
  "lambda1_true", "k_true", "accrual_true", "rr_true", "analysis_test", "strategy",
  "reject", "futility", "reject_stage", "futility_stage",
  "n_adapted", "adapted",
  "participants_with_events_stop", "events_observed_stop",
  "if_ia1", "if_ia2", "if_final",
  "ia1_time", "ia2_time", "final_time",
  "ia1_fallback", "ia2_fallback",
  "participants_with_events_ia1", "participants_with_events_ia2",
  "participants_with_events_final",
  "events_observed_ia1", "events_observed_ia2", "events_observed_final",
  "adapt_cut_time", "adapt_enroll_frac", "adapt_months_to_close_pred",
  "adapt_allowed", "adapt_applied",
  "ia2_adapt_cut_time",
  "ia2_enroll_frac", "ia2_months_to_close_pred",
  "ia2_adapt_allowed", "ia2_adapt_applied"
)
missing_plot_cols <- setdiff(plot_cols, names(all_results))
if (length(missing_plot_cols) > 0) {
  for (nm in missing_plot_cols) all_results[[nm]] <- NA
}

saveRDS(
  list(
    plot_data = as.data.frame(all_results[, ..plot_cols]),
    sim_runtime_seconds = sim_runtime_seconds,
    workers = workers,
    null_nonbinding_summary = null_nonbinding_by_test[[analysis_test_label]],
    null_nonbinding_by_test = null_nonbinding_by_test,
    null_nonbinding_n = null_n,
    null_nonbinding_runtime_seconds =
      null_nonbinding_runtime_by_test[[analysis_test_label]],
    null_nonbinding_runtime_by_test = null_nonbinding_runtime_by_test,
    generated_at = as.character(Sys.time()),
    settings = list(
      analysis_test_type = analysis_test_type,
      n_sims_initial = n_sims_initial,
      n_sims_production_power = n_sims_production_power,
      n_sims_production_type1 = n_sims_production_type1,
      n_sims_production_rr_gt1 = n_sims_production_rr_gt1,
      use_production = production,
      use_production_type1 = production,
      force_main = force_main,
      force_type1 = force_type1,
      cached_scenarios_at_start = cached_scenarios_at_start,
      cache_dir = cache_dir,
      design_note = design_note
    )
  ),
  output_path
)

message("Saved: ", output_path)
message(sprintf("Main simulation runtime: %.2f minutes", sim_runtime_seconds / 60))
message("Type I summaries:")
print(rbindlist(list(
  wald_type1$summary[, test := "Wald"],
  score_type1$summary[, test := "Score"]
), use.names = TRUE, fill = TRUE))
