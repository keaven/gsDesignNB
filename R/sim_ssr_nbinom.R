#' Simulate adaptive group sequential trials with sample size re-estimation
#'
#' Simulates recurrent-event group sequential trials with information-based
#' interim analyses and optional sample size re-estimation (SSR). Interim timing
#' follows the blinded-information targeting used in the SSR study vignette:
#' [get_cut_date()] searches for the earliest cut date where the blinded
#' information reaches the planned information fraction, with an unblinded
#' fallback if the blinded fit is unavailable.
#'
#' The function compares one or more strategies:
#' \describe{
#'   \item{"No adaptation"}{The trial keeps the planned sample size.}
#'   \item{"Blinded SSR"}{A blinded nuisance-parameter update is applied at
#'   `adapt_analysis` using [blinded_ssr()].}
#'   \item{"Unblinded SSR"}{An unblinded nuisance-parameter update is applied at
#'   `adapt_analysis` using [unblinded_ssr()].}
#' }
#'
#' @param n_sims Number of simulated trials.
#' @param enroll_rate Enrollment-rate data frame passed to [nb_sim()].
#' @param fail_rate Failure-rate data frame passed to [nb_sim()]. This defines
#'   the data-generating truth.
#' @param dropout_rate Optional dropout-rate data frame passed to [nb_sim()].
#' @param max_followup Maximum follow-up per patient in the simulated trial.
#' @param design A planning object of class `gsNB`, typically returned by
#'   [gsNBCalendar()] (optionally after [toInteger()]).
#' @param n_max Maximum total enrollment allowed after SSR. Defaults to 150% of
#'   the planned final enrollment, rounded up and constrained to be at least the
#'   planned sample size.
#' @param strategies Character vector of strategies to simulate. Must be chosen
#'   from `"No adaptation"`, `"Blinded SSR"`, and `"Unblinded SSR"`.
#' @param adapt_analysis Interim analysis index where SSR may be applied.
#'   Defaults to the last interim analysis (`design$k - 1`).
#' @param min_if_futility Minimum observed information fraction required before
#'   allowing a futility stop. Default is `0`.
#' @param max_enrollment_frac_for_adapt Maximum fraction of the planned
#'   enrollment already accrued at the adaptation cut for SSR to be allowed.
#'   Default is `1`.
#' @param min_months_to_close_for_adapt Minimum predicted months remaining to
#'   planned enrollment close required to allow SSR. Default is `0`.
#' @param analysis_lag_months Additional months of enrollment counted after a
#'   futility stop to approximate operational lag. Default is `0`.
#' @param event_gap Optional event-gap duration. If `NULL`, inherits
#'   `design$nb_design$inputs$event_gap` when available; otherwise defaults to
#'   `0`.
#' @param bound_info Information measure used for efficacy/futility bounds.
#'   Choices are `"unblinded_ml"` (default), `"blinded_ml"`,
#'   `"unblinded_mom"`, and `"blinded_mom"`.
#' @param first_min_time Minimum calendar time allowed for the first
#'   information-based interim search. Default is `1`.
#' @param min_analysis_gap Minimum gap between successive information-based
#'   interim searches. Default is `0.5`.
#' @param ignore_futility Logical. If `TRUE`, lower-bound crossings do not stop
#'   the trial.
#' @param metadata Optional named list or one-row data frame of scenario labels
#'   to repeat across the returned rows.
#' @param seed Random-seed control passed to `future.apply::future_lapply()`.
#'   Follows the same conventions as [sim_gs_nbinom()].
#'
#' @return An object of class `sim_ssr_nbinom` with components:
#'   \describe{
#'     \item{trial_results}{One row per simulation and strategy, containing
#'     trial-level outcomes and stage-specific summaries in wide format.}
#'     \item{analysis_results}{One row per simulation, strategy, and analysis in
#'     long format.}
#'     \item{settings}{A list of key design and simulation settings.}
#'   }
#'
#' @details
#' The returned `trial_results` include stage-specific columns such as
#' `z_ia1`, `if_ia1`, `ia1_time`, `participants_with_events_ia1`,
#' `events_observed_ia1`, and analogous columns for later analyses. For the
#' actual stopping stage, `participants_with_events_stop` and
#' `events_observed_stop` summarize the observed burden carried into the final
#' decision for each simulated trial.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' enroll_rate <- data.frame(rate = 12, duration = 6)
#' fail_rate <- data.frame(
#'   treatment = c("Control", "Experimental"),
#'   rate = c(0.6, 0.42),
#'   dispersion = 0.4
#' )
#' dropout_rate <- data.frame(
#'   treatment = c("Control", "Experimental"),
#'   rate = c(0.05, 0.05),
#'   duration = c(12, 12)
#' )
#' fixed_design <- sample_size_nbinom(
#'   lambda1 = 0.6,
#'   lambda2 = 0.42,
#'   dispersion = 0.4,
#'   power = 0.8,
#'   alpha = 0.025,
#'   accrual_rate = 12,
#'   accrual_duration = 6,
#'   trial_duration = 12,
#'   max_followup = 6
#' )
#' gs_design <- gsNBCalendar(
#'   fixed_design,
#'   k = 3,
#'   test.type = 4,
#'   alpha = 0.025,
#'   sfu = sfHSD,
#'   sfupar = -2,
#'   sfl = sfHSD,
#'   sflpar = 1,
#'   analysis_times = c(4, 8, 12)
#' )
#' sim_res <- sim_ssr_nbinom(
#'   n_sims = 2,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   dropout_rate = dropout_rate,
#'   max_followup = 6,
#'   design = gs_design,
#'   seed = 123
#' )
#' names(sim_res)
#' head(sim_res$trial_results)
sim_ssr_nbinom <- function(
  n_sims,
  enroll_rate,
  fail_rate,
  dropout_rate = NULL,
  max_followup,
  design,
  n_max = NULL,
  strategies = c("No adaptation", "Blinded SSR", "Unblinded SSR"),
  adapt_analysis = NULL,
  min_if_futility = 0,
  max_enrollment_frac_for_adapt = 1,
  min_months_to_close_for_adapt = 0,
  analysis_lag_months = 0,
  event_gap = NULL,
  bound_info = c("unblinded_ml", "blinded_ml", "unblinded_mom", "blinded_mom"),
  first_min_time = 1,
  min_analysis_gap = 0.5,
  ignore_futility = FALSE,
  metadata = NULL,
  seed = TRUE
) {
  bound_info <- match.arg(bound_info)

  if (!inherits(design, "gsNB")) {
    stop("design must be a 'gsNB' object created by gsNBCalendar().", call. = FALSE)
  }

  strategy_levels <- c("No adaptation", "Blinded SSR", "Unblinded SSR")
  if (!all(strategies %in% strategy_levels)) {
    stop(
      "strategies must be chosen from: ",
      paste(strategy_levels, collapse = ", "),
      call. = FALSE
    )
  }

  strategies <- unique(strategies)
  n_analyses <- design$k
  if (n_analyses < 2) {
    stop("design must include at least one interim and one final analysis.", call. = FALSE)
  }

  if (is.null(adapt_analysis)) {
    adapt_analysis <- n_analyses - 1L
  }
  if (!is.numeric(adapt_analysis) || length(adapt_analysis) != 1L ||
      adapt_analysis < 1 || adapt_analysis >= n_analyses) {
    stop("adapt_analysis must identify an interim analysis.", call. = FALSE)
  }

  inputs <- design$nb_design$inputs
  lambda1_plan <- inputs$lambda1
  lambda2_plan <- inputs$lambda2
  ratio_plan <- inputs$ratio
  power_plan <- 1 - design$beta
  alpha_plan <- design$alpha
  accrual_rate_plan <- inputs$accrual_rate
  accrual_dur_plan <- inputs$accrual_duration
  trial_dur_plan <- inputs$trial_duration
  dropout_rate_plan <- inputs$dropout_rate

  if (is.null(event_gap)) {
    if (!is.null(inputs$event_gap) && !is.na(inputs$event_gap)) {
      event_gap <- inputs$event_gap
    } else {
      event_gap <- 0
    }
  }

  planned_timing <- design$timing
  target_info <- design$n.I[design$k]
  n_planned <- ceiling(utils::tail(design$n_total, 1))
  fixed_n <- design$nb_design$n_total
  gs_inflation <- n_planned / fixed_n

  if (is.null(n_max)) {
    n_max <- ceiling(1.5 * n_planned)
  }
  n_max <- max(as.integer(ceiling(n_max)), as.integer(n_planned))

  metadata_df <- .ssr_metadata_df(metadata)
  stage_labels <- c(paste0("IA", seq_len(n_analyses - 1L)), "Final")
  stage_suffixes <- c(paste0("ia", seq_len(n_analyses - 1L)), "final")
  adapt_label <- stage_labels[adapt_analysis]
  adapt_suffix <- stage_suffixes[adapt_analysis]

  timing_method_at_t <- function(pool, t) {
    cut_dt <- cut_data_by_date(pool, cut_date = t, event_gap = event_gap)
    if (sum(cut_dt$events, na.rm = TRUE) < 2 || length(unique(cut_dt$id)) < 4) {
      return(list(info = 0, method = "insufficient"))
    }

    blind_res <- tryCatch(
      calculate_blinded_info(
        cut_dt,
        ratio = ratio_plan,
        lambda1_planning = lambda1_plan,
        lambda2_planning = lambda2_plan,
        event_gap = event_gap
      ),
      error = function(e) NULL
    )

    if (!is.null(blind_res) &&
        is.finite(blind_res$blinded_info) &&
        !is.na(blind_res$blinded_info) &&
        blind_res$blinded_info > 0) {
      return(list(info = blind_res$blinded_info, method = "blinded"))
    }

    mutze_res <- tryCatch(mutze_test(cut_dt), error = function(e) NULL)
    if (!is.null(mutze_res) &&
        is.finite(mutze_res$se) &&
        !is.na(mutze_res$se) &&
        mutze_res$se > 0) {
      return(list(info = 1 / mutze_res$se^2, method = "unblinded"))
    }

    list(info = 100, method = "fallback100")
  }

  run_one_sim <- function(sim_id) {
    sim_data <- nb_sim(
      enroll_rate = enroll_rate,
      fail_rate = fail_rate,
      dropout_rate = dropout_rate,
      max_followup = max_followup,
      n = n_max,
      event_gap = event_gap
    )

    id_enroll <- unique(sim_data[, c("id", "enroll_time")])
    id_enroll <- id_enroll[order(id_enroll$enroll_time, id_enroll$id), , drop = FALSE]
    all_ids <- id_enroll$id
    max_calendar <- max(sim_data$calendar_time, na.rm = TRUE) + 1

    planned_ids <- all_ids[seq_len(min(n_planned, length(all_ids)))]
    planned_pool <- sim_data[sim_data$id %in% planned_ids, , drop = FALSE]
    class(planned_pool) <- class(sim_data)

    ia_times <- rep(NA_real_, n_analyses)
    timing_methods <- rep(NA_character_, n_analyses)

    if (n_analyses > 1) {
      for (a in seq_len(n_analyses - 1L)) {
        min_t <- if (a == 1L) first_min_time else ia_times[a - 1L] + min_analysis_gap
        target_info_a <- planned_timing[a] * target_info
        ia_times[a] <- tryCatch(
          get_cut_date(
            data = planned_pool,
            target_info = target_info_a,
            event_gap = event_gap,
            lambda1 = lambda1_plan,
            lambda2 = lambda2_plan,
            ratio = ratio_plan,
            min_date = min_t,
            max_date = max_calendar
          ),
          error = function(e) max_calendar
        )
        timing_methods[a] <- timing_method_at_t(planned_pool, ia_times[a])$method
      }
    }

    strategy_results <- vector("list", length(strategies))
    analysis_rows <- list()
    names(strategy_results) <- strategies

    for (strat in strategies) {
      current_target_n <- n_planned
      stopped <- FALSE
      reject <- FALSE
      futility_stop <- FALSE
      reject_stage <- "No reject"
      futility_stage <- "No futility"
      adapt_cut_time <- NA_real_
      adapt_enroll_frac <- NA_real_
      adapt_months_to_close_pred <- NA_real_
      adapt_allowed <- FALSE
      adapt_applied <- FALSE

      z_values <- rep(NA_real_, n_analyses)
      info_values <- rep(NA_real_, n_analyses)
      if_values <- rep(NA_real_, n_analyses)
      cal_times <- rep(NA_real_, n_analyses)
      n_enrolled_values <- rep(NA_real_, n_analyses)
      participants_with_events_values <- rep(NA_real_, n_analyses)
      events_observed_values <- rep(NA_real_, n_analyses)
      stage_fallback <- rep(FALSE, n_analyses)

      for (a in seq_len(n_analyses)) {
        if (stopped) {
          break
        }

        if (a < n_analyses) {
          cut_time <- ia_times[a]
          interim_n <- min(n_planned, length(all_ids))
          analysis_ids <- all_ids[seq_len(interim_n)]
          enrolled_by_t <- id_enroll$id[
            id_enroll$id %in% analysis_ids & id_enroll$enroll_time <= cut_time
          ]
          analysis_pool <- sim_data[sim_data$id %in% enrolled_by_t, , drop = FALSE]
          class(analysis_pool) <- class(sim_data)

          if (a == adapt_analysis) {
            enrolled_n_stage <- length(enrolled_by_t)
            accrual_hat_stage <- if (cut_time > 0) enrolled_n_stage / cut_time else NA_real_
            predicted_enroll_close <- if (!is.na(accrual_hat_stage) && accrual_hat_stage > 0) {
              interim_n / accrual_hat_stage
            } else {
              NA_real_
            }

            adapt_cut_time <- if (!is.na(predicted_enroll_close)) {
              min(cut_time, predicted_enroll_close - min_months_to_close_for_adapt)
            } else {
              cut_time
            }
            if (adapt_analysis > 1) {
              adapt_cut_time <- max(adapt_cut_time, ia_times[adapt_analysis - 1L] + min_analysis_gap)
            }

            enrolled_n_adapt <- sum(
              id_enroll$id %in% analysis_ids &
                id_enroll$enroll_time <= adapt_cut_time
            )
            adapt_enroll_frac <- if (interim_n > 0) enrolled_n_adapt / interim_n else NA_real_
            accrual_hat <- if (adapt_cut_time > 0) {
              enrolled_n_adapt / adapt_cut_time
            } else {
              NA_real_
            }
            adapt_months_to_close_pred <- if (!is.na(accrual_hat) && accrual_hat > 0) {
              max(0, (interim_n - enrolled_n_adapt) / accrual_hat)
            } else {
              NA_real_
            }
            adapt_allowed <- !is.na(adapt_enroll_frac) &&
              is.finite(adapt_cut_time) &&
              adapt_cut_time <= cut_time &&
              adapt_enroll_frac <= max_enrollment_frac_for_adapt
          }
        } else {
          analysis_ids <- all_ids[seq_len(min(current_target_n, length(all_ids)))]
          analysis_pool <- sim_data[sim_data$id %in% analysis_ids, , drop = FALSE]
          class(analysis_pool) <- class(sim_data)
          last_enroll <- max(id_enroll$enroll_time[id_enroll$id %in% analysis_ids])
          cut_time <- max(
            max(ia_times, na.rm = TRUE),
            last_enroll + max_followup,
            max(analysis_pool$calendar_time, na.rm = TRUE)
          )
          timing_methods[a] <- "final_calendar"
        }

        cal_times[a] <- cut_time
        cut_data <- cut_data_by_date(analysis_pool, cut_date = cut_time, event_gap = event_gap)

        n_enrolled_values[a] <- nrow(cut_data)
        participants_with_events_values[a] <- sum(cut_data$events > 0, na.rm = TRUE)
        events_observed_values[a] <- sum(cut_data$events, na.rm = TRUE)

        metrics <- .ssr_extract_info_metrics(
          cut_data = cut_data,
          ratio_plan = ratio_plan,
          lambda1_plan = lambda1_plan,
          lambda2_plan = lambda2_plan,
          event_gap = event_gap
        )

        selected_info <- .ssr_select_info(metrics, bound_info)
        if (is.finite(selected_info) && !is.na(selected_info) && selected_info > 0) {
          info_values[a] <- selected_info
          if_values[a] <- selected_info / target_info
        }

        if (is.finite(metrics$z_value)) {
          z_values[a] <- metrics$z_value
        }

        stage_fallback[a] <- !is.na(timing_methods[a]) && timing_methods[a] != "blinded" &&
          a < n_analyses

        upper_bound <- NA_real_
        lower_bound <- NA_real_
        cross_upper <- FALSE
        cross_lower <- FALSE

        if (!is.na(z_values[a]) && !is.na(if_values[a])) {
          bounds <- .ssr_dynamic_bounds(
            analysis_index = a,
            n_analyses = n_analyses,
            observed_if = if_values,
            planned_timing = planned_timing,
            target_info = target_info,
            design = design
          )
          upper_bound <- bounds$upper_bound
          lower_bound <- bounds$lower_bound

          if (is.finite(upper_bound) && z_values[a] > upper_bound) {
            reject <- TRUE
            reject_stage <- stage_labels[a]
            cross_upper <- TRUE
            stopped <- TRUE
          } else if (a < n_analyses &&
                     !ignore_futility &&
                     if_values[a] >= min_if_futility &&
                     is.finite(lower_bound) &&
                     z_values[a] < lower_bound) {
            futility_stop <- TRUE
            futility_stage <- stage_labels[a]
            cross_lower <- TRUE
            stopped <- TRUE
          }
        }

        if (cross_lower && a < n_analyses) {
          n_futility_report <- sum(
            id_enroll$id %in% analysis_ids &
              id_enroll$enroll_time <= (cut_time + analysis_lag_months)
          )
          n_futility_report <- max(n_futility_report, length(unique(cut_data$id)))
          n_futility_report <- min(n_futility_report, length(analysis_ids))
          current_target_n <- n_futility_report
        }

        analysis_rows[[length(analysis_rows) + 1L]] <- data.frame(
          sim = sim_id,
          strategy = strat,
          analysis = a,
          analysis_stage = stage_labels[a],
          analysis_time = cut_time,
          n_target = current_target_n,
          n_enrolled = n_enrolled_values[a],
          participants_with_events = participants_with_events_values[a],
          events_observed = events_observed_values[a],
          z_value = z_values[a],
          info_value = info_values[a],
          info_fraction = if_values[a],
          info_unblinded_ml = metrics$info_unblinded_ml,
          info_blinded_ml = metrics$info_blinded_ml,
          info_unblinded_mom = metrics$info_unblinded_mom,
          info_blinded_mom = metrics$info_blinded_mom,
          timing_info_method = timing_methods[a],
          timing_fallback = stage_fallback[a],
          upper_bound = upper_bound,
          lower_bound = lower_bound,
          cross_upper = cross_upper,
          cross_lower = cross_lower,
          reject = reject,
          futility = futility_stop,
          adapted = current_target_n > n_planned,
          stringsAsFactors = FALSE
        )

        if (a == adapt_analysis &&
            !stopped &&
            strat != "No adaptation" &&
            adapt_allowed) {
          adapt_pool <- sim_data[sim_data$id %in% analysis_ids, , drop = FALSE]
          class(adapt_pool) <- class(sim_data)
          adapt_cut_data <- if (!is.na(adapt_cut_time)) {
            cut_data_by_date(adapt_pool, cut_date = adapt_cut_time, event_gap = event_gap)
          } else {
            cut_data
          }

          ssr_target_n <- NULL
          if (nrow(adapt_cut_data) >= 10) {
            if (strat == "Blinded SSR") {
              blind_ssr <- tryCatch(
                blinded_ssr(
                  data = adapt_cut_data,
                  ratio = ratio_plan,
                  lambda1_planning = lambda1_plan,
                  lambda2_planning = lambda2_plan,
                  power = power_plan,
                  alpha = alpha_plan,
                  accrual_rate = accrual_rate_plan,
                  accrual_duration = accrual_dur_plan,
                  trial_duration = trial_dur_plan,
                  dropout_rate = dropout_rate_plan,
                  max_followup = max_followup,
                  event_gap = event_gap
                ),
                error = function(e) NULL
              )
              if (!is.null(blind_ssr) && is.finite(blind_ssr$n_total_blinded)) {
                ssr_target_n <- blind_ssr$n_total_blinded
              }
            } else {
              unblind_ssr <- tryCatch(
                unblinded_ssr(
                  data = adapt_cut_data,
                  ratio = ratio_plan,
                  lambda1_planning = lambda1_plan,
                  lambda2_planning = lambda2_plan,
                  power = power_plan,
                  alpha = alpha_plan,
                  accrual_rate = accrual_rate_plan,
                  accrual_duration = accrual_dur_plan,
                  trial_duration = trial_dur_plan,
                  dropout_rate = dropout_rate_plan,
                  max_followup = max_followup,
                  event_gap = event_gap
                ),
                error = function(e) NULL
              )
              if (!is.null(unblind_ssr) && is.finite(unblind_ssr$n_total_unblinded)) {
                ssr_target_n <- unblind_ssr$n_total_unblinded
              }
            }
          }

          if (!is.null(ssr_target_n)) {
            n_new <- min(
              max(ceiling(ssr_target_n * gs_inflation), n_planned),
              n_max
            )
            adapt_applied <- n_new > n_planned
            current_target_n <- n_new
          }
        }
      }

      observed_analyses <- which(!is.na(cal_times))
      last_analysis <- if (length(observed_analyses) == 0L) NA_integer_ else max(observed_analyses)
      stop_stage <- if (length(last_analysis) == 0 || !is.finite(last_analysis)) {
        NA_character_
      } else if (reject) {
        reject_stage
      } else if (futility_stop) {
        futility_stage
      } else {
        stage_labels[last_analysis]
      }

      participants_with_events_stop <- if (!is.na(stop_stage) && length(last_analysis) == 1L) {
        participants_with_events_values[last_analysis]
      } else {
        NA_real_
      }
      events_observed_stop <- if (!is.na(stop_stage) && length(last_analysis) == 1L) {
        events_observed_values[last_analysis]
      } else {
        NA_real_
      }

      trial_row <- list(
        sim = sim_id,
        strategy = strat,
        reject = reject,
        futility = futility_stop,
        reject_stage = reject_stage,
        futility_stage = futility_stage,
        stop_stage = stop_stage,
        n_adapted = current_target_n,
        adapted = current_target_n > n_planned,
        adapt_stage = adapt_label,
        adapt_cut_time = adapt_cut_time,
        adapt_enroll_frac = adapt_enroll_frac,
        adapt_months_to_close_pred = adapt_months_to_close_pred,
        adapt_allowed = adapt_allowed,
        adapt_applied = adapt_applied,
        participants_with_events_stop = participants_with_events_stop,
        events_observed_stop = events_observed_stop
      )

      if (identical(adapt_suffix, "ia2")) {
        trial_row$ia2_adapt_cut_time <- adapt_cut_time
        trial_row$ia2_enroll_frac <- adapt_enroll_frac
        trial_row$ia2_months_to_close_pred <- adapt_months_to_close_pred
        trial_row$ia2_adapt_allowed <- adapt_allowed
        trial_row$ia2_adapt_applied <- adapt_applied
      }

      for (idx in seq_along(stage_suffixes)) {
        suffix <- stage_suffixes[idx]
        trial_row[[paste0("z_", suffix)]] <- z_values[idx]
        trial_row[[paste0("if_", suffix)]] <- if_values[idx]
        trial_row[[paste0("info_", suffix)]] <- info_values[idx]
        trial_row[[paste0(suffix, "_time")]] <- cal_times[idx]
        trial_row[[paste0("n_enrolled_", suffix)]] <- n_enrolled_values[idx]
        trial_row[[paste0("participants_with_events_", suffix)]] <-
          participants_with_events_values[idx]
        trial_row[[paste0("events_observed_", suffix)]] <- events_observed_values[idx]
        trial_row[[paste0(suffix, "_info_method")]] <- timing_methods[idx]
        trial_row[[paste0(suffix, "_fallback")]] <- stage_fallback[idx]
      }

      strategy_results[[strat]] <- data.frame(trial_row, stringsAsFactors = FALSE)
    }

    trial_results <- do.call(rbind, strategy_results)
    analysis_results <- do.call(rbind, analysis_rows)

    if (!is.null(metadata_df)) {
      trial_results <- cbind(metadata_df[rep(1, nrow(trial_results)), , drop = FALSE], trial_results)
      analysis_results <- cbind(metadata_df[rep(1, nrow(analysis_results)), , drop = FALSE], analysis_results)
    }

    list(
      trial_results = trial_results,
      analysis_results = analysis_results
    )
  }

  has_future_apply <- requireNamespace("future.apply", quietly = TRUE)
  if (has_future_apply) {
    results_list <- future.apply::future_lapply(
      seq_len(n_sims),
      run_one_sim,
      future.seed = seed
    )
  } else {
    if (is.numeric(seed)) {
      set.seed(seed)
    }
    results_list <- lapply(seq_len(n_sims), run_one_sim)
  }

  trial_results <- do.call(
    rbind,
    lapply(results_list, function(x) x$trial_results)
  )
  analysis_results <- do.call(
    rbind,
    lapply(results_list, function(x) x$analysis_results)
  )

  structure(
    list(
      trial_results = trial_results,
      analysis_results = analysis_results,
      settings = list(
        n_sims = n_sims,
        n_planned = n_planned,
        n_max = n_max,
        target_info = target_info,
        bound_info = bound_info,
        strategies = strategies,
        adapt_analysis = adapt_analysis,
        adapt_stage = adapt_label,
        planned_timing = planned_timing,
        event_gap = event_gap,
        min_if_futility = min_if_futility,
        max_enrollment_frac_for_adapt = max_enrollment_frac_for_adapt,
        min_months_to_close_for_adapt = min_months_to_close_for_adapt,
        analysis_lag_months = analysis_lag_months,
        gs_inflation = gs_inflation
      )
    ),
    class = "sim_ssr_nbinom"
  )
}

#' Summarize adaptive SSR simulation results
#'
#' Produces trial-level and analysis-level summaries from the output of
#' `sim_ssr_nbinom()`. The summary includes expected sample size, stopping
#' probabilities, expected participants with events, and expected events
#' observed.
#'
#' @param x A `sim_ssr_nbinom` object or a data frame containing the
#'   `trial_results` component returned by `sim_ssr_nbinom()`.
#' @param by Character vector of grouping columns for the trial-level summary.
#'   Defaults to `"strategy"`.
#'
#' @return A list with:
#'   \describe{
#'     \item{trial_summary}{Trial-level grouped summary.}
#'     \item{analysis_summary}{Analysis-level grouped summary.}
#'   }
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' enroll_rate <- data.frame(rate = 10, duration = 4)
#' fail_rate <- data.frame(
#'   treatment = c("Control", "Experimental"),
#'   rate = c(0.5, 0.35),
#'   dispersion = 0.3
#' )
#' fixed_design <- sample_size_nbinom(
#'   lambda1 = 0.5,
#'   lambda2 = 0.35,
#'   dispersion = 0.3,
#'   power = 0.8,
#'   alpha = 0.025,
#'   accrual_rate = 10,
#'   accrual_duration = 4,
#'   trial_duration = 8,
#'   max_followup = 4
#' )
#' gs_design <- gsNBCalendar(
#'   fixed_design,
#'   k = 3,
#'   test.type = 4,
#'   alpha = 0.025,
#'   analysis_times = c(3, 5, 8)
#' )
#' sim_res <- sim_ssr_nbinom(
#'   n_sims = 2,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   max_followup = 4,
#'   design = gs_design,
#'   strategies = "No adaptation",
#'   seed = 321
#' )
#' summarize_ssr_sim(sim_res)$trial_summary
summarize_ssr_sim <- function(x, by = "strategy") {
  analysis_time <- participants_with_events <- events_observed <- NULL
  info_value <- info_fraction <- timing_fallback <- NULL

  if (inherits(x, "sim_ssr_nbinom")) {
    trial_results <- x$trial_results
    analysis_results <- x$analysis_results
  } else {
    trial_results <- x
    analysis_results <- .ssr_long_from_trial_results(trial_results)
  }

  trial_dt <- data.table::as.data.table(trial_results)
  by_cols <- intersect(by, names(trial_dt))
  if (length(by_cols) == 0L) {
    by_cols <- NULL
  }

  trial_summary_expr <- quote(.(
    n_sims = .N,
    rejection_rate = mean(reject, na.rm = TRUE),
    futility_rate = mean(futility, na.rm = TRUE),
    mean_n = mean(n_adapted, na.rm = TRUE),
    sd_n = stats::sd(n_adapted, na.rm = TRUE),
    pct_adapted = 100 * mean(adapted, na.rm = TRUE),
    expected_participants_with_events = if ("participants_with_events_stop" %in% names(trial_dt)) {
      mean(participants_with_events_stop, na.rm = TRUE)
    } else {
      NA_real_
    },
    expected_events_observed = if ("events_observed_stop" %in% names(trial_dt)) {
      mean(events_observed_stop, na.rm = TRUE)
    } else {
      NA_real_
    }
  ))
  trial_summary <- if (is.null(by_cols)) {
    trial_dt[, eval(trial_summary_expr)]
  } else {
    trial_dt[, eval(trial_summary_expr), by = by_cols]
  }

  stage_suffixes <- .ssr_stage_suffixes_from_trial_results(trial_dt)
  for (suffix in stage_suffixes) {
    stage_expr <- quote(.(
      pct_reach = 100 * mean(!is.na(get(paste0(suffix, "_time"))), na.rm = TRUE),
      mean_if = mean(get(paste0("if_", suffix)), na.rm = TRUE),
      mean_time = mean(get(paste0(suffix, "_time")), na.rm = TRUE),
      mean_participants_with_events = if (paste0("participants_with_events_", suffix) %in% names(trial_dt)) {
        mean(get(paste0("participants_with_events_", suffix)), na.rm = TRUE)
      } else {
        NA_real_
      },
      mean_events_observed = if (paste0("events_observed_", suffix) %in% names(trial_dt)) {
        mean(get(paste0("events_observed_", suffix)), na.rm = TRUE)
      } else {
        NA_real_
      },
      n_fallback = if (paste0(suffix, "_fallback") %in% names(trial_dt)) {
        sum(get(paste0(suffix, "_fallback")), na.rm = TRUE)
      } else {
        NA_real_
      },
      prob_reject = mean(reject_stage == .ssr_stage_label_from_suffix(suffix), na.rm = TRUE),
      prob_futility = mean(futility_stage == .ssr_stage_label_from_suffix(suffix), na.rm = TRUE)
    ))
    stage_summary <- if (is.null(by_cols)) {
      trial_dt[, eval(stage_expr)]
    } else {
      trial_dt[, eval(stage_expr), by = by_cols]
    }

    data.table::setnames(
      stage_summary,
      old = c(
        "pct_reach", "mean_if", "mean_time", "mean_participants_with_events",
        "mean_events_observed", "n_fallback", "prob_reject", "prob_futility"
      ),
      new = c(
        paste0("pct_reach_", suffix),
        paste0("mean_if_", suffix),
        paste0("mean_", suffix, "_time"),
        paste0("mean_participants_with_events_", suffix),
        paste0("mean_events_observed_", suffix),
        paste0("n_fallback_", suffix),
        paste0("cross_", suffix),
        paste0("pct_futility_", suffix)
      )
    )

    trial_summary <- if (is.null(by_cols)) {
      cbind(trial_summary, stage_summary)
    } else {
      merge(trial_summary, stage_summary, by = by_cols, all = TRUE)
    }
  }

  if ("adapt_cut_time" %in% names(trial_dt)) {
    adapt_expr <- quote(.(
      mean_adapt_cut_time = mean(adapt_cut_time, na.rm = TRUE),
      mean_adapt_enroll_pct = 100 * mean(adapt_enroll_frac, na.rm = TRUE),
      mean_adapt_months_to_close = mean(adapt_months_to_close_pred, na.rm = TRUE),
      pct_adapt_allowed = 100 * mean(adapt_allowed, na.rm = TRUE),
      pct_adapt_applied = 100 * mean(adapt_applied, na.rm = TRUE)
    ))
    adapt_summary <- if (is.null(by_cols)) {
      trial_dt[, eval(adapt_expr)]
    } else {
      trial_dt[, eval(adapt_expr), by = by_cols]
    }

    trial_summary <- if (is.null(by_cols)) {
      cbind(trial_summary, adapt_summary)
    } else {
      merge(trial_summary, adapt_summary, by = by_cols, all = TRUE)
    }
  }

  analysis_dt <- data.table::as.data.table(analysis_results)
  if (nrow(analysis_dt) > 0) {
    analysis_by <- c(by_cols, "analysis", "analysis_stage")
    analysis_summary <- analysis_dt[, .(
      n_rows = .N,
      pct_reach = 100 * mean(!is.na(analysis_time), na.rm = TRUE),
      mean_time = mean(analysis_time, na.rm = TRUE),
      mean_n_enrolled = mean(n_enrolled, na.rm = TRUE),
      mean_participants_with_events = mean(participants_with_events, na.rm = TRUE),
      mean_events_observed = mean(events_observed, na.rm = TRUE),
      mean_info = mean(info_value, na.rm = TRUE),
      mean_if = mean(info_fraction, na.rm = TRUE),
      prob_cross_upper = mean(cross_upper, na.rm = TRUE),
      prob_cross_lower = mean(cross_lower, na.rm = TRUE),
      n_fallback = sum(timing_fallback, na.rm = TRUE)
    ), by = analysis_by]
  } else {
    analysis_summary <- data.frame()
  }

  list(
    trial_summary = as.data.frame(trial_summary),
    analysis_summary = as.data.frame(analysis_summary)
  )
}

.ssr_metadata_df <- function(metadata) {
  if (is.null(metadata)) {
    return(NULL)
  }

  if (is.data.frame(metadata)) {
    if (nrow(metadata) != 1L) {
      stop("metadata data frame must have exactly one row.", call. = FALSE)
    }
    return(metadata)
  }

  if (is.list(metadata) && !is.null(names(metadata))) {
    return(as.data.frame(metadata, stringsAsFactors = FALSE))
  }

  stop("metadata must be a named list or one-row data frame.", call. = FALSE)
}

.ssr_extract_info_metrics <- function(
  cut_data,
  ratio_plan,
  lambda1_plan,
  lambda2_plan,
  event_gap
) {
  res <- list(
    z_value = NA_real_,
    method_used = NA_character_,
    dispersion_model = NA_real_,
    info_unblinded_ml = NA_real_,
    info_blinded_ml = NA_real_,
    info_unblinded_mom = NA_real_,
    info_blinded_mom = NA_real_
  )

  if (nrow(cut_data) < 4 || sum(cut_data$events, na.rm = TRUE) < 2) {
    return(res)
  }

  calc_info <- function(dat, lam1, lam2, dispersion) {
    if (!is.finite(lam1) || !is.finite(lam2) || !is.finite(dispersion) || dispersion < 0) {
      return(NA_real_)
    }

    dat_c <- dat[dat$treatment == "Control", , drop = FALSE]
    dat_e <- dat[dat$treatment == "Experimental", , drop = FALSE]

    mu_c <- lam1 * dat_c$tte
    mu_e <- lam2 * dat_e$tte
    term_c <- sum(mu_c / (1 + dispersion * mu_c))
    term_e <- sum(mu_e / (1 + dispersion * mu_e))

    if (!is.finite(term_c) || !is.finite(term_e) || term_c <= 0 || term_e <= 0) {
      return(NA_real_)
    }

    1 / (1 / term_c + 1 / term_e)
  }

  test_res <- tryCatch(mutze_test(cut_data), error = function(e) NULL)
  if (!is.null(test_res) && is.finite(test_res$se) && !is.na(test_res$se) && test_res$se > 0) {
    res$z_value <- -test_res$z
    res$method_used <- test_res$method
    res$dispersion_model <- test_res$dispersion
    res$info_unblinded_ml <- 1 / test_res$se^2
  }

  blinded_res <- tryCatch(
    calculate_blinded_info(
      cut_data,
      ratio = ratio_plan,
      lambda1_planning = lambda1_plan,
      lambda2_planning = lambda2_plan,
      event_gap = event_gap
    ),
    error = function(e) NULL
  )
  if (!is.null(blinded_res)) {
    res$info_blinded_ml <- blinded_res$blinded_info
  }

  mom_unblinded <- tryCatch(
    estimate_nb_mom(cut_data, group = "treatment"),
    error = function(e) NULL
  )
  if (!is.null(mom_unblinded) && !any(is.na(mom_unblinded$lambda))) {
    lambda1_mom <- mom_unblinded$lambda["Control"]
    lambda2_mom <- mom_unblinded$lambda["Experimental"]
    if (is.na(lambda1_mom)) {
      lambda1_mom <- mom_unblinded$lambda[1]
    }
    if (is.na(lambda2_mom)) {
      lambda2_mom <- mom_unblinded$lambda[2]
    }
    res$info_unblinded_mom <- calc_info(
      cut_data,
      lambda1_mom,
      lambda2_mom,
      mom_unblinded$dispersion
    )
  }

  mom_blinded <- tryCatch(estimate_nb_mom(cut_data), error = function(e) NULL)
  if (!is.null(mom_blinded) &&
      is.finite(mom_blinded$lambda) &&
      is.finite(mom_blinded$dispersion) &&
      mom_blinded$dispersion >= 0) {
    p1 <- 1 / (1 + ratio_plan)
    p2 <- ratio_plan / (1 + ratio_plan)
    rr_plan <- lambda2_plan / lambda1_plan
    lambda1_adj <- mom_blinded$lambda / (p1 + p2 * rr_plan)
    lambda2_adj <- lambda1_adj * rr_plan

    mu1 <- lambda1_adj * cut_data$tte
    mu2 <- lambda2_adj * cut_data$tte
    w1 <- p1 * sum(mu1 / (1 + mom_blinded$dispersion * mu1))
    w2 <- p2 * sum(mu2 / (1 + mom_blinded$dispersion * mu2))
    if (is.finite(w1) && is.finite(w2) && w1 > 0 && w2 > 0) {
      res$info_blinded_mom <- 1 / (1 / w1 + 1 / w2)
    }
  }

  res
}

.ssr_select_info <- function(metrics, key) {
  value <- switch(
    key,
    unblinded_ml = metrics$info_unblinded_ml,
    blinded_ml = metrics$info_blinded_ml,
    unblinded_mom = metrics$info_unblinded_mom,
    blinded_mom = metrics$info_blinded_mom
  )

  if (!is.finite(value) || is.na(value) || value <= 0) {
    metrics$info_unblinded_ml
  } else {
    value
  }
}

.ssr_dynamic_bounds <- function(
  analysis_index,
  n_analyses,
  observed_if,
  planned_timing,
  target_info,
  design
) {
  timing_vec <- c(observed_if[seq_len(analysis_index)], if (analysis_index < n_analyses) 1)
  spend_vec <- c(
    pmin(planned_timing[seq_len(analysis_index)], observed_if[seq_len(analysis_index)]),
    if (analysis_index < n_analyses) 1
  )

  if (analysis_index == n_analyses) {
    timing_vec[analysis_index] <- 1
    spend_vec[analysis_index] <- 1
  }

  for (j in 2:length(timing_vec)) {
    if (!is.na(timing_vec[j]) && !is.na(timing_vec[j - 1L]) &&
        timing_vec[j] <= timing_vec[j - 1L]) {
      timing_vec[j] <- timing_vec[j - 1L] + 0.001
    }
    if (!is.na(spend_vec[j]) && !is.na(spend_vec[j - 1L]) &&
        spend_vec[j] <= spend_vec[j - 1L]) {
      spend_vec[j] <- spend_vec[j - 1L] + 0.001
    }
  }

  if (any(is.na(timing_vec))) {
    return(list(upper_bound = NA_real_, lower_bound = NA_real_))
  }

  gs_args <- list(
    k = length(timing_vec),
    test.type = design$test.type,
    alpha = design$alpha,
    beta = design$beta,
    sfu = design$upper$sf,
    sfupar = design$upper$param,
    sfl = design$lower$sf,
    sflpar = design$lower$param,
    n.fix = target_info,
    timing = timing_vec,
    usTime = spend_vec,
    lsTime = spend_vec
  )

  if (design$test.type %in% c(7, 8) && !is.null(design$harm)) {
    gs_args$sfharm <- design$harm$sf
    gs_args$sfharmparam <- design$harm$param
  }

  gs_updated <- tryCatch(
    do.call(gsDesign::gsDesign, gs_args),
    error = function(e) NULL
  )

  if (is.null(gs_updated)) {
    return(list(upper_bound = NA_real_, lower_bound = NA_real_))
  }

  list(
    upper_bound = gs_updated$upper$bound[analysis_index],
    lower_bound = gs_updated$lower$bound[analysis_index]
  )
}

.ssr_stage_suffixes_from_trial_results <- function(dt) {
  cols <- names(dt)
  unique(sub("^if_", "", cols[grepl("^if_", cols)]))
}

.ssr_stage_label_from_suffix <- function(suffix) {
  if (identical(suffix, "final")) {
    "Final"
  } else {
    paste0("IA", sub("^ia", "", suffix))
  }
}

.ssr_long_from_trial_results <- function(trial_results) {
  dt <- data.table::as.data.table(trial_results)
  suffixes <- .ssr_stage_suffixes_from_trial_results(dt)
  if (length(suffixes) == 0L) {
    return(data.frame())
  }

  rows <- lapply(seq_along(suffixes), function(idx) {
    suffix <- suffixes[idx]
    required <- c(
      paste0(suffix, "_time"),
      paste0("n_enrolled_", suffix),
      paste0("participants_with_events_", suffix),
      paste0("events_observed_", suffix),
      paste0("if_", suffix),
      paste0("info_", suffix),
      paste0("z_", suffix)
    )
    if (!all(required %in% names(dt))) {
      return(NULL)
    }

    stage_specific_cols <- c(
      paste0("z_", suffixes),
      paste0("if_", suffixes),
      paste0("info_", suffixes),
      paste0("n_enrolled_", suffixes),
      paste0("participants_with_events_", suffixes),
      paste0("events_observed_", suffixes),
      paste0(suffixes, "_time"),
      paste0(suffixes, "_fallback"),
      paste0(suffixes, "_info_method")
    )
    base_cols <- setdiff(names(dt), stage_specific_cols)

    data.frame(
      as.data.frame(dt[, base_cols, with = FALSE]),
      analysis = idx,
      analysis_stage = .ssr_stage_label_from_suffix(suffix),
      analysis_time = dt[[paste0(suffix, "_time")]],
      n_enrolled = dt[[paste0("n_enrolled_", suffix)]],
      participants_with_events = dt[[paste0("participants_with_events_", suffix)]],
      events_observed = dt[[paste0("events_observed_", suffix)]],
      z_value = dt[[paste0("z_", suffix)]],
      info_value = dt[[paste0("info_", suffix)]],
      info_fraction = dt[[paste0("if_", suffix)]],
      timing_fallback = if (paste0(suffix, "_fallback") %in% names(dt)) {
        dt[[paste0(suffix, "_fallback")]]
      } else {
        FALSE
      },
      cross_upper = dt$reject_stage == .ssr_stage_label_from_suffix(suffix),
      cross_lower = dt$futility_stage == .ssr_stage_label_from_suffix(suffix),
      stringsAsFactors = FALSE
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(data.frame())
  }

  do.call(rbind, rows)
}
