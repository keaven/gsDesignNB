#' Simulate group sequential clinical trial for negative binomial outcomes
#'
#' Simulates multiple replicates of a group sequential clinical trial with negative
#' binomial outcomes, performing interim analyses at specified calendar times.
#' Supports parallel execution via the \pkg{future} framework for faster simulation
#' with reproducible random number generation.
#'
#' @param n_sims Number of simulations to run.
#' @param enroll_rate Enrollment rates (data frame with `rate` and `duration`).
#' @param fail_rate Failure rates (data frame with `treatment`, `rate`, `dispersion`).
#' @param dropout_rate Dropout rates (data frame with `treatment`, `rate`, `duration`).
#' @param max_followup Maximum follow-up time.
#' @param event_gap Event gap duration. If `NULL`, inherits
#'   `design$inputs$event_gap` when available; otherwise defaults to `0`.
#' @param analysis_times Vector of calendar times for interim and final analyses.
#'   Optional if `cuts` is provided.
#' @param n_target Total sample size to enroll (optional, if not defined by `enroll_rate`).
#' @param design An object of class `gsNB` or `sample_size_nbinom_result`.
#'   Used to extract planning parameters (`lambda1`, `lambda2`, `ratio`) for blinded
#'   information estimation.
#' @param data_cut Function to cut data for analysis. Defaults to [cut_data_by_date()].
#'   The function must accept `sim_data`, `cut_date`, and `event_gap` as arguments.
#' @param cuts A list of cutting criteria for each analysis. Each element of the list
#'   should be a list of arguments for [get_cut_date()] (e.g., `planned_calendar`,
#'   `target_events`, `target_info`). If provided, `analysis_times` is ignored
#'   (or used as a fallback if `planned_calendar` is missing in a cut).
#' @param seed Random seed for reproducible simulations. Controls the
#'   `future.seed` argument of [future.apply::future_lapply()]:
#'   - `TRUE` (default): Automatically generates parallel-safe L'Ecuyer-CMRG
#'     random number streams. Results are reproducible when preceded by
#'     [set.seed()] regardless of the number of workers.
#'   - An integer: Used as the seed for L'Ecuyer-CMRG streams directly
#'     (equivalent to calling [set.seed()] with this value before the run).
#'   - `FALSE` or `NULL`: No special RNG handling (not recommended; results
#'     may not be reproducible in parallel).
#'
#'   When \pkg{future.apply} is not installed, `seed` is used with
#'   [set.seed()] for sequential execution. See **Details** for parallel usage.
#'
#' @details
#' ## Parallel execution
#'
#' This function uses [future.apply::future_lapply()] to distribute simulation
#' replicates across workers. By default, simulations run sequentially
#' (equivalent to [lapply()]). To enable parallel execution, set a
#' \pkg{future} plan before calling this function:
#'
#' ```
#' library(future)
#' plan(multisession, workers = 4)   # use 4 parallel workers
#' results <- sim_gs_nbinom(...)
#' plan(sequential)                  # restore default
#' ```
#'
#' ## Reproducibility
#'
#' The default `seed = TRUE` ensures that results are fully reproducible
#' regardless of the \pkg{future} plan (sequential or parallel) and regardless
#' of the number of workers. This is achieved via the L'Ecuyer-CMRG algorithm
#' which generates statistically independent random number streams for each
#' simulation replicate. To obtain the same results across runs:
#'
#' ```
#' set.seed(42)
#' res1 <- sim_gs_nbinom(n_sims = 100, ..., seed = TRUE)
#'
#' set.seed(42)
#' res2 <- sim_gs_nbinom(n_sims = 100, ..., seed = TRUE)
#'
#' identical(res1, res2)  # TRUE, even with different plan()
#' ```
#'
#' @return A data frame containing simulation results for each analysis of each trial.
#'   Columns include:
#'   \describe{
#'     \item{sim}{Simulation ID}
#'     \item{analysis}{Analysis index}
#'     \item{analysis_time}{Calendar time of analysis}
#'     \item{n_enrolled}{Number of subjects enrolled}
#'     \item{n_ctrl}{Number of subjects in control group}
#'     \item{n_exp}{Number of subjects in experimental group}
#'     \item{events_total}{Total events observed}
#'     \item{events_ctrl}{Events in control group}
#'     \item{events_exp}{Events in experimental group}
#'     \item{exposure_at_risk_ctrl}{Exposure at risk in control group (adjusted for event gaps)}
#'     \item{exposure_at_risk_exp}{Exposure at risk in experimental group (adjusted for event gaps)}
#'     \item{exposure_total_ctrl}{Total exposure in control group (calendar follow-up)}
#'     \item{exposure_total_exp}{Total exposure in experimental group (calendar follow-up)}
#'     \item{z_stat}{Z-statistic from the Wald test (positive favors experimental if rate ratio < 1)}
#'     \item{estimate}{Estimated log rate ratio from the model}
#'     \item{se}{Standard error of the estimate}
#'     \item{method_used}{Method used for inference ("nb" or "poisson")}
#'     \item{dispersion}{Estimated dispersion parameter from the model}
#'     \item{blinded_info}{Estimated blinded statistical information (ML)}
#'     \item{unblinded_info}{Observed unblinded statistical information (ML)}
#'     \item{info_unblinded_ml}{Observed unblinded statistical information (ML)}
#'     \item{info_blinded_ml}{Estimated blinded statistical information (ML)}
#'     \item{info_unblinded_mom}{Observed unblinded statistical information (Method of Moments)}
#'     \item{info_blinded_mom}{Estimated blinded statistical information (Method of Moments)}
#'   }
#'
#' @importFrom data.table as.data.table
#'
#' @export
#'
#' @examples
#' # Basic sequential usage with reproducible seed
#' set.seed(123)
#' enroll_rate <- data.frame(rate = 10, duration = 3)
#' fail_rate <- data.frame(
#'   treatment = c("Control", "Experimental"),
#'   rate = c(0.6, 0.4),
#'   dispersion = 0.2
#' )
#' dropout_rate <- data.frame(
#'   treatment = c("Control", "Experimental"),
#'   rate = c(0.05, 0.05),
#'   duration = c(6, 6)
#' )
#' design <- sample_size_nbinom(
#'   lambda1 = 0.6, lambda2 = 0.4, dispersion = 0.2, power = 0.8,
#'   accrual_rate = enroll_rate$rate, accrual_duration = enroll_rate$duration,
#'   trial_duration = 6
#' )
#' cuts <- list(
#'   list(planned_calendar = 2),
#'   list(planned_calendar = 4)
#' )
#' sim_results <- sim_gs_nbinom(
#'   n_sims = 2,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   dropout_rate = dropout_rate,
#'   max_followup = 4,
#'   n_target = 30,
#'   design = design,
#'   cuts = cuts,
#'   seed = TRUE
#' )
#' head(sim_results)
#'
#' \dontrun{
#' # Parallel execution (requires future and future.apply)
#' library(future)
#' plan(multisession, workers = 4)
#' set.seed(42)
#' sim_results <- sim_gs_nbinom(
#'   n_sims = 1000,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   dropout_rate = dropout_rate,
#'   max_followup = 4,
#'   n_target = 30,
#'   design = design,
#'   cuts = cuts,
#'   seed = TRUE
#' )
#' plan(sequential)
#' }
sim_gs_nbinom <- function(
  n_sims, enroll_rate, fail_rate, dropout_rate = NULL,
  max_followup, event_gap = NULL, analysis_times = NULL,
  n_target = NULL, design = NULL,
  data_cut = cut_data_by_date, cuts = NULL,
  seed = TRUE
) {
  # Validate inputs
  if (is.null(design)) {
    stop("design object must be provided to extract planning parameters.")
  }

  if (is.null(cuts)) {
    if (is.null(analysis_times)) {
      stop("Either analysis_times or cuts must be provided.")
    }
    # Create default cuts based on calendar times
    cuts <- lapply(analysis_times, function(t) list(planned_calendar = t))
  }

  n_analyses <- length(cuts)

  # Extract planning parameters for blinded info estimation
  inputs <- if (inherits(design, "gsNB")) design$nb_design$inputs else design$inputs
  lambda1_plan <- inputs$lambda1
  lambda2_plan <- inputs$lambda2
  ratio_plan <- inputs$ratio

  if (is.null(event_gap)) {
    if (!is.null(inputs$event_gap) && !is.na(inputs$event_gap)) {
      event_gap <- inputs$event_gap
    } else {
      event_gap <- 0
    }
  }

  # Function to run one simulation
  run_one_sim <- function(sim_id) {
    # Generate trial data
    sim_data <- nb_sim(
      enroll_rate = enroll_rate,
      fail_rate = fail_rate,
      dropout_rate = dropout_rate,
      max_followup = max_followup,
      n = n_target,
      event_gap = event_gap
    )

    res_list <- vector("list", n_analyses)

    # Keep track of previous cut time to ensure monotonicity if needed
    last_cut_time <- 0

    for (k in seq_len(n_analyses)) {
      # Determine cut date based on criteria
      cut_args <- cuts[[k]]

      # Inject data and parameters if not present (though get_cut_date expects them)
      cut_args$data <- sim_data
      cut_args$event_gap <- event_gap
      cut_args$ratio <- ratio_plan
      cut_args$lambda1 <- lambda1_plan
      cut_args$lambda2 <- lambda2_plan
      cut_args$min_date <- last_cut_time
      # Default max_date to something reasonable if not provided?
      # get_cut_date defaults to Inf.

      cut_time <- do.call(get_cut_date, cut_args)
      last_cut_time <- cut_time

      # Cut data
      cut_data <- data_cut(sim_data, cut_date = cut_time, event_gap = event_gap)

      # Filter enrolled subjects
      enrolled <- unique(sim_data$id[sim_data$enroll_time <= cut_time])
      cut_data <- cut_data[cut_data$id %in% enrolled, ]

      # Summarize counts
      dt <- data.table::as.data.table(cut_data)
      counts <- dt[, .(
        events = sum(events),
        exposure_at_risk = sum(tte),
        exposure_total = sum(tte_total),
        n = .N
      ), by = treatment]

      # Helper to safely extract values (returns 0 if not found)
      get_val <- function(col, trt) {
        val <- counts[treatment == trt][[col]]
        if (length(val) == 0) 0 else val
      }

      n_enrolled <- nrow(cut_data)
      n_ctrl <- get_val("n", "Control")
      n_exp <- get_val("n", "Experimental")
      events_ctrl <- get_val("events", "Control")
      events_exp <- get_val("events", "Experimental")
      exp_at_risk_ctrl <- get_val("exposure_at_risk", "Control")
      exp_at_risk_exp <- get_val("exposure_at_risk", "Experimental")
      exp_total_ctrl <- get_val("exposure_total", "Control")
      exp_total_exp <- get_val("exposure_total", "Experimental")
      events_total <- events_ctrl + events_exp

      z_stat <- NA_real_
      estimate <- NA_real_
      se <- NA_real_
      method_used <- NA_character_
      dispersion <- NA_real_
      
      # Information estimates
      info_unblinded_ml <- NA_real_
      info_blinded_ml <- NA_real_
      info_unblinded_mom <- NA_real_
      info_blinded_mom <- NA_real_
      
      # Legacy columns (mapped to ML)
      blinded_info <- NA_real_
      unblinded_info <- NA_real_

      # Run analysis if sufficient data
      if (n_enrolled >= 4 && events_total >= 2) {
        
        # --- 1. Unblinded ML ---
        test_res <- tryCatch(mutze_test(cut_data), error = function(e) NULL)

        if (!is.null(test_res)) {
          z_stat <- test_res$z
          estimate <- test_res$estimate
          se <- test_res$se
          method_used <- test_res$method
          dispersion <- test_res$dispersion
          
          info_unblinded_ml <- 1 / test_res$se^2
          unblinded_info <- info_unblinded_ml
        }

        # --- 2. Blinded ML ---
        blinded_res <- calculate_blinded_info(
            cut_data,
            ratio = ratio_plan,
            lambda1_planning = lambda1_plan,
            lambda2_planning = lambda2_plan,
            event_gap = event_gap
        )
        info_blinded_ml <- blinded_res$blinded_info
        blinded_info <- info_blinded_ml
        
        # --- Helper for Info Calculation ---
        # Calculates Fisher info for log-rate-ratio given specific parameters and data
        calc_info <- function(dat, lam1, lam2, k) {
          # Unblinded calculation: sum per group
          dat_c <- dat[dat$treatment == "Control", ]
          dat_e <- dat[dat$treatment == "Experimental", ]
          
          # Info = 1 / (1/Info_C + 1/Info_E)
          # Info_group = sum( mu / (1 + k*mu) )
          
          mu_c <- lam1 * dat_c$tte
          term_c <- sum(mu_c / (1 + k * mu_c))
          
          mu_e <- lam2 * dat_e$tte
          term_e <- sum(mu_e / (1 + k * mu_e))
          
          if(term_c <= 0 || term_e <= 0) return(0)
          
          1 / (1/term_c + 1/term_e)
        }
        
        # --- 3. Unblinded MoM ---
        mom_unblinded <- tryCatch(estimate_nb_mom(cut_data, group = "treatment"), error = function(e) NULL)
        if (!is.null(mom_unblinded) && !any(is.na(mom_unblinded$lambda))) {
             # lambda is named vector c("Control" = ..., "Experimental" = ...)
             # Ensure we map correctly
             l1_mom <- mom_unblinded$lambda["Control"]
             l2_mom <- mom_unblinded$lambda["Experimental"]
             k_mom <- mom_unblinded$dispersion
             
             # If names missing or different, try positional (usually factor order Control, Exp)
             if(is.na(l1_mom)) l1_mom <- mom_unblinded$lambda[1]
             if(is.na(l2_mom)) l2_mom <- mom_unblinded$lambda[2]
             
             info_unblinded_mom <- calc_info(cut_data, l1_mom, l2_mom, k_mom)
        }
        
        # --- 4. Blinded MoM ---
        mom_blinded <- tryCatch(estimate_nb_mom(cut_data), error = function(e) NULL)
        if (!is.null(mom_blinded)) {
             l_agg <- mom_blinded$lambda
             k_agg <- mom_blinded$dispersion
             
             # Split lambda using planning ratio (similar to calculate_blinded_info)
             # Assumption: lambda_agg is mixture of l1 and l2 weighted by sample size/allocation
             p1 <- 1 / (1 + ratio_plan)
             p2 <- ratio_plan / (1 + ratio_plan)
             rr_plan <- lambda2_plan / lambda1_plan
             
             l1_adj <- l_agg / (p1 + p2 * rr_plan)
             l2_adj <- l1_adj * rr_plan
             
             # For Blinded Info, we treat the data as if we don't know the treatment,
             # effectively calculating the expected information.
             # Using the .blinded_info_from_tte logic locally:
             
             mu1 <- l1_adj * cut_data$tte
             mu2 <- l2_adj * cut_data$tte
             
             w1 <- p1 * sum(mu1 / (1 + k_agg * mu1))
             w2 <- p2 * sum(mu2 / (1 + k_agg * mu2))
             
             if (w1 > 0 && w2 > 0) {
                var_log_rr <- 1/w1 + 1/w2
                info_blinded_mom <- 1 / var_log_rr
             } else {
                info_blinded_mom <- 0
             }
        }
      }

      res_list[[k]] <- data.frame(
        sim = sim_id,
        analysis = k,
        analysis_time = cut_time,
        n_enrolled = n_enrolled,
        n_ctrl = n_ctrl,
        n_exp = n_exp,
        events_total = events_total,
        events_ctrl = events_ctrl,
        events_exp = events_exp,
        exposure_at_risk_ctrl = exp_at_risk_ctrl,
        exposure_at_risk_exp = exp_at_risk_exp,
        exposure_total_ctrl = exp_total_ctrl,
        exposure_total_exp = exp_total_exp,
        z_stat = z_stat,
        estimate = estimate,
        se = se,
        method_used = method_used,
        dispersion = dispersion,
        blinded_info = blinded_info,
        unblinded_info = unblinded_info, # Legacy
        info_unblinded_ml = info_unblinded_ml,
        info_blinded_ml = info_blinded_ml,
        info_unblinded_mom = info_unblinded_mom,
        info_blinded_mom = info_blinded_mom
      )
    }
    do.call(rbind, res_list)
  }

  # Run simulations (parallel-safe via future.apply if available)
  has_future_apply <- requireNamespace("future.apply", quietly = TRUE)

  if (has_future_apply) {
    results_list <- future.apply::future_lapply(
      seq_len(n_sims),
      run_one_sim,
      future.seed = seed
    )
  } else {
    if (is.numeric(seed)) set.seed(seed)
    results_list <- lapply(seq_len(n_sims), run_one_sim)
  }

  do.call(rbind, results_list)
}
