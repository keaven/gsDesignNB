#' Simulate Group Sequential Clinical Trial for Negative Binomial Outcomes
#'
#' Simulates multiple replicates of a group sequential clinical trial with negative
#' binomial outcomes, performing interim analyses at specified calendar times.
#'
#' @param n_sims Number of simulations to run.
#' @param enroll_rate Enrollment rates (data frame with `rate` and `duration`).
#' @param fail_rate Failure rates (data frame with `treatment`, `rate`, `dispersion`).
#' @param dropout_rate Dropout rates (data frame with `treatment`, `rate`, `duration`).
#' @param max_followup Maximum follow-up time.
#' @param event_gap Event gap duration.
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
#'
#' @return A data frame containing simulation results for each analysis of each trial.
#'   Columns include:
#'   \describe{
#'     \item{sim}{Simulation ID}
#'     \item{analysis}{Analysis index}
#'     \item{analysis_time}{Calendar time of analysis}
#'     \item{n_enrolled}{Number of subjects enrolled}
#'     \item{events_total}{Total events observed}
#'     \item{events_ctrl}{Events in control group}
#'     \item{events_exp}{Events in experimental group}
#'     \item{exposure_ctrl}{Total exposure in control group}
#'     \item{exposure_exp}{Total exposure in experimental group}
#'     \item{z_stat}{Z-statistic from the Wald test (positive favors experimental if rate ratio < 1)}
#'     \item{blinded_info}{Estimated blinded statistical information}
#'     \item{unblinded_info}{Observed unblinded statistical information}
#'   }
#'
#' @export
#' @importFrom data.table as.data.table
sim_gs_nbinom <- function(n_sims, enroll_rate, fail_rate, dropout_rate = NULL,
                          max_followup, event_gap = 0, analysis_times = NULL,
                          n_target = NULL, design = NULL,
                          data_cut = cut_data_by_date, cuts = NULL) {
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
      counts <- dt[, .(events = sum(events), exposure = sum(tte)), by = treatment]
      
      n_enrolled <- nrow(cut_data)
      events_ctrl <- sum(counts[treatment == "Control"]$events)
      events_exp <- sum(counts[treatment == "Experimental"]$events)
      exp_ctrl <- sum(counts[treatment == "Control"]$exposure)
      exp_exp <- sum(counts[treatment == "Experimental"]$exposure)
      events_total <- events_ctrl + events_exp

      z_stat <- NA_real_
      blinded_info <- NA_real_
      unblinded_info <- NA_real_

      # Run analysis if sufficient data
      if (n_enrolled >= 4 && events_total >= 2) {
        test_res <- tryCatch(mutze_test(cut_data), error = function(e) NULL)
        
        if (!is.null(test_res)) {
          z_stat <- test_res$z
          unblinded_info <- 1 / test_res$se^2
          
          # Blinded info estimation
          blinded_res <- calculate_blinded_info(
            cut_data,
            ratio = ratio_plan,
            lambda1_planning = lambda1_plan,
            lambda2_planning = lambda2_plan,
            event_gap = event_gap
          )
          blinded_info <- blinded_res$blinded_info
        }
      }

      res_list[[k]] <- data.frame(
        sim = sim_id,
        analysis = k,
        analysis_time = cut_time,
        n_enrolled = n_enrolled,
        events_total = events_total,
        events_ctrl = events_ctrl,
        events_exp = events_exp,
        exposure_ctrl = exp_ctrl,
        exposure_exp = exp_exp,
        z_stat = z_stat,
        blinded_info = blinded_info,
        unblinded_info = unblinded_info
      )
    }
    do.call(rbind, res_list)
  }

  # Run simulations
  results_list <- lapply(seq_len(n_sims), run_one_sim)
  do.call(rbind, results_list)
}
