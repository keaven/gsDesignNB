#' Check group sequential bounds
#'
#' Updates the group sequential design boundaries based on observed information
#' and checks if boundaries have been crossed.
#'
#' @param sim_results Data frame of simulation results (from [sim_gs_nbinom()]).
#' @param design The planning `gsNB` object.
#' @param info_scale Character. Legacy selector for `"blinded"` (default) or
#'   `"unblinded"` information. Ignored when `info_col` is supplied.
#' @param info_col Optional explicit column name containing the information
#'   metric to use for bounds, e.g. `"info_unblinded_ml"` or
#'   `"info_blinded_mom"`.
#' @param boundary_update `"fixed"` or `"recompute"`. `NULL` preserves the legacy path.
#' @param spending Spending strategy as in [gsNBPower()]. `NULL` preserves legacy
#'   behavior unless another monitoring option is supplied.
#' @param full_spending_at_final Complete error spending at the final analysis.
#' @param final_info_col Projected-final information column for `"expected_final"`.
#' @param ignore_futility Ignore futility stops when evaluating nonbinding Type I error.
#'
#' @return A data frame with added columns:
#'   \describe{
#'     \item{cross_upper}{Logical, true if upper bound crossed (efficacy)}
#'     \item{cross_lower}{Logical, true if lower bound crossed (futility)}
#'     \item{cross_harm}{Logical, true if harm bound crossed (test.type 7 or 8)}
#'   }
#'
#' @export
#'
#' @examples
#' design <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 2, timing = c(0.5, 1))
#' sim_df <- data.frame(
#'   sim = c(1, 1, 2, 2),
#'   analysis = c(1, 2, 1, 2),
#'   z_stat = c(2.5, NA, -0.2, 2.2),
#'   blinded_info = c(50, 100, 50, 100),
#'   unblinded_info = c(50, 100, 50, 100)
#' )
#' check_gs_bound(sim_df, design)
#' check_gs_bound(sim_df, design, info_col = "unblinded_info")
check_gs_bound <- function(
  sim_results,
  design,
  info_scale = c("blinded", "unblinded"),
  info_col = NULL, boundary_update = NULL, spending = NULL,
  full_spending_at_final = TRUE,
  final_info_col = "info_blinded_final_projected", ignore_futility = FALSE
) {
  cross_harm <- NULL
  info_scale <- match.arg(info_scale)

  if (!inherits(design, "gsNB") && !inherits(design, "gsDesign")) {
    stop("design must be a gsDesign or gsNB object")
  }

  dt <- data.table::as.data.table(sim_results)

  # Ensure we have z_stat
  if (!"z_stat" %in% names(dt)) stop("sim_results must contain 'z_stat'")

  # Identify information column
  if (is.null(info_col)) {
    info_col <- if (info_scale == "blinded") "blinded_info" else "unblinded_info"
  }
  if (!info_col %in% names(dt)) stop(paste("sim_results must contain", info_col))
  if (!is.null(spending) || !is.null(boundary_update) || ignore_futility ||
      !full_spending_at_final) {
    return(.gs_nb_check_bounds(
      as.data.frame(dt), design, info_col,
      if (is.null(boundary_update)) "recompute" else boundary_update,
      if (is.null(spending)) "min_planned_actual" else spending,
      full_spending_at_final, final_info_col, ignore_futility
    ))
  }

  # Process each simulation separately
  process_sim <- function(sub_dt) {
    k_max <- max(sub_dt$analysis)
    max_info <- design$n.fix

    cross_upper <- rep(FALSE, k_max)
    cross_lower <- rep(FALSE, k_max)
    cross_harm <- rep(FALSE, k_max)
    stopped <- FALSE
    has_harm <- design$test.type %in% c(7, 8)

    # Current timing vector (starts with design timing)
    current_timing <- design$timing

    for (k in seq_len(k_max)) {
      if (stopped) {
        # Already stopped
        break
      }

      obs_info <- sub_dt[analysis == k][[info_col]]
      z_val <- sub_dt[analysis == k]$z_stat

      if (is.na(z_val)) next

      # Flip Z if needed?
      # In vignette: z_eff <- -test_result$z.
      # mutze_test returns Z for treatment coefficient.
      # If lambda2 < lambda1 (benefit), coeff is negative, Z is negative.
      # gsDesign assumes positive Z for efficacy.
      # So we should flip Z if the design anticipates negative effect.
      # The user vignette flips it. I will assume we flip it.
      # Ideally, we check delta1.
      if (!is.null(design$delta1) && design$delta1 < 0) {
        z_val <- -z_val
      }

      # Calculate information fraction
      frac <- obs_info / max_info

      # Cap fraction at 1 if overrunning
      if (frac > 1) frac <- 1

      # Update timing
      current_timing[k] <- frac

      # Ensure monotonicity
      if (k > 1 && current_timing[k] <= current_timing[k - 1]) {
        current_timing[k] <- current_timing[k - 1] + 0.001
      }
      if (k < k_max && current_timing[k + 1] <= current_timing[k]) {
        # Push future timings if needed (though they will be overwritten later)
        current_timing[k + 1] <- min(current_timing[k] + 0.001, 0.999)
      }

      # Handle final analysis logic
      if (k == k_max) current_timing[k] <- 1

      # Update design bounds
      # We need to recreate the design with the new timing
      # We extract parameters from the original design
      # Note: usTime is stored in gsNB objects (from my recent fix)
      # If present, we might want to respect it, BUT the vignette logic
      # explicitly overwrites timing based on information fraction.
      # Standard error spending usually adapts to information.
      # If usTime was used, gsDesign usually ignores timing?
      # Wait, if usTime is provided, gsDesign uses it for spending.
      # If we want to use Information Fraction spending (Lan-DeMets), we should NOT pass usTime?
      # The vignette recalculates bounds:
      # gsDesign(..., timing = current_timing, ...)
      # It passes usTime too!
      # If usTime is passed, gsDesign uses usTime for spending calculation.
      # So changing `timing` (information fraction) only affects the drift/B-value scaling?
      # Actually, boundaries depend on Information Fraction for covariance.
      # Spending depends on Spending Function argument `t`.
      # If `usTime` is provided, `t = usTime`. If not, `t = timing`.

      # To follow vignette logic:

      # It passes `usTime` AND `timing`.
      # Use tryCatch to handle edge cases where spending function can't
      # achieve target beta with the observed information fraction.
      gs_args <- list(
          k = design$k,
          test.type = design$test.type,
          alpha = design$alpha,
          beta = design$beta,
          astar = design$astar,
          delta = design$delta,
          sfu = design$upper$sf,
          sfupar = design$upper$param,
          sfl = design$lower$sf,
          sflpar = design$lower$param,
          tol = design$tol,
          r = design$r,
          n.fix = max_info,
          timing = current_timing,
          usTime = design$usTime,
          lsTime = design$lsTime
      )
      if (design$test.type %in% c(7, 8) && !is.null(design$harm)) {
          gs_args$sfharm <- design$harm$sf
          gs_args$sfharmparam <- design$harm$param
      }
      if (!is.null(design$testUpper)) gs_args$testUpper <- design$testUpper
      if (!is.null(design$testLower)) gs_args$testLower <- design$testLower
      if (!is.null(design$testHarm)) gs_args$testHarm <- design$testHarm

      temp_gs <- tryCatch(
        do.call(gsDesign::gsDesign, gs_args),
        error = function(e) {
          # If the spending function fails, fall back to original design bounds
          # This can happen with extreme information fractions
          design
        }
      )

      upper_bound <- if (length(temp_gs$upper$bound) >= k) {
        temp_gs$upper$bound[k]
      } else {
        Inf
      }
      lower_bound <- if (length(temp_gs$lower$bound) >= k) {
        temp_gs$lower$bound[k]
      } else {
        -Inf
      }

      if (is.finite(upper_bound) && z_val > upper_bound) {
        cross_upper[k] <- TRUE
        stopped <- TRUE
      } else if (has_harm && !is.null(temp_gs$harm) && z_val < temp_gs$harm$bound[k]) {
        cross_harm[k] <- TRUE
        stopped <- TRUE
      } else if (is.finite(lower_bound) && z_val < lower_bound) {
        cross_lower[k] <- TRUE
        stopped <- TRUE
      }
    }

    return(list(cross_upper = cross_upper, cross_lower = cross_lower,
                cross_harm = cross_harm))
  }

  # Apply to all sims
  # This might be slow for many sims if done in R loop.
  # But bounds update is fast.

  results <- dt[, process_sim(.SD), by = sim]

  # Merge back
  dt[, cross_upper := results$cross_upper]
  dt[, cross_lower := results$cross_lower]
  dt[, cross_harm := results$cross_harm]

  as.data.frame(dt)
}

.gs_nb_check_bounds <- function(
  data, design, info_col, policy, spending, full_final, final_info_col, ignore_futility
) {
  policy <- match.arg(policy, c("fixed", "recompute"))
  spending <- match.arg(spending, c("min_planned_actual", "planned_final", "expected_final", "calendar"))
  if (!all(c("sim", "analysis") %in% names(data))) {
    stop("sim_results must contain sim and analysis.", call. = FALSE)
  }
  if (any(!is.finite(data$analysis)) || any(data$analysis < 1) ||
      any(data$analysis > design$k) || any(data$analysis != floor(data$analysis))) {
    stop("analysis must identify a planned look.", call. = FALSE)
  }
  if (spending == "expected_final" && policy != "fixed" && !final_info_col %in% names(data)) {
    stop("sim_results must contain ", final_info_col, ".", call. = FALSE)
  }
  data$cross_upper <- data$cross_lower <- data$cross_harm <- FALSE
  data$upper_bound <- data$lower_bound <- data$harm_bound <- NA_real_
  data$spending_time <- data$lower_spending_time <- NA_real_
  data$bound_status <- "not_evaluated"
  data$bound_error <- ""
  groups <- split(seq_len(nrow(data)), data$sim)
  for (rows in groups) {
    rows <- rows[order(data$analysis[rows])]
    if (anyDuplicated(data$analysis[rows])) stop("Each simulation needs unique analysis indices.", call. = FALSE)
    stopped <- FALSE
    previous_ids <- previous_info <- previous_dates <- upper_history <- lower_history <- numeric()
    skipped_ids <- if ("status" %in% names(data)) data$analysis[rows[data$status[rows] == "skipped"]] else numeric()
    for (row in rows) {
      look <- data$analysis[row]
      if (look %in% skipped_ids) {
        data$bound_status[row] <- "skipped"
        next
      }
      if (stopped) next
      observed <- data[[info_col]][row]
      z_value <- data$z_stat[row]
      if (!is.finite(z_value) || (policy != "fixed" && (!is.finite(observed) || observed <= 0))) {
        data$bound_status[row] <- "missing_information_or_statistic"
        next
      }
      date <- if ("analysis_time" %in% names(data)) data$analysis_time[row] else
        if (!is.null(design$T)) design$T[look] else look
      future_ids <- seq_len(design$k)[seq_len(design$k) > look]
      ids <- c(previous_ids, look, future_ids)
      planned_final <- utils::tail(design$n.I, 1)
      final_estimate <- if (spending == "expected_final" && policy != "fixed") {
        data[[final_info_col]][row]
      } else max(planned_final, observed + planned_final - design$n.I[look])
      remaining <- planned_final - design$n.I[look]
      future_info <- if (length(future_ids) && remaining > 0) {
        observed + (final_estimate - observed) *
          (design$n.I[future_ids] - design$n.I[look]) / remaining
      } else numeric()
      information <- c(previous_info, observed, future_info)
      future_dates <- if (is.null(design$T)) future_ids else design$T[future_ids]
      dates <- c(previous_dates, date, future_dates)
      bounds <- tryCatch({
        if (policy == "fixed") {
          list(upper = design$upper$bound[ids],
               lower = if (is.null(design$lower)) rep(-Inf, length(ids)) else design$lower$bound[ids],
               harm = if (is.null(design$harm)) rep(-Inf, length(ids)) else design$harm$bound[ids],
               us_time = rep(NA_real_, length(ids)), ls_time = rep(NA_real_, length(ids)))
        } else {
          if (any(!is.finite(information)) || any(information <= 0) ||
              any(diff(information) <= 0)) {
            stop("Monitoring information is not positive and increasing.")
          }
          do.call(.gs_nb_bounds, list(
            design = design, info = information, dates = dates, look_ids = ids,
            policy = policy, spending = spending, full_final = full_final,
            final_info = final_estimate,
            history = list(upper = upper_history, lower = lower_history)
          ))
        }
      }, error = function(error) error)
      if (inherits(bounds, "error")) {
        data$bound_status[row] <- "failed"
        data$bound_error[row] <- conditionMessage(bounds)
        next
      }
      position <- length(previous_ids) + 1L
      data$upper_bound[row] <- bounds$upper[position]
      data$lower_bound[row] <- bounds$lower[position]
      data$harm_bound[row] <- bounds$harm[position]
      data$spending_time[row] <- bounds$us_time[position]
      data$lower_spending_time[row] <- bounds$ls_time[position]
      data$bound_status[row] <- "ok"
      if (!is.null(design$delta1) && design$delta1 < 0) z_value <- -z_value
      if (z_value > data$upper_bound[row]) {
        data$cross_upper[row] <- TRUE
        stopped <- TRUE
      } else if (design$test.type %in% c(7, 8) && z_value < data$harm_bound[row]) {
        data$cross_harm[row] <- TRUE
        stopped <- TRUE
      } else if (!ignore_futility && z_value < data$lower_bound[row]) {
        data$cross_lower[row] <- TRUE
        stopped <- TRUE
      }
      previous_ids <- c(previous_ids, look)
      previous_info <- c(previous_info, observed)
      previous_dates <- c(previous_dates, date)
      upper_history <- c(upper_history, data$spending_time[row])
      lower_history <- c(lower_history, data$lower_spending_time[row])
    }
  }
  data
}
