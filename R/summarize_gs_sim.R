#' Summarize group sequential simulation results
#'
#' Provides a summary of the operating characteristics of the group sequential design
#' based on simulation results.
#'
#' @param x A data frame returned by [check_gs_bound()] (or [sim_gs_nbinom()] if
#'   bounds are manually checked). Must contain columns `cross_upper`, `cross_lower`.
#' @param info_trim Proportion of observations trimmed from each tail when
#'   summarizing information estimates. Defaults to `0.01` to reduce sensitivity
#'   to occasional unstable fitted information estimates in small interim data
#'   sets.
#'
#' @return A list containing:
#'   \describe{
#'     \item{n_sim}{Number of simulations}
#'     \item{power}{Overall power (probability of crossing upper bound)}
#'     \item{futility}{Overall futility rate (probability of crossing lower bound and not upper)}
#'     \item{analysis_summary}{Data frame with per-analysis statistics (sample
#'       size, events, information, crossings, and optional exposure columns when
#'       present in `x`).}
#'   }
#'
#' @importFrom data.table as.data.table
#'
#' @export
#'
#' @examples
#' design <- gsDesign::gsDesign(k = 2, n.fix = 80, test.type = 2, timing = c(0.5, 1))
#' sim_df <- data.frame(
#'   sim = c(1, 1, 2, 2),
#'   analysis = c(1, 2, 1, 2),
#'   z_stat = c(2.4, NA, -0.5, 1.9),
#'   blinded_info = c(40, 80, 40, 80),
#'   unblinded_info = c(40, 80, 40, 80),
#'   n_enrolled = c(30, 60, 30, 60),
#'   events_total = c(12, 25, 10, 22)
#' )
#' bounds_checked <- check_gs_bound(sim_df, design)
#' summarize_gs_sim(bounds_checked)
summarize_gs_sim <- function(x, info_trim = 0.01) {
  dt <- data.table::as.data.table(x)

  if (!all(c("cross_upper", "cross_lower") %in% names(dt))) {
    stop("Input must contain 'cross_upper' and 'cross_lower' columns. Run check_gs_bound() first.")
  }

  if (!is.numeric(info_trim) || length(info_trim) != 1 ||
      is.na(info_trim) || info_trim < 0 || info_trim >= 0.5) {
    stop("info_trim must be a single number in [0, 0.5).")
  }

  finite_mean <- function(x, trim = 0) {
    x <- x[is.finite(x) & !is.na(x)]
    if (length(x) == 0) {
      return(NA_real_)
    }
    mean(x, trim = trim)
  }

  n_sims <- length(unique(dt$sim))

  # Overall Power: Did any analysis cross upper?
  power_by_sim <- dt[, .(success = any(cross_upper)), by = sim]
  overall_power <- mean(power_by_sim$success)

  # Overall Futility: Did any analysis cross lower (and NO upper)?
  # Note: usually crossing lower stops the trial, so subsequent upper is impossible.
  # But if simulation continued, we check.
  futility_by_sim <- dt[, .(futility = any(cross_lower) & !any(cross_upper)), by = sim]
  overall_futility <- mean(futility_by_sim$futility)

  # Per-Analysis Summary
  analysis_summary <- dt[, .(
    n_enrolled = finite_mean(n_enrolled),
    events = finite_mean(events_total),
    info_blinded = finite_mean(blinded_info, trim = info_trim),
    info_unblinded = finite_mean(unblinded_info, trim = info_trim),
    n_cross_upper = sum(cross_upper, na.rm = TRUE),
    n_cross_lower = sum(cross_lower, na.rm = TRUE),
    prob_cross_upper = sum(cross_upper, na.rm = TRUE) / n_sims,
    prob_cross_lower = sum(cross_lower, na.rm = TRUE) / n_sims
  ), by = analysis]

  optional_mean_cols <- c(
    "n_ctrl", "n_exp",
    "events_ctrl", "events_exp",
    "exposure_total_ctrl", "exposure_total_exp",
    "exposure_at_risk_ctrl", "exposure_at_risk_exp",
    "analysis_time"
  )

  for (col in optional_mean_cols[optional_mean_cols %in% names(dt)]) {
    col_summary <- dt[, .(value = finite_mean(get(col))), by = analysis]
    data.table::setnames(col_summary, "value", col)
    analysis_summary <- merge(
      analysis_summary,
      col_summary,
      by = "analysis",
      all.x = TRUE,
      sort = FALSE
    )
  }

  # Cumulative crossing probabilities
  analysis_summary[, cum_prob_upper := cumsum(prob_cross_upper)]

  preferred_order <- c(
    "analysis", "analysis_time", "n_enrolled", "n_ctrl", "n_exp",
    "exposure_total_ctrl", "exposure_total_exp",
    "exposure_at_risk_ctrl", "exposure_at_risk_exp",
    "events", "events_ctrl", "events_exp",
    "info_blinded", "info_unblinded",
    "n_cross_upper", "n_cross_lower",
    "prob_cross_upper", "prob_cross_lower", "cum_prob_upper"
  )
  data.table::setcolorder(
    analysis_summary,
    c(intersect(preferred_order, names(analysis_summary)),
      setdiff(names(analysis_summary), preferred_order))
  )

  list(
    n_sim = n_sims,
    power = overall_power,
    futility = overall_futility,
    analysis_summary = as.data.frame(analysis_summary)
  )
}
