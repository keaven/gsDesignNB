#' Summarize group sequential simulation results
#'
#' Provides a summary of the operating characteristics of the group sequential design
#' based on simulation results.
#'
#' @param x A data frame returned by [check_gs_bound()] (or [sim_gs_nbinom()] if
#'   bounds are manually checked). Must contain columns `cross_upper`, `cross_lower`.
#'
#' @return A list containing:
#'   \describe{
#'     \item{n_sim}{Number of simulations}
#'     \item{power}{Overall power (probability of crossing upper bound)}
#'     \item{futility}{Overall futility rate (probability of crossing lower bound and not upper)}
#'     \item{analysis_summary}{Data frame with per-analysis statistics (events, crossings)}
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
summarize_gs_sim <- function(x) {
  dt <- data.table::as.data.table(x)

  if (!all(c("cross_upper", "cross_lower") %in% names(dt))) {
    stop("Input must contain 'cross_upper' and 'cross_lower' columns. Run check_gs_bound() first.")
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
    n_enrolled = mean(n_enrolled, na.rm = TRUE),
    events = mean(events_total, na.rm = TRUE),
    info_blinded = mean(blinded_info, na.rm = TRUE),
    info_unblinded = mean(unblinded_info, na.rm = TRUE),
    n_cross_upper = sum(cross_upper, na.rm = TRUE),
    n_cross_lower = sum(cross_lower, na.rm = TRUE),
    prob_cross_upper = sum(cross_upper, na.rm = TRUE) / n_sims,
    prob_cross_lower = sum(cross_lower, na.rm = TRUE) / n_sims
  ), by = analysis]

  # Cumulative crossing probabilities
  analysis_summary[, cum_prob_upper := cumsum(prob_cross_upper)]

  list(
    n_sim = n_sims,
    power = overall_power,
    futility = overall_futility,
    analysis_summary = as.data.frame(analysis_summary)
  )
}
