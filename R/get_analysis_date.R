#' Find calendar date for target event count
#'
#' Finds the calendar time (since start of randomization) at which a specified
#' total number of events is reached in the simulated dataset.
#'
#' @param data A data frame of simulated data, typically from [nb_sim()].
#' @param planned_events Integer. The target number of events.
#' @param event_gap Gap duration after each event during which no new events are counted.
#'   Can be a numeric value (default `5 / 365.25`) or a function returning a numeric value.
#'
#' @return Numeric. The calendar date when `planned_events` is achieved.
#'   If the dataset contains fewer than `planned_events`, returns the maximum
#'   calendar time in the dataset and prints a message.
#'
#' @examples
#' enroll_rate <- data.frame(rate = 20 / (5 / 12), duration = 5 / 12)
#' fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3))
#' dropout_rate <- data.frame(
#'   treatment = c("Control", "Experimental"),
#'   rate = c(0.1, 0.05), duration = c(100, 100)
#' )
#' sim <- nb_sim(enroll_rate, fail_rate, dropout_rate, max_followup = 2, n = 40)
#' get_analysis_date(sim, planned_events = 15)
#' @export
get_analysis_date <- function(data, planned_events, event_gap = 5 / 365.25) {
  dt <- data.table::as.data.table(data)

  # Filter for actual events (event == 1)
  events_dt <- dt[event == 1]

  # If no events at all
  if (nrow(events_dt) == 0) {
    message(sprintf("Only 0 events in trial"))
    return(max(dt$calendar_time))
  }

  # Identify valid events respecting the gap
  get_valid_indices <- function(tte_vec, gap_rule) {
    if (length(tte_vec) == 0) {
      return(integer(0))
    }

    # Sort by tte (time relative to randomization)
    ord <- order(tte_vec)
    sorted_tte <- tte_vec[ord]

    keep_logical <- logical(length(sorted_tte))
    last_valid_end <- -Inf

    for (i in seq_along(sorted_tte)) {
      t <- sorted_tte[i]
      if (t < last_valid_end) next # Skip if in gap

      keep_logical[i] <- TRUE

      g <- if (is.function(gap_rule)) gap_rule() else gap_rule
      last_valid_end <- t + g
    }

    ord[keep_logical]
  }

  valid_rows <- events_dt[, .SD[get_valid_indices(tte, event_gap)], by = id]
  total_events <- nrow(valid_rows)

  if (total_events < planned_events) {
    message(sprintf("Only %d events in trial", total_events))
    # Return max calendar time in the entire dataset (including censoring times)
    return(max(data$calendar_time))
  }

  # Sort valid events by calendar time
  valid_times <- sort(valid_rows$calendar_time)

  # The calendar time of the planned_events-th event
  return(valid_times[planned_events])
}
