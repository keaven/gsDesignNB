#' Determine analysis date based on criteria
#'
#' Finds the earliest calendar date at which all specified criteria are met.
#' Criteria can include a specific calendar date, a target number of events,
#' a target number of completers, or a target amount of blinded information.
#'
#' @param data A data frame of simulated data (from [nb_sim()]).
#' @param planned_calendar Numeric. Target calendar time.
#' @param target_events Integer. Target number of observed events.
#' @param target_completers Integer. Target number of subjects with complete follow-up.
#' @param target_info Numeric. Target blinded information.
#' @param event_gap Numeric. Gap duration for event counting and info calculation.
#' @param ratio Numeric. Randomization ratio (experimental/control) for info calculation.
#' @param lambda1 Numeric. Planned control rate for info calculation.
#' @param lambda2 Numeric. Planned experimental rate for info calculation.
#' @param min_date Numeric. Minimum possible date (e.g., 0 or previous analysis time).
#' @param max_date Numeric. Maximum possible date (e.g., trial duration).
#'
#' @return Numeric. The calendar date satisfying the criteria. If criteria cannot be met
#'   within `max_date` (or data limits), returns `max_date` (or max data time).
#'
#' @importFrom stats uniroot
#'
#' @export
#'
#' @examples
#' set.seed(456)
#' enroll_rate <- data.frame(rate = 15, duration = 1)
#' fail_rate <- data.frame(
#'   treatment = c("Control", "Experimental"),
#'   rate = c(0.6, 0.4)
#' )
#' sim_data <- nb_sim(enroll_rate, fail_rate, max_followup = 1, n = 20)
#' get_cut_date(sim_data, planned_calendar = 0.5, target_events = 5, event_gap = 0)
get_cut_date <- function(
  data, planned_calendar = NULL, target_events = NULL,
  target_completers = NULL, target_info = NULL,
  event_gap = 0, ratio = 1, lambda1 = NULL, lambda2 = NULL,
  min_date = 0, max_date = Inf
) {
  dates <- numeric(0)

  # 1. Calendar Time
  if (!is.null(planned_calendar)) {
    dates <- c(dates, planned_calendar)
  }

  # 2. Events
  if (!is.null(target_events)) {
    d_events <- get_analysis_date(data, planned_events = target_events, event_gap = event_gap)
    dates <- c(dates, d_events)
  }

  # 3. Completers
  if (!is.null(target_completers)) {
    d_completers <- cut_date_for_completers(data, target_completers = target_completers)
    dates <- c(dates, d_completers)
  }

  # 4. Information
  if (!is.null(target_info)) {
    if (is.null(lambda1) || is.null(lambda2)) {
      stop("lambda1 and lambda2 must be provided for information targeting.")
    }

    # Define function to calculate info at time t
    calc_info_at_t <- function(t) {
      cut_dt <- cut_data_by_date(data, cut_date = t, event_gap = event_gap)
      # Check if sufficient data for estimation
      if (sum(cut_dt$events) < 2) {
        return(0)
      }

      res <- tryCatch(
        calculate_blinded_info(
          cut_dt,
          ratio = ratio, lambda1_planning = lambda1, lambda2_planning = lambda2,
          event_gap = event_gap
        ),
        error = function(e) list(blinded_info = 0)
      )
      res$blinded_info
    }

    # Search for time t where info >= target_info
    # Bound search by min_date and max_date
    # Check max first
    info_max <- calc_info_at_t(max_date)
    if (info_max < target_info) {
      d_info <- max_date
    } else {
      # Check min
      info_min <- calc_info_at_t(min_date)
      if (info_min >= target_info) {
        d_info <- min_date
      } else {
        # Root finding
        # Function: calc_info_at_t(t) - target_info
        # Note: Information is not strictly continuous/monotonic due to discreteness of events,
        # but broadly monotonic. uniroot should work ok.
        root_res <- tryCatch(
          uniroot(
            function(t) calc_info_at_t(t) - target_info,
            interval = c(min_date, max_date),
            extendInt = "upX" # Allow extending up if needed (though we checked max)
          ),
          error = function(e) NULL
        )

        if (!is.null(root_res)) {
          d_info <- root_res$root
        } else {
          # Fallback: linear interpolation or grid search?
          # Since we checked bounds, uniroot failing implies non-monotonicity or error.
          # Return max_date conservatively.
          d_info <- max_date
        }
      }
    }
    dates <- c(dates, d_info)
  }

  if (length(dates) == 0) {
    return(max_date)
  }

  # Take the maximum of required dates (must satisfy ALL criteria)
  final_date <- max(dates, na.rm = TRUE)

  # Constrain to limits
  final_date <- max(min_date, min(final_date, max_date))

  return(final_date)
}
