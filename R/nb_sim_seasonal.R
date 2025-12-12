#' Simulate recurrent events with seasonal rates
#'
#' Simulates recurrent events where event rates depend on the season.
#'
#' @param enroll_rate A data frame with columns `rate` and `duration`.
#' @param fail_rate A data frame with columns `treatment`, `season`, `rate`, and optionally `dispersion`.
#'   Seasons should be "Spring", "Summer", "Fall", "Winter".
#' @param dropout_rate A data frame with columns `treatment`, `rate`, `duration`.
#' @param max_followup Numeric. Max follow-up duration (years).
#' @param randomization_start_date Date. Start of randomization.
#' @param n Integer. Total sample size.
#' @param block Character vector for block randomization.
#'
#' @return A data frame of class `nb_sim_seasonal` with columns:
#'   `id`, `treatment`, `season`, `enroll_time`, `start`, `end`, `event`, `calendar_start`, `calendar_end`.
#'   Rows represent intervals of risk or events. `event=1` indicates an event at `end`.
#'   `event=0` indicates censoring or end of a seasonal interval at `end`.
#'
#' @import data.table
#' @export
nb_sim_seasonal <- function(enroll_rate, fail_rate, dropout_rate = NULL, max_followup = NULL,
                            randomization_start_date = NULL, n = NULL,
                            block = c(rep("Control", 2), rep("Experimental", 2))) {
  if (is.null(randomization_start_date)) stop("randomization_start_date is required")
  start_date <- as.Date(randomization_start_date)

  # Reuse nb_sim logic for enrollment and randomization
  # (Copying relevant parts or refactoring would be better, but for now we inline/copy to ensure independence)

  if (is.null(enroll_rate) || !is.data.frame(enroll_rate)) stop("enroll_rate must be a data frame")
  if (is.null(fail_rate) || !is.data.frame(fail_rate)) stop("fail_rate must be a data frame")
  if (is.null(max_followup)) stop("max_followup must be provided")

  if (is.null(n)) {
    n <- round(sum(enroll_rate$rate * enroll_rate$duration))
  }

  enroll_times <- simtrial::rpwexp_enroll(n = n, enroll_rate = enroll_rate)
  enroll_times <- sort(enroll_times)[seq_len(n)]

  treatments <- unique(fail_rate$treatment)

  # Block randomization
  if (is.null(block)) {
    block <- rep(treatments, each = 2)
  }
  if (!all(block %in% treatments)) {
    # Fallback if default block doesn't match treatments
    block <- rep(treatments, each = 2)
  }

  n_blocks <- ceiling(n / length(block))
  full_blocks <- replicate(n_blocks, sample(block))
  assigned_trt <- as.vector(full_blocks)[seq_len(n)]

  # Helper to get dispersion by treatment
  # Assumes dispersion is constant per treatment (not seasonal)
  # Look for dispersion in fail_rate.
  get_dispersion <- function(trt) {
    sub <- fail_rate[fail_rate$treatment == trt, ]
    if ("dispersion" %in% names(sub)) {
      d <- unique(sub$dispersion)
      if (length(d) > 1) stop("Dispersion must be constant within treatment")
      if (length(d) == 1) {
        return(d)
      }
    }
    return(0)
  }

  # Pre-calculate dispersions
  dispersions <- sapply(unique(assigned_trt), get_dispersion)
  names(dispersions) <- unique(assigned_trt)

  # Generate frailty
  # frailty_i ~ Gamma(1/k, 1/k) so mean=1, var=k.
  # Rate(t) = base_rate(season(t)) * frailty_i
  dt_subjects <- data.table(
    id = seq_len(n),
    treatment = assigned_trt,
    enroll_time = enroll_times
  )

  dt_subjects[, dispersion := dispersions[treatment]]
  dt_subjects[, frailty := 1.0]

  # Apply gamma frailty
  has_disp <- dt_subjects$dispersion > 0
  if (any(has_disp)) {
    k <- dt_subjects$dispersion[has_disp]
    dt_subjects$frailty[has_disp] <- rgamma(sum(has_disp), shape = 1 / k, scale = k)
  }

  # Helper: Dropout time
  # Reuse simplified logic or same as nb_sim
  dropout_dt <- if (!is.null(dropout_rate)) data.table(dropout_rate) else NULL

  compute_dropout_time <- function(trt) {
    if (is.null(dropout_dt) || nrow(dropout_dt) == 0) {
      return(Inf)
    }
    dr_sub <- dropout_dt[is.na(treatment) | treatment == trt]
    if (nrow(dr_sub) == 0) dr_sub <- dropout_dt[is.na(treatment)]
    if (nrow(dr_sub) == 0) {
      return(Inf)
    }

    # Piecewise exponential
    t_curr <- 0
    for (j in seq_len(nrow(dr_sub))) {
      rate_j <- dr_sub$rate[j]
      dur_j <- dr_sub$duration[j]
      if (!is.na(rate_j) && rate_j > 0) {
        e <- rexp(1, rate_j)
        if (e < dur_j) {
          return(t_curr + e)
        }
      }
      t_curr <- t_curr + dur_j
    }
    last_rate <- tail(dr_sub$rate, 1)
    if (!is.na(last_rate) && last_rate > 0) {
      return(t_curr + rexp(1, last_rate))
    }
    return(Inf)
  }

  # Helper: Get season for a date
  # Northern hemisphere standard
  get_season <- function(dates) {
    m <- as.integer(format(dates, "%m"))
    d <- as.integer(format(dates, "%d"))

    # Simple fixed dates for seasons
    # Spring: Mar 20 - Jun 20
    # Summer: Jun 21 - Sep 21
    # Fall: Sep 22 - Dec 20
    # Winter: Dec 21 - Mar 19

    # Or simpler:
    # Winter: Dec, Jan, Feb
    # Spring: Mar, Apr, May
    # ...
    # User said: "split at each season". This implies exact boundaries.
    # Let's use meteorological seasons (1st of month) for simplicity and stability,
    # unless exact astronomical is needed. Meteorological is standard in stats usually.
    # Winter: Dec 1 to Feb 28/29
    # Spring: Mar 1 to May 31
    # Summer: Jun 1 to Aug 31
    # Fall: Sep 1 to Nov 30

    m <- as.integer(format(dates, "%m"))
    s <- rep("Winter", length(dates))
    s[m %in% 3:5] <- "Spring"
    s[m %in% 6:8] <- "Summer"
    s[m %in% 9:11] <- "Fall"

    s
  }

  # Helper: Get next season start date
  get_next_season_start <- function(date) {
    y <- as.integer(format(date, "%Y"))
    m <- as.integer(format(date, "%m"))

    candidates <- as.Date(paste0(
      c(y, y, y, y, y + 1),
      c("-03-01", "-06-01", "-09-01", "-12-01", "-03-01")
    ))
    candidates <- candidates[candidates > date]
    min(candidates)
  }

  simulate_subject_seasonal <- function(id, treatment, enroll_time, frailty) {
    d_time <- compute_dropout_time(treatment)
    stop_time <- min(d_time, max_followup)

    # Calendar times
    cal_enroll <- start_date + enroll_time * 365.25 # approx days
    cal_stop <- cal_enroll + stop_time * 365.25

    # Generate intervals based on season changes
    curr_date <- cal_enroll
    intervals <- list()

    # Look up rates
    trt_rates <- fail_rate[fail_rate$treatment == treatment, ]
    # Ensure all seasons present or handle missing

    t_curr <- 0 # relative time

    while (t_curr < stop_time) {
      curr_season <- get_season(curr_date)

      # Find end of this season segment
      next_season_date <- get_next_season_start(curr_date)

      # Limit by stop time
      dist_to_season_change <- as.numeric(next_season_date - curr_date) / 365.25 # in years (or units of rate)
      # Assuming rate units match max_followup units.
      # Usually max_followup is months or years. If rates are per year, we use years.
      # nb_sim usually assumes consistent units.
      # Let's assume input units are consistent (e.g. years).

      # If max_followup is 1 (year), and rate is per year.
      # If max_followup is 12 (months), rate is per month.
      # We need to know the unit of 'rate' vs calendar days.
      # User didn't specify units. "max_followup 1 year".
      # Let's assume time unit is YEAR for rates and max_followup.
      # Conversion from Date diff to Time unit: Days / 365.25.

      time_to_season_end <- as.numeric(next_season_date - curr_date) / 365.25
      remaining_followup <- stop_time - t_curr

      segment_dur <- min(time_to_season_end, remaining_followup)

      # Get rate for current season
      r <- trt_rates$rate[trt_rates$season == curr_season]
      if (length(r) == 0) r <- 0 # Should warn?

      lambda <- r * frailty

      # Simulate events in this segment
      # Poisson process with rate lambda
      # We record events relative to start of segment
      segment_events <- numeric(0)
      if (lambda > 0 && segment_dur > 0) {
        t_evt <- 0
        repeat {
          gap <- rexp(1, lambda)
          t_evt <- t_evt + gap
          if (t_evt <= segment_dur) {
            segment_events <- c(segment_events, t_evt)
          } else {
            break
          }
        }
      }

      # Create rows
      # Start of segment (relative)
      seg_start_rel <- t_curr

      if (length(segment_events) > 0) {
        # Events
        # We create a row for each event
        # Format: start, end, event=1
        # start is previous event time or segment start
        # Wait, standard counting process format: (t_start, t_stop, status)
        # Interval 1: (seg_start, seg_start + evt1, 1)
        # Interval 2: (seg_start + evt1, seg_start + evt2, 1)
        # ...
        # Final: (seg_start + last_evt, seg_start + segment_dur, 0)

        last_t <- 0
        for (evt_t in segment_events) {
          intervals[[length(intervals) + 1]] <- list(
            season = curr_season,
            start = seg_start_rel + last_t,
            end = seg_start_rel + evt_t,
            event = 1,
            cal_start = curr_date + last_t * 365.25,
            cal_end = curr_date + evt_t * 365.25
          )
          last_t <- evt_t
        }
        # Remainder of segment
        if (last_t < segment_dur) {
          intervals[[length(intervals) + 1]] <- list(
            season = curr_season,
            start = seg_start_rel + last_t,
            end = seg_start_rel + segment_dur,
            event = 0,
            cal_start = curr_date + last_t * 365.25,
            cal_end = curr_date + segment_dur * 365.25
          )
        }
      } else {
        # No events, just the segment
        intervals[[length(intervals) + 1]] <- list(
          season = curr_season,
          start = seg_start_rel,
          end = seg_start_rel + segment_dur,
          event = 0,
          cal_start = curr_date,
          cal_end = curr_date + segment_dur * 365.25
        )
      }

      # Advance
      t_curr <- t_curr + segment_dur
      curr_date <- curr_date + segment_dur * 365.25

      # Precision fix to align dates
      if (abs(t_curr - stop_time) < 1e-6) break
    }

    # Combine intervals
    rbindlist(intervals)
  }

  results <- dt_subjects[, simulate_subject_seasonal(id, treatment, enroll_time, frailty), by = .(id, treatment, enroll_time)]

  setDF(results)
  class(results) <- c("nb_sim_seasonal", "data.frame")
  results
}
