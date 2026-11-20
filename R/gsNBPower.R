#' Asymptotic power for a negative binomial group sequential design
#'
#' Evaluates a fixed-size design under alternative rate, dispersion, enrollment,
#' dropout, and follow-up assumptions without resizing the trial.
#'
#' @param design A `gsNB` design, or `NULL` for a fixed-size specification.
#' @param lambda1,lambda2 Control and experimental rates; omitted values inherit
#'   the reference design.
#' @param dispersion NB dispersion, scalar or arm-specific vector.
#' @param accrual_rate,accrual_duration Piecewise enrollment rates and durations.
#' @param trial_duration Final calendar date when no design is supplied.
#' @param dropout_rate Dropout hazards, as in [sample_size_nbinom()].
#' @param max_followup Maximum patient follow-up.
#' @param event_gap Minimum gap between counted events.
#' @param ratio Experimental-to-control allocation ratio.
#' @param rr0 Null rate ratio.
#' @param alpha One-sided efficacy error budget.
#' @param test_type `"score"` or `"wald"`; defaults to the fixed design's test.
#' @param planned_calendar Calendar analysis dates; defaults to `design$T`.
#'   Use `NA` for an interim determined only by other criteria.
#' @param target_info,target_events Expected information or event targets at each
#'   look. Scalars recycle; `NA` omits a requirement.
#' @param min_n,min_followup Minimum enrollment followed by a minimum time since
#'   randomization. These counts do not exclude patients who subsequently drop out.
#' @param max_calendar Hard calendar caps; defaults to the reference final date.
#' @param boundary_update Recompute spending boundaries or retain numerical bounds.
#' @param spending `"min_planned_actual"`, `"planned_final"`, `"expected_final"`,
#'   or `"calendar"`.
#' @param full_spending_at_final Spend the full efficacy budget at the final look.
#' @param usTime,lsTime Optional upper and lower spending-time overrides.
#' @param ... Boundary settings passed to [gsDesign::gsDesign()] when `design`
#'   is `NULL`.
#'
#' @details
#' Power uses a canonical normal approximation with signed log-rate-ratio drift.
#' The score approximation uses separate null and alternative variances and
#' rescales each numerical Z boundary before canonical integration. These are
#' planning approximations, not guarantees for estimated information or random
#' analysis dates. Expected resources are measured at the stopping decision.
#'
#' @return A `gsNBPower` list with `power`, `futility_probability`,
#'   `harm_probability`, `look_results`, expected resources, and `assumptions`.
#' @export
#' @examples
#' fixed <- sample_size_nbinom(
#'   lambda1 = 0.5, lambda2 = 0.35, dispersion = 0.5, power = 0.9,
#'   accrual_rate = 10, accrual_duration = 12, trial_duration = 24
#' )
#' design <- gsDesignNB::toInteger(gsNB(fixed))
#' gsNBPower(design, lambda2 = 0.4)
gsNBPower <- function(
  design = NULL, lambda1 = NULL, lambda2 = NULL, dispersion = NULL,
  accrual_rate = NULL, accrual_duration = NULL, trial_duration = NULL,
  dropout_rate = NULL, max_followup = NULL, event_gap = NULL,
  ratio = NULL, rr0 = NULL, alpha = NULL, test_type = NULL,
  planned_calendar = NULL, target_info = NULL, target_events = NULL,
  min_n = NULL, min_followup = NULL, max_calendar = NULL,
  boundary_update = c("recompute", "fixed"),
  spending = c("min_planned_actual", "planned_final", "expected_final", "calendar"),
  full_spending_at_final = TRUE, usTime = NULL, lsTime = NULL, ...
) {
  boundary_update <- match.arg(boundary_update)
  spending <- match.arg(spending)
  if (!is.null(design) && !inherits(design, "gsNB")) {
    stop("design must be a gsNB object.", call. = FALSE)
  }
  inherit <- function(value, fallback) if (is.null(value)) fallback else value
  inputs <- if (is.null(design)) list() else design$nb_design$inputs
  lambda1 <- inherit(lambda1, inputs$lambda1)
  lambda2 <- inherit(lambda2, inputs$lambda2)
  dispersion <- inherit(dispersion, inputs$dispersion)
  ratio <- inherit(ratio, inherit(inputs$ratio, 1))
  rr0 <- inherit(rr0, inherit(inputs$rr0, 1))
  dropout_rate <- inherit(dropout_rate, inherit(inputs$dropout_rate, 0))
  max_followup <- inherit(max_followup, inherit(inputs$max_followup, Inf))
  event_gap <- inherit(event_gap, inherit(inputs$event_gap, 0))
  alpha <- inherit(alpha, if (is.null(design)) 0.025 else design$alpha)
  test_type <- match.arg(inherit(test_type, inherit(inputs$test_type, "wald")),
                         c("score", "wald"))
  if (is.null(lambda1) || is.null(lambda2) || is.null(dispersion)) {
    stop("lambda1, lambda2, and dispersion are required without a design.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) stop("alpha must be in (0, 1).", call. = FALSE)
  if (!is.logical(full_spending_at_final) || length(full_spending_at_final) != 1L ||
      is.na(full_spending_at_final)) {
    stop("full_spending_at_final must be TRUE or FALSE.", call. = FALSE)
  }
  planned_calendar <- inherit(planned_calendar, if (is.null(design)) trial_duration else design$T)
  if (!is.numeric(planned_calendar) || !length(planned_calendar) ||
      any(!is.finite(planned_calendar[!is.na(planned_calendar)])) ||
      any(planned_calendar[!is.na(planned_calendar)] <= 0) ||
      any(diff(planned_calendar[!is.na(planned_calendar)]) <= 0)) {
    stop("planned_calendar must contain positive increasing dates or NA interims.", call. = FALSE)
  }
  if (!is.null(design) && length(planned_calendar) != design$k) {
    stop("planned_calendar must have design$k elements.", call. = FALSE)
  }
  if (!is.null(design)) {
    planned_accrual <- .gs_nb_accrual(design)
    accrual_rate <- inherit(accrual_rate, planned_accrual$rate)
    accrual_duration <- inherit(accrual_duration, planned_accrual$duration)
    capped <- .gs_nb_cap_accrual(accrual_rate, accrual_duration,
                                utils::tail(design$n_total, 1))
    accrual_rate <- capped$rate
    accrual_duration <- capped$duration
  }
  if (is.null(accrual_rate) || is.null(accrual_duration)) {
    stop("accrual_rate and accrual_duration are required without a design.", call. = FALSE)
  }
  assumptions <- list(
    lambda1 = lambda1, lambda2 = lambda2, dispersion = dispersion,
    accrual_rate = accrual_rate, accrual_duration = accrual_duration,
    dropout_rate = dropout_rate, max_followup = max_followup, event_gap = event_gap,
    ratio = ratio, rr0 = rr0, alpha = alpha, test_type = test_type
  )
  schedule <- .gs_nb_schedule(
    planned_calendar, target_info, target_events, min_n, min_followup,
    max_calendar, assumptions,
    if (is.null(design)) trial_duration else utils::tail(design$T, 1)
  )
  look_ids <- which(schedule$status == "analyzed")
  analysis_dates <- schedule$analysis_time[look_ids]
  summaries <- lapply(analysis_dates, function(time) .gs_nb_expectation(time, assumptions))
  info <- vapply(summaries, function(stage) 1 / stage$variance, numeric(1))
  info_null <- vapply(summaries, function(stage) 1 / stage$variance_null, numeric(1))
  if (any(!is.finite(info)) || any(info <= 0) || any(diff(info) <= 0)) {
    stop("Expected information must be positive, finite, and increasing.", call. = FALSE)
  }
  n_looks <- length(info)
  if (is.null(design)) {
    settings <- list(...)
    settings$k <- n_looks
    settings$alpha <- alpha
    settings$n.I <- info
    settings$maxn.IPlan <- utils::tail(info, 1)
    settings$delta <- abs(log(lambda1 * rr0 / lambda2))
    if (is.null(settings$test.type)) settings$test.type <- 1
    if (n_looks == 1L) {
      design <- list(k = 1L, alpha = alpha, beta = 0.1, n.I = info,
                     test.type = settings$test.type, delta = settings$delta,
                     upper = list(bound = stats::qnorm(1 - alpha)), lower = NULL,
                     tol = 1e-6, r = 18)
    } else {
      design <- do.call(gsDesign::gsDesign, settings)
    }
    design$T <- analysis_dates
    design$delta1 <- log(lambda2 / lambda1)
    design$delta0 <- log(rr0)
    if (spending == "min_planned_actual") spending <- "expected_final"
    look_ids_for_bounds <- seq_len(n_looks)
  } else {
    look_ids_for_bounds <- look_ids
  }
  orientation <- if (!is.null(design$delta1) && design$delta1 > 0) -1 else 1
  theta <- orientation * log(lambda1 * rr0 / lambda2)
  bounds <- .gs_nb_bounds(
    design, info, analysis_dates, look_ids_for_bounds, boundary_update,
    spending, full_spending_at_final, usTime, lsTime, alpha = alpha
  )
  probabilities <- .gs_nb_probabilities(
    info, if (test_type == "score") info_null else info,
    theta, bounds$upper, bounds$lower, bounds$harm, design$r
  )
  look_results <- data.frame(
    look = look_ids, analysis_time = analysis_dates,
    n_total = vapply(summaries, function(stage) stage$n_total, numeric(1)),
    information = info, information_null = info_null,
    information_fraction = info / utils::tail(design$n.I, 1),
    expected_events = vapply(summaries, function(stage) stage$total_events, numeric(1)),
    exposure_ctrl = vapply(summaries, function(stage) stage$exposure[1], numeric(1)),
    exposure_exp = vapply(summaries, function(stage) stage$exposure[2], numeric(1)),
    upper_bound = bounds$upper, lower_bound = bounds$lower, harm_bound = bounds$harm,
    spending_time = bounds$us_time, lower_spending_time = bounds$ls_time,
    efficacy_prob = probabilities$upper, futility_prob = probabilities$lower,
    harm_prob = probabilities$harm,
    cumulative_efficacy = cumsum(probabilities$upper),
    cumulative_futility = cumsum(probabilities$lower),
    continuation_prob = probabilities$continuation, status = "analyzed"
  )
  decision_mass <- probabilities$upper + probabilities$lower + probabilities$harm
  decision_mass[n_looks] <- decision_mass[n_looks] + utils::tail(probabilities$continuation, 1)
  expected_n <- sum(decision_mass * look_results$n_total)
  expected_events <- sum(decision_mass * look_results$expected_events)
  expected_analysis_time <- sum(decision_mass * analysis_dates)
  if (length(look_ids) < nrow(schedule)) {
    all_looks <- look_results[rep(NA_integer_, nrow(schedule)), ]
    all_looks$look <- seq_len(nrow(schedule))
    all_looks$status <- schedule$status
    all_looks$efficacy_prob <- 0
    all_looks$futility_prob <- 0
    all_looks$harm_prob <- 0
    all_looks[look_ids, ] <- look_results
    look_results <- all_looks
  }
  look_results$target_met <- schedule$target_met
  rownames(look_results) <- NULL
  structure(list(
    power = sum(probabilities$upper), futility_probability = sum(probabilities$lower),
    harm_probability = sum(probabilities$harm), look_results = look_results,
    expected_n = expected_n, expected_events = expected_events,
    expected_analysis_time = expected_analysis_time,
    assumptions = c(assumptions, list(boundary_update = boundary_update, spending = spending,
                    full_spending_at_final = full_spending_at_final,
                    planned_calendar = planned_calendar,
                    target_info = target_info, target_events = target_events,
                    min_n = min_n, min_followup = min_followup,
                    max_calendar = max_calendar)),
    reference_design = design
  ), class = "gsNBPower")
}

.gs_nb_expectation <- function(time, assumptions) {
  do.call(sample_size_nbinom, c(assumptions, list(power = NULL, trial_duration = time)))
}

.gs_nb_schedule <- function(calendar, target_info, target_events, min_n,
                            min_followup, max_calendar, assumptions, final_date) {
  n_looks <- length(calendar)
  expand <- function(value, name, fallback = NA_real_) {
    if (is.null(value)) value <- fallback
    if (!is.numeric(value) || !length(value) ||
        !length(value) %in% c(1L, n_looks) ||
        any(!is.finite(value[!is.na(value)])) || any(value[!is.na(value)] < 0)) {
      stop(name, " must have one or k nonnegative finite values or NA.", call. = FALSE)
    }
    rep(value, length.out = n_looks)
  }
  if (is.null(final_date) || !is.finite(final_date)) final_date <- utils::tail(calendar, 1)
  if (length(final_date) != 1L || !is.finite(final_date) || final_date <= 0) {
    stop("A positive finite final calendar date is required.", call. = FALSE)
  }
  info_targets <- expand(target_info, "target_info")
  event_targets <- expand(target_events, "target_events")
  enrollment_targets <- expand(min_n, "min_n")
  followup <- expand(min_followup, "min_followup", 0)
  if (any(!is.na(followup) & followup > 0 & is.na(enrollment_targets))) {
    stop("min_followup requires min_n at the same look.", call. = FALSE)
  }
  caps <- expand(max_calendar, "max_calendar", final_date)
  caps[is.na(caps)] <- final_date
  caps <- pmin(caps, utils::tail(caps, 1))
  dates <- rep(NA_real_, n_looks)
  status <- rep("skipped", n_looks)
  target_met <- rep(FALSE, n_looks)
  last_date <- 0
  value_at <- function(time, metric) {
    if (time <= 0) return(0)
    summary <- .gs_nb_expectation(time, assumptions)
    if (metric == "information") 1 / summary$variance else summary$total_events
  }
  for (look in seq_len(n_looks)) {
    cap <- caps[look]
    floor <- if (is.na(calendar[look])) last_date else max(last_date, calendar[look])
    if (!is.na(enrollment_targets[look])) {
      enrolled_time <- .gs_nb_enrollment_time(
        enrollment_targets[look], assumptions$accrual_rate, assumptions$accrual_duration
      )
      floor <- max(floor, enrolled_time + if (is.na(followup[look])) 0 else followup[look])
    }
    candidate <- floor
    for (metric in c("information", "events")) {
      target <- if (metric == "information") info_targets[look] else event_targets[look]
      if (is.na(target)) next
      if (value_at(cap, metric) < target) {
        candidate <- Inf
      } else if (is.finite(candidate) && value_at(candidate, metric) < target) {
        candidate <- stats::uniroot(
          function(time) value_at(time, metric) - target,
          interval = c(candidate, cap), tol = 1e-7
        )$root
      }
    }
    target_met[look] <- is.finite(candidate) && candidate <= cap
    if (look == n_looks) {
      dates[look] <- cap
      status[look] <- "analyzed"
    } else if (target_met[look] && candidate > last_date && candidate < caps[n_looks]) {
      dates[look] <- candidate
      status[look] <- "analyzed"
      last_date <- candidate
    }
  }
  data.frame(analysis_time = dates, status = status, target_met = target_met)
}

.gs_nb_enrollment_time <- function(target, rate, duration) {
  if (target <= 0) return(0)
  cumulative <- cumsum(rate * duration)
  period <- which(cumulative >= target)[1]
  if (is.na(period)) return(Inf)
  prior_n <- if (period == 1L) 0 else cumulative[period - 1L]
  prior_time <- if (period == 1L) 0 else sum(duration[seq_len(period - 1L)])
  prior_time + (target - prior_n) / rate[period]
}

.gs_nb_accrual <- function(design) {
  nb <- design$nb_design
  base <- .gs_nb_expectation(utils::tail(design$T, 1), list(
    lambda1 = nb$inputs$lambda1, lambda2 = nb$inputs$lambda2,
    dispersion = nb$inputs$dispersion, accrual_rate = nb$accrual_rate,
    accrual_duration = nb$accrual_duration, dropout_rate = nb$inputs$dropout_rate,
    max_followup = nb$inputs$max_followup, event_gap = nb$inputs$event_gap,
    ratio = nb$inputs$ratio
  ))
  list(rate = nb$accrual_rate * utils::tail(design$n_total, 1) / base$n_total,
       duration = nb$accrual_duration)
}

.gs_nb_cap_accrual <- function(rate, duration, cap) {
  if (!is.numeric(rate) || !is.numeric(duration) || length(rate) != length(duration) ||
      !length(rate) || any(!is.finite(rate)) || any(rate <= 0) ||
      any(!is.finite(duration)) || any(duration <= 0)) {
    stop("Enrollment rates and durations must be positive finite equal-length vectors.", call. = FALSE)
  }
  totals <- rate * duration
  endpoint <- which(cumsum(totals) >= cap)[1]
  if (!is.na(endpoint)) {
    prior <- if (endpoint == 1L) 0 else sum(totals[seq_len(endpoint - 1L)])
    duration[endpoint] <- (cap - prior) / rate[endpoint]
    rate <- rate[seq_len(endpoint)]
    duration <- duration[seq_len(endpoint)]
  }
  list(rate = rate, duration = duration)
}

.gs_nb_spending <- function(info, planned_info, dates, final_info, strategy,
                            override = NULL, full_final = TRUE, history = NULL) {
  time <- if (!is.null(override)) override else switch(
    strategy,
    min_planned_actual = pmin(info, planned_info) / utils::tail(planned_info, 1),
    planned_final = info / utils::tail(planned_info, 1),
    expected_final = info / final_info,
    calendar = dates / utils::tail(dates, 1)
  )
  if (length(time) != length(info) || any(!is.finite(time)) || any(time <= 0)) {
    stop("Spending times must be positive finite values, one per active look.", call. = FALSE)
  }
  time <- pmin(time, 1)
  if (length(time) > 1L) time[-length(time)] <- pmin(time[-length(time)], 0.999)
  if (length(history)) {
    time[seq_along(history)] <- history
    future <- seq.int(length(history) + 1L, length(time))
    time[future] <- pmax(time[future], utils::tail(history, 1))
  }
  if (full_final) time[length(time)] <- 1
  if (any(diff(time) < 0)) stop("Spending times must be nondecreasing.", call. = FALSE)
  time
}

.gs_nb_bounds <- function(design, info, dates, look_ids, policy, spending,
                          full_final, usTime = NULL, lsTime = NULL,
                          final_info = utils::tail(info, 1), alpha = design$alpha,
                          history = NULL) {
  planned_info <- design$n.I[look_ids]
  upper_time <- if (is.null(usTime)) design$usTime else usTime
  lower_time <- if (is.null(lsTime)) design$lsTime else lsTime
  if (!is.null(upper_time) && length(upper_time) == design$k) upper_time <- upper_time[look_ids]
  if (!is.null(lower_time) && length(lower_time) == design$k) lower_time <- lower_time[look_ids]
  planned_info[length(planned_info)] <- utils::tail(design$n.I, 1)
  us_time <- .gs_nb_spending(info, planned_info, dates, final_info, spending,
                            upper_time, full_final, history$upper)
  ls_time <- .gs_nb_spending(info, planned_info, dates, final_info, spending,
                            lower_time, full_final, history$lower)
  if (policy == "fixed") {
    upper <- design$upper$bound[look_ids]
    lower <- if (is.null(design$lower)) rep(-Inf, length(info)) else design$lower$bound[look_ids]
    harm <- if (is.null(design$harm)) rep(-Inf, length(info)) else design$harm$bound[look_ids]
  } else if (length(info) == 1L) {
    budget <- if (is.null(design$upper$sf)) alpha else
      .gs_nb_spending_function(design$upper$sf)(alpha, us_time, design$upper$param)$spend
    upper <- stats::qnorm(1 - utils::tail(budget, 1))
    lower <- if (design$test.type == 1L) -Inf else if (design$test.type == 2L) -upper else upper
    harm <- -Inf
  } else {
    settings <- list(
      k = length(info), test.type = design$test.type, alpha = alpha,
      beta = design$beta, astar = design$astar, delta = design$delta,
      n.I = info, maxn.IPlan = utils::tail(design$n.I, 1),
      sfu = .gs_nb_spending_function(design$upper$sf), sfupar = design$upper$param,
      tol = design$tol, r = design$r, usTime = us_time, lsTime = ls_time
    )
    if (!is.null(design$lower)) {
      settings$sfl <- .gs_nb_spending_function(design$lower$sf)
      settings$sflpar <- design$lower$param
    }
    if (!is.null(design$harm)) {
      settings$sfharm <- .gs_nb_spending_function(design$harm$sf)
      settings$sfharmparam <- design$harm$param
    }
    for (mask in c("testUpper", "testLower", "testHarm")) {
      if (!is.null(design[[mask]])) settings[[mask]] <- rep(design[[mask]], length.out = design$k)[look_ids]
    }
    updated <- do.call(gsDesign::gsDesign, settings)
    upper <- updated$upper$bound
    lower <- if (is.null(updated$lower)) rep(-Inf, length(info)) else updated$lower$bound
    harm <- if (is.null(updated$harm)) rep(-Inf, length(info)) else updated$harm$bound
  }
  list(upper = upper, lower = lower, harm = harm, us_time = us_time, ls_time = ls_time)
}

.gs_nb_spending_function <- function(fun) {
  if (is.function(fun) && "checkScalar" %in% all.names(body(fun)) &&
      !exists("checkScalar", envir = environment(fun), inherits = TRUE)) {
    environment(fun) <- asNamespace("gsDesign")
  }
  fun
}

.gs_nb_probabilities <- function(info, info_null, theta, upper, lower, harm, r = 18) {
  scaling <- sqrt(info / info_null)
  lower_operational <- pmax(lower, harm)
  integrate <- function(n_looks, lower_last = NULL) {
    active <- seq_len(n_looks)
    lower_scaled <- lower_operational[active] * scaling[active]
    if (!is.null(lower_last)) lower_scaled[n_looks] <- lower_last * scaling[n_looks]
    upper_scaled <- upper[active] * scaling[active]
    if (n_looks == 1L) {
      mean_z <- theta * sqrt(info[1])
      return(list(upper = stats::pnorm(upper_scaled - mean_z, lower.tail = FALSE),
                  lower = stats::pnorm(lower_scaled - mean_z)))
    }
    probability <- gsDesign::gsProbability(
      k = n_looks, theta = theta, n.I = info[active],
      a = lower_scaled, b = upper_scaled, r = r
    )
    list(upper = as.numeric(probability$upper$prob), lower = as.numeric(probability$lower$prob))
  }
  probability <- integrate(length(info))
  harm_prob <- vapply(seq_along(info), function(look) {
    if (!is.finite(harm[look])) 0 else utils::tail(integrate(look, harm[look])$lower, 1)
  }, numeric(1))
  list(upper = probability$upper, lower = pmax(0, probability$lower - harm_prob),
       harm = harm_prob, continuation = pmax(0, 1 - cumsum(probability$upper + probability$lower)))
}

#' @param x A `gsNBPower` result.
#' @rdname gsNBPower
#' @export
print.gsNBPower <- function(x, ...) {
  cat("Asymptotic NB group sequential power:", formatC(x$power, digits = 4, format = "f"), "\n")
  print(x$look_results, row.names = FALSE)
  invisible(x)
}