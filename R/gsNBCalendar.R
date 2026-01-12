#' Group sequential design for negative binomial outcomes
#'
#' Creates a group sequential design for negative binomial outcomes based on
#' sample size calculations from [sample_size_nbinom()].
#'
#' @param x An object of class `sample_size_nbinom_result` from
#'   [sample_size_nbinom()].
#' @param k Number of analyses (interim + final). Default is 3.
#' @param test.type Test type as in [gsDesign::gsDesign()]:
#'   \describe{
#'     \item{1}{One-sided}
#'     \item{2}{Two-sided symmetric}
#'     \item{3}{Two-sided, asymmetric, binding futility bound, beta-spending}
#'     \item{4}{Two-sided, asymmetric, non-binding futility bound, beta-spending}
#'     \item{5}{Two-sided, asymmetric, binding futility bound, lower spending}
#'     \item{6}{Two-sided, asymmetric, non-binding futility bound, lower spending}
#'   }
#'   Default is 4.
#' @param alpha Type I error (one-sided). Default is 0.025.
#' @param beta Type II error (1 - power). Default is 0.1.
#' @param astar Allocated Type I error for lower bound for test.type = 5 or 6.
#'   Default is 0.
#' @param delta Standardized effect size. Default is 0 (computed from design).
#' @param sfu Spending function for upper bound. Default is `gsDesign::sfHSD`.
#' @param sfupar Parameter for upper spending function. Default is -4.
#' @param sfl Spending function for lower bound. Default is `gsDesign::sfHSD`.
#' @param sflpar Parameter for lower spending function. Default is -2.
#' @param tol Tolerance for convergence. Default is 1e-06.
#' @param r Integer controlling grid size for numerical integration.
#'   Default is 18.
#' @param usTime Spending time for upper bound (optional).
#' @param lsTime Spending time for lower bound (optional).
#' @param analysis_times Vector of calendar times for each analysis.
#'   Must have length k. These times are stored in the `T`
#'   element and displayed by [gsDesign::gsBoundSummary()].
#'
#' @return An object of class `gsNB` which inherits from `gsDesign`
#'   and `sample_size_nbinom_result`. 
#'   While the final sample size would be planned total enrollment, interim analysis
#'   sample sizes are the expected number enrolled at the times specified in `analysis_times`.
#' Output value contains all elements from
#'   [gsDesign::gsDesign()] plus:
#'   \describe{
#'     \item{nb_design}{The original `sample_size_nbinom_result` object}
#'     \item{n1}{A vector with sample size per analysis for group 1}
#'     \item{n2}{A vector with sample size per analysis for group 2}
#'     \item{n_total}{A vector with total sample size per analysis}
#'     \item{events}{A vector with expected total events per analysis}
#'     \item{events1}{A vector with expected events per analysis for group 1}
#'     \item{events2}{A vector with expected events per analysis for group 2}
#'     \item{T}{Calendar time at each analysis (if `analysis_times` provided)}
#'   }
#'   Note that `n.I` in the returned object represents the statistical information
#'   at each analysis.
#'
#' @importFrom gsDesign gsDesign sfHSD
#'
#' @export
#'
#' @references
#' Jennison, C. and Turnbull, B.W. (2000),
#' _Group Sequential Methods with Applications to Clinical Trials_.
#' Boca Raton: Chapman and Hall.
#'
#' @examples
#' # First create a sample size calculation
#' nb_ss <- sample_size_nbinom(
#'   lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
#'   accrual_rate = 10, accrual_duration = 20, trial_duration = 24
#' )
#'
#' # Then create a group sequential design with analysis times
#' gs_design <- gsNBCalendar(nb_ss,
#'   k = 3, test.type = 4,
#'   analysis_times = c(10, 18, 24)
#' )
gsNBCalendar <- function(
  x,
  k = 3,
  test.type = 4,
  alpha = 0.025,
  beta = 0.1,
  astar = 0,
  delta = 0,
  sfu = gsDesign::sfHSD,
  sfupar = -4,
  sfl = gsDesign::sfHSD,
  sflpar = -2,
  tol = 1e-06,
  r = 18,
  usTime = NULL,
  lsTime = NULL,
  analysis_times = NULL
) {
  # Validate input

  if (!inherits(x, "sample_size_nbinom_result")) {
    stop("x must be an object of class 'sample_size_nbinom_result'")
  }

  if (is.null(analysis_times)) {
    stop("analysis_times must be provided")
  }

  if (length(analysis_times) != k) {
    stop("analysis_times must have length k (", k, ")")
  }

  # For NB, the test statistic is based on log(lambda2/lambda1)
  # RR = lambda2/lambda1, so delta1 = log(RR)
  # When treatment is beneficial (lambda2 < lambda1), RR < 1 and delta1 < 0
  # gsBoundSummary with logdelta=TRUE will show exp(delta1) = RR
  risk_ratio <- x$inputs$lambda2 / x$inputs$lambda1
  delta1 <- log(risk_ratio)

  # Statistical information is the inverse of the variance of the log rate ratio
  # This is used as n.fix so that gs$n.I represents actual statistical information
  info_fixed <- 1 / x$variance

  # Call gsDesign with the provided parameters
  gs <- gsDesign::gsDesign(
    k = k,
    test.type = test.type,
    alpha = alpha,
    beta = beta,
    astar = astar,
    delta = delta,
    delta1 = delta1,
    delta0 = 0,
    n.fix = info_fixed,
    sfu = sfu,
    sfupar = sfupar,
    sfl = sfl,
    sflpar = sflpar,
    tol = tol,
    r = r,
    usTime = usTime,
    lsTime = lsTime
  )

  # Calculate sample sizes per analysis based on information fraction
  # gs$n.I contains the statistical information at each analysis
  # Scale to get sample sizes using the relationship: info / info_fixed = n / n_fixed
  
  # Inflate accrual rate to match the max information required by the group sequential design
  inflation_factor <- utils::tail(gs$n.I, 1) / info_fixed
  # Use the computed accrual rate from the fixed design (which might be scaled)
  # and inflate it further.
  # Note: x$accrual_rate is the rate used in the fixed design calculation.
  base_accrual_rate <- x$accrual_rate
  new_accrual_rate <- base_accrual_rate * inflation_factor
  
  # Calculate expected enrollment and information at each analysis time
  n_total_at_analysis <- numeric(k)
  info_at_analysis <- numeric(k)
  events_at_analysis <- numeric(k)
  events1_at_analysis <- numeric(k)
  events2_at_analysis <- numeric(k)
  
  for (i in 1:k) {
    # We use sample_size_nbinom to calculate the properties at time t
    # We pass power = NULL to force calculation based on accrual/duration
    # sample_size_nbinom will handle truncation of accrual if analysis_times[i] < total accrual duration
    res_i <- sample_size_nbinom(
      lambda1 = x$inputs$lambda1,
      lambda2 = x$inputs$lambda2,
      rr0 = x$inputs$rr0,
      dispersion = x$inputs$dispersion,
      power = NULL, # Calculate based on fixed N/duration
      alpha = x$inputs$alpha,
      sided = x$inputs$sided,
      ratio = x$inputs$ratio,
      accrual_rate = new_accrual_rate,
      accrual_duration = x$inputs$accrual_duration,
      trial_duration = analysis_times[i],
      dropout_rate = x$inputs$dropout_rate,
      max_followup = x$inputs$max_followup,
      event_gap = x$inputs$event_gap
    )
    
    n_total_at_analysis[i] <- res_i$n_total
    info_at_analysis[i] <- 1 / res_i$variance
    events_at_analysis[i] <- res_i$total_events
    events1_at_analysis[i] <- res_i$events_n1
    events2_at_analysis[i] <- res_i$events_n2
  }

  # Update the gsDesign object with the actual information and timing
  gs$n.I <- info_at_analysis
  gs$timing <- info_at_analysis / utils::tail(info_at_analysis, 1)
  
  # Per-group sample sizes
  ratio <- x$inputs$ratio
  n1_at_analysis <- n_total_at_analysis / (1 + ratio)
  n2_at_analysis <- n_total_at_analysis * ratio / (1 + ratio)

  # Build result object

  # Start with the gsDesign object
  result <- gs

  # Add negative binomial specific components
  result$nb_design <- x
  result$n1 <- n1_at_analysis
  result$n2 <- n2_at_analysis
  result$n_total <- n_total_at_analysis
  result$events <- events_at_analysis
  result$events1 <- events1_at_analysis
  result$events2 <- events2_at_analysis
  result$ratio <- ratio
  result$usTime <- usTime
  result$lsTime <- lsTime
  result$T <- analysis_times

  # Set the class to inherit from both gsDesign and sample_size_nbinom_result
  class(result) <- c("gsNB", "gsDesign", "sample_size_nbinom_result")

  return(result)
}


#' Compute statistical information at analysis time
#'
#' Computes the statistical information for the log rate ratio at a given
#' analysis time, accounting for staggered enrollment and varying exposure times.
#'
#' @param analysis_time The calendar time of the analysis.
#' @param accrual_rate The enrollment rate (subjects per time unit).
#' @param accrual_duration The duration of the enrollment period.
#' @param lambda1 Event rate for group 1 (control).
#' @param lambda2 Event rate for group 2 (treatment).
#' @param dispersion The negative binomial dispersion parameter.
#' @param ratio Allocation ratio (n2/n1). Default is 1.
#' @param dropout_rate Dropout rate (hazard rate). Default is 0.
#' @param event_gap Gap duration after each event during which no new events are counted.
#'   Default is 0.
#' @param max_followup Maximum follow-up time per subject. Exposure time is
#'   truncated at this value. Default is `Inf` (no truncation).
#'
#' @return The statistical information (inverse of variance) at the analysis time.
#'
#' @export
#'
#' @examples
#' compute_info_at_time(
#'   analysis_time = 12,
#'   accrual_rate = 10,
#'   accrual_duration = 10,
#'   lambda1 = 0.5,
#'   lambda2 = 0.3,
#'   dispersion = 0.1
#' )
compute_info_at_time <- function(
  analysis_time, accrual_rate, accrual_duration,
  lambda1, lambda2, dispersion, ratio = 1,
  dropout_rate = 0, event_gap = 0, max_followup = Inf
) {
  # Number of subjects enrolled by analysis_time
  enrollment_time <- min(analysis_time, accrual_duration)
  n_total <- accrual_rate * enrollment_time
  n1 <- n_total / (1 + ratio)
  n2 <- n_total * ratio / (1 + ratio)

  # Average exposure calculation must account for dropout
  # We integrate dropout probability over the enrollment period and follow-up

  # Helper for expected exposure for a subject enrolled at time t
  # Exposure is min(analysis_time - t, dropout_time, max_followup)
  # If analysis_time - t > 0

  # Simplified average exposure calculation with dropout
  # 1. Without dropout:
  #    If analysis_time <= accrual_duration: avg = analysis_time / 2
  #    If analysis_time > accrual_duration: avg = analysis_time - accrual_duration/2

  # 2. With dropout, we need to integrate:
  #    E[min(T-t, E)] where E ~ Exp(dropout_rate)
  #    This is (1 - exp(-lambda * (T-t))) / lambda

  # We need the average of this over t in [0, enrollment_time]
  # Let tau = T - t. As t goes 0 -> enrollment_time, tau goes T -> T-enrollment_time
  # Range of follow-up times: [min_f, max_f]
  max_f <- analysis_time
  min_f <- analysis_time - enrollment_time

  # Cap follow-up at max_followup
  # We need to calculate average exposure over [min_f, max_f]
  # but exposure is capped at max_followup

  calc_avg_exposure <- function(min_t, max_t, rate, max_fu) {
    if (min_t >= max_t) {
      return(0)
    }

    # If entire range is above max_fu, exposure is constant (capped)
    if (min_t >= max_fu) {
      if (rate <= 0) {
        return(max_fu)
      }
      return((1 - exp(-rate * max_fu)) / rate)
    }

    # If entire range is below max_fu, standard calculation
    if (max_t <= max_fu) {
      if (rate <= 0) {
        return((min_t + max_t) / 2)
      } else {
        term <- (exp(-rate * min_t) - exp(-rate * max_t)) / (rate * (max_t - min_t))
        return((1 - term) / rate)
      }
    }

    # Split range: [min_t, max_fu] and [max_fu, max_t]
    w1 <- (max_fu - min_t) / (max_t - min_t)
    w2 <- (max_t - max_fu) / (max_t - min_t)

    e1 <- calc_avg_exposure(min_t, max_fu, rate, max_fu)
    e2 <- calc_avg_exposure(max_fu, max_t, rate, max_fu) # This will hit the first case

    return(w1 * e1 + w2 * e2)
  }

  # Handle vector parameters
  if (length(dropout_rate) == 1) dropout_rate <- rep(dropout_rate, 2)
  if (length(max_followup) == 1) max_followup <- rep(max_followup, 2)
  if (length(dispersion) == 1) dispersion <- rep(dispersion, 2)

  avg_exposure1 <- calc_avg_exposure(min_f, max_f, dropout_rate[1], max_followup[1])
  avg_exposure2 <- calc_avg_exposure(min_f, max_f, dropout_rate[2], max_followup[2])

  # Adjust rates for event gap
  if (!is.null(event_gap) && event_gap > 0) {
    lambda1_eff <- lambda1 / (1 + lambda1 * event_gap)
    lambda2_eff <- lambda2 / (1 + lambda2 * event_gap)
  } else {
    lambda1_eff <- lambda1
    lambda2_eff <- lambda2
  }

  # Expected events per subject
  mu1 <- lambda1_eff * avg_exposure1
  mu2 <- lambda2_eff * avg_exposure2

  # Variance of log rate ratio
  # Var(log(lambda2/lambda1)) = (1/mu1 + k)/n1 + (1/mu2 + k)/n2
  k1 <- dispersion[1]
  k2 <- dispersion[2]
  variance <- (1 / mu1 + k1) / n1 + (1 / mu2 + k2) / n2


  # Information is inverse of variance
  info <- 1 / variance

  return(info)
}


#' Summary for gsNB objects
#'
#' Provides a textual summary of a group sequential design for negative binomial
#' outcomes, similar to the summary provided by [gsDesign::gsDesign()].
#' For tabular output, use [gsDesign::gsBoundSummary()] directly on
#' the gsNB object.
#'
#' @param object An object of class `gsNB`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A character string summarizing the design (invisibly). The summary
#'   is also printed to the console.
#'
#' @export
#'
#' @examples
#' nb_ss <- sample_size_nbinom(
#'   lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
#'   accrual_rate = 10, accrual_duration = 20, trial_duration = 24
#' )
#' gs_design <- gsNBCalendar(nb_ss, k = 3)
#' summary(gs_design)
#'
#' # For tabular bounds summary, use gsBoundSummary() directly:
#' gsBoundSummary(gs_design)
summary.gsNB <- function(object, ...) {
  # Extract negative binomial design information
  nb <- object$nb_design
  inputs <- nb$inputs

  # Calculate risk ratio (lambda2/lambda1)
  risk_ratio <- inputs$lambda2 / inputs$lambda1

  # Determine test type description
  test_type_desc <- switch(as.character(object$test.type),
    "1" = "One-sided",
    "2" = "Two-sided symmetric",
    "3" = "Asymmetric two-sided with binding futility bound",
    "4" = "Asymmetric two-sided with non-binding futility bound",
    "5" = "Asymmetric two-sided with binding futility bound (lower spending)",
    "6" = "Asymmetric two-sided with non-binding futility bound (lower spending)",
    "Unknown test type"
  )

  # Format exposure text
  exposure_text <- sprintf("average exposure %.2f", nb$exposure)
  if (!is.null(inputs$event_gap) && inputs$event_gap > 0 && !is.null(nb$exposure_at_risk_n1)) {
    exposure_text <- sprintf(
      "average exposure (calendar) %.2f, (at-risk n1=%.2f, n2=%.2f)",
      nb$exposure, nb$exposure_at_risk_n1, nb$exposure_at_risk_n2
    )
  }

  # Build the summary text
  max_followup_str <- ""
  if (!is.null(inputs$max_followup)) {
    # max_followup is stored as a vector of length 2, but enforced to be scalar
    max_followup_str <- sprintf(", max follow-up %.1f", inputs$max_followup[1])
  }

  event_gap_str <- ""
  if (!is.null(inputs$event_gap) && any(inputs$event_gap > 0)) {
    event_gap_str <- sprintf(", event gap %.2f", inputs$event_gap[1])
  }

  dropout_str <- ""
  if (!is.null(inputs$dropout_rate)) {
    if (length(inputs$dropout_rate) == 1) {
      if (inputs$dropout_rate > 0) {
        dropout_str <- sprintf(", dropout rate %.4f", inputs$dropout_rate)
      }
    } else {
      if (any(inputs$dropout_rate > 0)) {
        dropout_str <- sprintf(", dropout rates (%.4f, %.4f)", inputs$dropout_rate[1], inputs$dropout_rate[2])
      }
    }
  }

  # Spending function descriptions
  upper_spend <- if (!is.null(object$upper$name)) {
    paste0(
      "Upper spending: ", object$upper$name,
      if (!is.null(object$upper$parname)) paste0(" (", object$upper$parname, " = ", object$upper$param, ")") else ""
    )
  } else {
    "Upper spending: Custom"
  }

  lower_spend <- if (!is.null(object$lower$name)) {
    paste0(
      "Lower spending: ", object$lower$name,
      if (!is.null(object$lower$parname)) paste0(" (", object$lower$parname, " = ", object$lower$param, ")") else ""
    )
  } else {
    "Lower spending: Custom"
  }

  summary_text <- sprintf(
    paste0(
      "%s group sequential design for negative binomial outcomes, ",
      "%d analyses, ",
      "total sample size %.1f, ",
      "%.0f percent power, ",
      "%.1f percent (1-sided) Type I error. ",
      "Control rate %.4f, treatment rate %.4f, ",
      "risk ratio %.4f, dispersion %.4f. ",
      "Accrual duration %.1f, trial duration %.1f%s%s%s, ",
      "%s. ",
      "Randomization ratio %.0f:1.\n",
      "%s\n%s"
    ),
    test_type_desc,
    object$k,
    object$n_total[object$k],
    (1 - object$beta) * 100,
    object$alpha * 100,
    inputs$lambda1,
    inputs$lambda2,
    risk_ratio,
    inputs$dispersion,
    sum(inputs$accrual_duration),
    inputs$trial_duration,
    max_followup_str,
    event_gap_str,
    dropout_str,
    exposure_text,
    inputs$ratio,
    upper_spend,
    lower_spend
  )

  class(summary_text) <- "gsNBsummary"
  summary_text
}


#' Print method for gsNBsummary objects
#'
#' @param x An object of class `gsNBsummary`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
#'
#' @examples
#' nb_ss <- sample_size_nbinom(
#'   lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
#'   accrual_rate = 10, accrual_duration = 20, trial_duration = 24
#' )
#' gs_design <- gsNBCalendar(nb_ss, k = 3)
#' s <- summary(gs_design)
#' print(s)
print.gsNBsummary <- function(x, ...) {
  cat(strwrap(x, width = 80), sep = "\n")
  cat("\n")
  invisible(x)
}


#' Convert group sequential design to integer sample sizes
#'
#' Generic function to round sample sizes in a group sequential design to integers.
#' This extends the [gsDesign::toInteger()] function from the gsDesign
#' package to work with `gsNB` objects.
#'
#' @param x An object of class `gsNB` or `gsDesign`.
#' @param ... Additional arguments passed to methods.
#'
#' @return An object of the same class as input with integer sample sizes.
#'
#' @export
#'
#' @examples
#' nb_ss <- sample_size_nbinom(
#'   lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
#'   accrual_rate = 10, accrual_duration = 20, trial_duration = 24
#' )
#' gs_design <- gsNBCalendar(nb_ss, k = 3)
#' gs_integer <- toInteger(gs_design)
toInteger <- function(x, ...) {
  UseMethod("toInteger")
}


#' @describeIn toInteger Method for `gsDesign` objects (calls [gsDesign::toInteger()]).
#' @param ratio Randomization ratio (n2/n1).
#' @param roundUpFinal Logical flag indicating whether to round the final analysis
#'   sample size up to meet or exceed the target size.
#'
#' @export
toInteger.gsDesign <- function(x, ratio = x$ratio, roundUpFinal = TRUE, ...) {
  gsDesign::toInteger(x, ratio = ratio, roundUpFinal = roundUpFinal)
}


#' @describeIn toInteger Method for `gsNB` objects.
#'
#' Rounds sample sizes in a group sequential negative binomial design to integers,
#' respecting the randomization ratio.
#'
#' @param ratio Randomization ratio (n2/n1). If an integer is provided, rounding
#'   is done to a multiple of `ratio + 1`. If `ratio < 1` and `1/ratio` is an
#'   integer (e.g., 1:2 allocation, ratio = 0.5), rounding is done to a multiple
#'   of `1/ratio + 1`. Default uses the ratio from the original design.
#' @param roundUpFinal If `TRUE` (default), the final sample size is rounded
#'   up to ensure the target is met. If `FALSE`, rounding is to the nearest
#'   integer.
#'
#' @details
#' This function rounds sample sizes at each analysis to integers while maintaining
#' the randomization ratio and ensuring monotonically increasing sample sizes across
#' analyses. Only the final analysis sample size is rounded to an integer;
#' interim sample sizes remain as expected (non-integer) values based on
#' the information fraction.
#'
#' When `analysis_times` were provided to [gsNBCalendar()],
#' the statistical information (`n.I`) is recomputed at each analysis
#' time based on the new sample size and expected exposures.
#'
#' @export
toInteger.gsNB <- function(x, ratio = x$nb_design$inputs$ratio, roundUpFinal = TRUE, ...) {
  # Make a copy of the object
  result <- x

  k <- x$k

  # Only round the final sample size to an integer
  # Interim sample sizes remain as expected values (non-integer)
  n_total_new <- x$n_total

  # Round final analysis sample size
  group_size <- 1
  if (is.numeric(ratio) && ratio > 0) {
    if (ratio >= 1 && abs(ratio - round(ratio)) < 1e-8) {
      group_size <- ratio + 1
    } else if (ratio < 1 && abs(1 / ratio - round(1 / ratio)) < 1e-8) {
      group_size <- 1 / ratio + 1
    }
  }

  if (roundUpFinal) {
    n_total_new[k] <- ceiling(x$n_total[k] / group_size) * group_size
  } else {
    n_total_new[k] <- round(x$n_total[k] / group_size) * group_size
  }

  # Recalculate interim sample sizes based on timing and new final sample size
  for (i in seq_len(k - 1)) {
    n_total_new[i] <- x$timing[i] * n_total_new[k]
  }

  # Calculate per-group sample sizes
  result$n_total <- n_total_new
  result$n1 <- n_total_new / (1 + ratio)
  result$n2 <- n_total_new * ratio / (1 + ratio)

  # Update n.I (statistical information) at each analysis

  # If analysis_times (T) are available, compute information at each time
  # accounting for enrollment and exposure
  if (!is.null(x$T)) {
    # Get design parameters from nb_design
    nb <- x$nb_design
    lambda1 <- nb$inputs$lambda1
    lambda2 <- nb$inputs$lambda2
    dispersion <- nb$inputs$dispersion
    max_followup <- nb$inputs$max_followup
    if (is.null(max_followup)) max_followup <- Inf

    # Compute accrual rate based on new final sample size
    # accrual_rate = n_total_final / accrual_duration
    accrual_duration <- nb$accrual_duration
    new_accrual_rate <- n_total_new[k] / accrual_duration

    # Compute information at each analysis time
    info_at_analyses <- numeric(k)
    for (i in seq_len(k)) {
      info_at_analyses[i] <- compute_info_at_time(
        analysis_time = x$T[i],
        accrual_rate = new_accrual_rate,
        accrual_duration = accrual_duration,
        lambda1 = lambda1,
        lambda2 = lambda2,
        dispersion = dispersion,
        ratio = ratio,
        dropout_rate = nb$inputs$dropout_rate,
        event_gap = nb$inputs$event_gap,
        max_followup = max_followup
      )
    }
    result$n.I <- info_at_analyses

    # Update n.fix to be the final analysis information
    result$n.fix <- info_at_analyses[k]

    # Update timing based on new information fractions
    result$timing <- info_at_analyses / info_at_analyses[k]
  } else {
    # If no calendar times, just scale n.I by the sample size increase
    # This assumes information is proportional to sample size (approx true)
    scaling_factor <- n_total_new[k] / x$n_total[k]
    result$n.I <- x$n.I * scaling_factor
    result$n.fix <- result$n.I[k]
  }

  # Update gsDesign object properties
  # We need to update the gsDesign object to reflect the new n.I
  # This is crucial for gsBoundSummary to work correctly

  # Create a new gsDesign object with updated n.I
  # We use the original parameters but new n.I

  # Note: gsDesign() function recalculates bounds based on n.I
  gs_updated <- gsDesign::gsDesign(
    k = k,
    test.type = x$test.type,
    alpha = x$alpha,
    beta = x$beta,
    astar = x$astar,
    delta = x$delta,
    delta1 = x$delta1,
    delta0 = x$delta0,
    n.I = result$n.I,
    maxn.IPlan = result$n.I[k],
    sfu = x$upper$sf,
    sfupar = x$upper$param,
    sfl = x$lower$sf,
    sflpar = x$lower$param,
    tol = x$tol,
    r = x$r,
    usTime = x$usTime,
    lsTime = x$lsTime
  )

  # Copy updated gsDesign slots to result
  result$upper <- gs_updated$upper
  result$lower <- gs_updated$lower
  result$theta <- gs_updated$theta
  result$falseposnb <- gs_updated$falseposnb
  result$en <- gs_updated$en

  return(result)
}
