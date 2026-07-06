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
#'     \item{7}{Two-sided, asymmetric, binding futility and binding harm bounds}
#'     \item{8}{Two-sided, asymmetric, non-binding futility and non-binding harm bounds}
#'   }
#'   Default is 4.
#' @param alpha Type I error (one-sided). Default is 0.025.
#' @param beta Type II error (1 - power). Default is 0.1.
#' @param astar For test.type 5 or 6, allocated Type I error for the lower
#'   bound. For test.type 7 or 8, total probability of crossing the harm bound
#'   under the null. Default is 0 (set to `1 - alpha` internally when needed).
#' @param delta Standardized effect size. Default is 0 (computed from design).
#' @param sfu Spending function for upper bound. Default is `gsDesign::sfHSD`.
#' @param sfupar Parameter for upper spending function. Default is -4.
#' @param sfl Spending function for lower bound. Default is `gsDesign::sfHSD`.
#' @param sflpar Parameter for lower spending function. Default is -2.
#' @param sfharm Spending function for the harm bound (test.type 7 or 8).
#'   Default is `gsDesign::sfHSD`.
#' @param sfharmparam Parameter for the harm spending function. Default is -2.
#' @param testUpper Logical scalar or vector of length `k` specifying which
#'   analyses include an upper (efficacy) bound. `TRUE` (default) means all
#'   analyses. Where `FALSE`, the upper bound is set to `+20` (effectively
#'   `Inf`) and displayed as `NA`. Must be `TRUE` at the final analysis.
#' @param testLower Logical scalar or vector of length `k` specifying which
#'   analyses include a lower (futility) bound. `TRUE` (default) means all
#'   analyses. Where `FALSE`, the lower bound is set to `-20` (effectively
#'   `-Inf`) and displayed as `NA`. Ignored for test.type 1.
#' @param testHarm Logical scalar or vector of length `k` specifying which
#'   analyses include a harm bound (test.type 7 or 8 only). `TRUE` (default)
#'   means all analyses. Where `FALSE`, the harm bound is set to `-20` and
#'   displayed as `NA`.
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
#'     \item{exposure}{A vector with expected average calendar exposure per analysis}
#'     \item{exposure_at_risk1}{A vector with expected at-risk exposure per analysis for group 1}
#'     \item{exposure_at_risk2}{A vector with expected at-risk exposure per analysis for group 2}
#'     \item{variance}{A vector with variance of log rate ratio per analysis}
#'     \item{T}{Calendar time at each analysis (if `analysis_times` provided)}
#'     \item{testUpper}{Logical vector indicating which analyses have an efficacy bound}
#'     \item{testLower}{Logical vector indicating which analyses have a futility bound}
#'     \item{testHarm}{Logical vector indicating which analyses have a harm bound
#'       (test.type 7 or 8 only)}
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
#'
#' # Selective bound testing: futility only at IA1, efficacy deferred to IA2+
#' gs_selective <- gsNBCalendar(nb_ss,
#'   k = 3, test.type = 4,
#'   analysis_times = c(10, 18, 24),
#'   testUpper = c(FALSE, TRUE, TRUE),
#'   testLower = c(TRUE, TRUE, FALSE)
#' )
#' gs_selective
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
  sfharm = gsDesign::sfHSD,
  sfharmparam = -2,
  testUpper = TRUE,
  testLower = TRUE,
  testHarm = TRUE,
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

  # Build the argument list for gsDesign::gsDesign(), conditionally including
  # parameters introduced in newer gsDesign releases so older installations can
  # still use non-harm-bound designs.
  gs_args <- list(
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

  gs_formals <- names(formals(gsDesign::gsDesign))
  if ("sfharm" %in% gs_formals) {
    gs_args$sfharm <- sfharm
    gs_args$sfharmparam <- sfharmparam
  }
  if ("testUpper" %in% gs_formals) {
    gs_args$testUpper <- testUpper
    gs_args$testLower <- testLower
  }
  if ("testHarm" %in% gs_formals) {
    gs_args$testHarm <- testHarm
  }

  gs <- do.call(gsDesign::gsDesign, gs_args)

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
  exposure_at_analysis <- numeric(k)
  exposure_at_risk1_at_analysis <- numeric(k)
  exposure_at_risk2_at_analysis <- numeric(k)
  variance_at_analysis <- numeric(k)
  
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
    exposure_at_analysis[i] <- res_i$exposure[1]
    exposure_at_risk1_at_analysis[i] <- res_i$exposure_at_risk_n1
    exposure_at_risk2_at_analysis[i] <- res_i$exposure_at_risk_n2
    variance_at_analysis[i] <- res_i$variance
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
  result$exposure <- exposure_at_analysis
  result$exposure_at_risk1 <- exposure_at_risk1_at_analysis
  result$exposure_at_risk2 <- exposure_at_risk2_at_analysis
  result$variance <- variance_at_analysis
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
#' Computes the statistical information \eqn{\mathcal{I}} for the log rate
#' ratio \eqn{\theta = \log(\lambda_2/\lambda_1)} at a given calendar analysis
#' time, accounting for staggered enrollment, dropout, maximum follow-up, and
#' event gaps.
#'
#' @param analysis_time Calendar time of the analysis.
#' @param accrual_rate Enrollment rate (subjects per time unit).
#' @param accrual_duration Duration of the enrollment period.
#' @param lambda1 Event rate \eqn{\lambda_1} for group 1 (control).
#' @param lambda2 Event rate \eqn{\lambda_2} for group 2 (treatment).
#' @param dispersion Dispersion parameter \eqn{k} such that
#'   \eqn{\mathrm{Var}(Y) = \mu + k\mu^2}. Can be a vector of length 2.
#' @param ratio Allocation ratio \eqn{r = n_2/n_1}. Default is 1.
#' @param dropout_rate Dropout hazard rate. Default is 0. Can be a vector of
#'   length 2 for group-specific rates (control, treatment).
#' @param event_gap Gap duration after each event. Default is 0.
#' @param max_followup Maximum follow-up time per subject. Default is `Inf`.
#'   Can be a vector of length 2.
#'
#' @details
#' This function delegates to [sample_size_nbinom()] with `power = NULL` and
#' returns \eqn{\mathcal{I} = 1/\mathrm{Var}(\hat\theta)} from the resulting
#' variance. This ensures full consistency with package design calculations,
#' including piecewise accrual, dropout, max follow-up truncation, event-gap
#' correction, and follow-up variability inflation (\eqn{Q_g}).
#'
#' @return The statistical information \eqn{\mathcal{I}} (inverse of variance)
#'   at the analysis time.
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
  ss <- sample_size_nbinom(
    lambda1 = lambda1,
    lambda2 = lambda2,
    rr0 = 1,
    dispersion = dispersion,
    power = NULL,
    alpha = 0.025,
    sided = 1,
    ratio = ratio,
    accrual_rate = accrual_rate,
    accrual_duration = accrual_duration,
    trial_duration = analysis_time,
    dropout_rate = dropout_rate,
    max_followup = max_followup,
    event_gap = event_gap
  )

  1 / ss$variance
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
#' gs_design <- gsNBCalendar(nb_ss, k = 3, analysis_times = c(12, 18, 24))
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
    "7" = "Asymmetric two-sided with binding futility and harm bounds",
    "8" = "Asymmetric two-sided with non-binding futility and harm bounds",
    "Unknown test type"
  )

  # Format exposure text
  exposure_text <- sprintf("average exposure %.2f", nb$exposure[1])
  if (!is.null(inputs$event_gap) && inputs$event_gap > 0 && !is.null(nb$exposure_at_risk_n1)) {
    exposure_text <- sprintf(
      "average exposure (calendar) %.2f, (at-risk n1=%.2f, n2=%.2f)",
      nb$exposure[1], nb$exposure_at_risk_n1, nb$exposure_at_risk_n2
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
    if (is.data.frame(inputs$dropout_rate)) {
      if (any(inputs$dropout_rate$rate > 0)) {
        dropout_str <- ", piecewise dropout"
      }
    } else if (length(inputs$dropout_rate) == 1) {
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
    param_str <- if (!is.null(object$upper$parname)) {
      paste0(" (", paste(object$upper$parname, "=", round(object$upper$param, 2), collapse = ", "), ")")
    } else ""
    paste0("Upper spending: ", object$upper$name, param_str)
  } else {
    "Upper spending: Custom"
  }

  lower_spend <- if (!is.null(object$lower$name)) {
    param_str <- if (!is.null(object$lower$parname)) {
      paste0(" (", paste(object$lower$parname, "=", round(object$lower$param, 2), collapse = ", "), ")")
    } else ""
    paste0("Lower spending: ", object$lower$name, param_str)
  } else {
    "Lower spending: Custom"
  }

  harm_spend <- ""
  if (object$test.type %in% c(7, 8) && !is.null(object$harm)) {
    harm_spend <- if (!is.null(object$harm$name)) {
      param_str <- if (!is.null(object$harm$parname)) {
        paste0(" (", paste(object$harm$parname, "=", round(object$harm$param, 2), collapse = ", "), ")")
      } else ""
      paste0("\nHarm spending: ", object$harm$name, param_str)
    } else {
      "\nHarm spending: Custom"
    }
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
      "%s\n%s%s"
    ),
    test_type_desc,
    object$k,
    object$n_total[object$k],
    (1 - object$beta) * 100,
    object$alpha * 100,
    inputs$lambda1,
    inputs$lambda2,
    risk_ratio,
    inputs$dispersion[1],
    sum(inputs$accrual_duration),
    inputs$trial_duration,
    max_followup_str,
    event_gap_str,
    dropout_str,
    exposure_text,
    inputs$ratio,
    upper_spend,
    lower_spend,
    harm_spend
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
#' gs_design <- gsNBCalendar(nb_ss, k = 3, analysis_times = c(12, 18, 24))
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
#' gs_design <- gsNBCalendar(nb_ss, k = 3, analysis_times = c(12, 18, 24))
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
#' This function rounds the final sample size while maintaining the
#' randomization ratio. When calendar analysis times are available, interim
#' sample sizes remain expected enrollment counts at those calendar times after
#' rescaling the accrual rate to the rounded final sample size.
#'
#' When `analysis_times` were provided to [gsNBCalendar()],
#' expected events, exposure, and statistical information (`n.I`) are recomputed
#' at each analysis time based on the new sample size and expected exposures.
#'
#' @export
toInteger.gsNB <- function(x, ratio = x$nb_design$inputs$ratio, roundUpFinal = TRUE, ...) {
  # Make a copy of the object
  result <- x

  k <- x$k

  group_size <- 1
  if (is.numeric(ratio) && ratio > 0) {
    if (ratio >= 1 && abs(ratio - round(ratio)) < 1e-8) {
      group_size <- ratio + 1
    } else if (ratio < 1 && abs(1 / ratio - round(1 / ratio)) < 1e-8) {
      group_size <- 1 / ratio + 1
    }
  }

  # Round final analysis sample size.
  if (roundUpFinal) {
    n_total_final <- ceiling(x$n_total[k] / group_size) * group_size
  } else {
    n_total_final <- round(x$n_total[k] / group_size) * group_size
  }

  if (!is.null(x$T)) {
    nb <- x$nb_design

    accrual_duration <- nb$accrual_duration
    info_fixed <- 1 / nb$variance
    information_scale <- x$n.I[k] / info_fixed
    rounding_scale <- n_total_final / x$n_total[k]
    new_accrual_rate <- nb$accrual_rate * information_scale * rounding_scale

    n_total_new <- numeric(k)
    info_at_analyses <- numeric(k)
    events_at_analyses <- numeric(k)
    events1_at_analyses <- numeric(k)
    events2_at_analyses <- numeric(k)
    exposure_at_analyses <- numeric(k)
    exposure_at_risk1_at_analyses <- numeric(k)
    exposure_at_risk2_at_analyses <- numeric(k)
    variance_at_analyses <- numeric(k)

    for (i in seq_len(k)) {
      res_i <- sample_size_nbinom(
        lambda1 = nb$inputs$lambda1,
        lambda2 = nb$inputs$lambda2,
        rr0 = nb$inputs$rr0,
        dispersion = nb$inputs$dispersion,
        power = NULL,
        alpha = nb$inputs$alpha,
        sided = nb$inputs$sided,
        ratio = ratio,
        accrual_rate = new_accrual_rate,
        accrual_duration = accrual_duration,
        trial_duration = x$T[i],
        dropout_rate = nb$inputs$dropout_rate,
        max_followup = nb$inputs$max_followup,
        event_gap = nb$inputs$event_gap
      )

      n_total_new[i] <- res_i$n_total
      info_at_analyses[i] <- 1 / res_i$variance
      events_at_analyses[i] <- res_i$total_events
      events1_at_analyses[i] <- res_i$events_n1
      events2_at_analyses[i] <- res_i$events_n2
      exposure_at_analyses[i] <- res_i$exposure[1]
      exposure_at_risk1_at_analyses[i] <- res_i$exposure_at_risk_n1
      exposure_at_risk2_at_analyses[i] <- res_i$exposure_at_risk_n2
      variance_at_analyses[i] <- res_i$variance
    }

    result$n_total <- n_total_new
    result$n1 <- n_total_new / (1 + ratio)
    result$n2 <- n_total_new * ratio / (1 + ratio)
    result$events <- events_at_analyses
    result$events1 <- events1_at_analyses
    result$events2 <- events2_at_analyses
    result$exposure <- exposure_at_analyses
    result$exposure_at_risk1 <- exposure_at_risk1_at_analyses
    result$exposure_at_risk2 <- exposure_at_risk2_at_analyses
    result$variance <- variance_at_analyses
    result$n.I <- info_at_analyses
    result$n.fix <- info_at_analyses[k]
    result$timing <- info_at_analyses / info_at_analyses[k]
  } else {
    # If no calendar times, just scale n.I by the sample size increase
    # This assumes information is proportional to sample size (approx true)
    n_total_new <- x$n_total
    n_total_new[k] <- n_total_final
    for (i in seq_len(k - 1)) {
      n_total_new[i] <- x$timing[i] * n_total_new[k]
    }

    result$n_total <- n_total_new
    result$n1 <- n_total_new / (1 + ratio)
    result$n2 <- n_total_new * ratio / (1 + ratio)

    scaling_factor <- n_total_final / x$n_total[k]
    result$n.I <- x$n.I * scaling_factor
    result$n.fix <- result$n.I[k]
  }

  # Update gsDesign object properties
  # We need to update the gsDesign object to reflect the new n.I
  # This is crucial for gsBoundSummary to work correctly

  # Create a new gsDesign object with updated n.I
  # We use the original parameters but new n.I

  # Note: gsDesign() function recalculates bounds based on n.I
  gs_args <- list(
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
  if (x$test.type %in% c(7, 8) && !is.null(x$harm)) {
    gs_args$sfharm <- x$harm$sf
    gs_args$sfharmparam <- x$harm$param
  }
  if (!is.null(x$testUpper)) gs_args$testUpper <- x$testUpper
  if (!is.null(x$testLower)) gs_args$testLower <- x$testLower
  if (!is.null(x$testHarm)) gs_args$testHarm <- x$testHarm

  gs_updated <- do.call(gsDesign::gsDesign, gs_args)

  # Copy updated gsDesign slots to result
  result$upper <- gs_updated$upper
  result$lower <- gs_updated$lower
  result$theta <- gs_updated$theta
  result$falseposnb <- gs_updated$falseposnb
  result$en <- gs_updated$en
  if (x$test.type %in% c(7, 8)) {
    result$harm <- gs_updated$harm
    result$testHarm <- gs_updated$testHarm
  }
  result$testUpper <- gs_updated$testUpper
  result$testLower <- gs_updated$testLower

  return(result)
}


#' Update group sequential bounds with observed information
#'
#' Given a planned `gsNB` (or `gsDesign`) object and observed statistical
#' information at one or more analyses, recompute the group sequential
#' boundaries and return an updated design object together with a
#' [gsDesign::gsBoundSummary()]-style table.
#'
#' The observed information determines the covariance structure of the test
#' statistics (via the information fraction `timing`), while `spending_time`
#' controls how much of the error-spending budget has been used.
#' When `spending_time` is `NULL` (the default), spending is driven by the
#' observed information fraction.  Supplying an explicit `spending_time` is
#' useful when the monitoring charter specifies calendar-driven spending that
#' differs from the observed information fraction.
#'
#' @param design A `gsNB` or `gsDesign` object produced by
#'   [gsNBCalendar()] (or [gsDesign::gsDesign()]).
#' @param observed_info Numeric vector of observed statistical information at
#'   each analysis conducted so far.
#'   Its length must be between 1 and `design$k`.
#'   If shorter than `design$k`, information at future analyses is projected
#'   from the planned design.
#' @param spending_time Optional numeric vector the same length as
#'   `observed_info` giving the spending time (between 0 and 1) for each
#'   analysis.  When `NULL`, spending time equals the information fraction.
#'   If shorter than `design$k`, future spending times are taken from the
#'   planned design.
#'
#' @return A list with components:
#'   \describe{
#'     \item{design}{The updated `gsDesign` object with recalculated
#'       boundaries.}
#'     \item{bounds}{A data frame from [gsDesign::gsBoundSummary()]
#'       showing Z-boundaries, nominal p-values, approximate treatment effects
#'       at the boundary, and cumulative crossing probabilities at each
#'       analysis.}
#'     \item{info}{A data frame with one row per analysis containing the
#'       information fraction (`IF`), spending time (`spending_time`),
#'       upper and lower Z-boundaries, and cumulative upper and lower
#'       spending.}
#'   }
#'
#' @export
#'
#' @examples
#' library(gsDesign)
#' nb_ss <- sample_size_nbinom(
#'   lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
#'   accrual_rate = 10, accrual_duration = 20, trial_duration = 24
#' )
#' gs <- gsNBCalendar(nb_ss, k = 3, analysis_times = c(12, 18, 24))
#'
#' # After observing information at the first interim
#' upd <- update_gsNB(gs, observed_info = gs$n.I[1] * 0.95)
#' upd$bounds
#' upd$info
update_gsNB <- function(design, observed_info, spending_time = NULL) {
  # --- input validation -------------------------------------------------------
  if (!inherits(design, "gsDesign")) {
    stop("`design` must be a gsDesign or gsNB object.", call. = FALSE)
  }
  n_obs <- length(observed_info)
  if (n_obs < 1 || n_obs > design$k) {
    stop(
      sprintf(
        "`observed_info` must have between 1 and %d elements (design$k).",
        design$k
      ),
      call. = FALSE
    )
  }
  if (!is.null(spending_time)) {
    if (length(spending_time) != n_obs) {
      stop(
        "`spending_time` must be the same length as `observed_info`.",
        call. = FALSE
      )
    }
    if (any(spending_time <= 0 | spending_time > 1)) {
      stop("`spending_time` values must be in (0, 1].", call. = FALSE)
    }
  }

  # --- build timing vectors ---------------------------------------------------
  max_info <- design$n.fix

  # Information fractions: observed for conducted analyses, planned for future

  timing <- design$timing
  timing[seq_len(n_obs)] <- observed_info / max_info

  # If the last observation is the final analysis, set to 1
  if (n_obs == design$k) timing[n_obs] <- max(timing[n_obs], 1)

  # Enforce monotonicity
  for (i in seq_along(timing)[-1]) {
    if (timing[i] <= timing[i - 1]) {
      timing[i] <- timing[i - 1] + 1e-4
    }
  }

  # Spending time
  us_time <- design$usTime
  ls_time <- design$lsTime
  if (!is.null(spending_time)) {
    if (is.null(us_time)) us_time <- design$timing
    if (is.null(ls_time)) ls_time <- design$timing
    us_time[seq_len(n_obs)] <- spending_time
    ls_time[seq_len(n_obs)] <- spending_time
  }

  # --- rebuild design ---------------------------------------------------------
  gs_args <- list(
    k         = design$k,
    test.type = design$test.type,
    alpha     = design$alpha,
    beta      = design$beta,
    astar     = design$astar,
    delta     = design$delta,
    sfu       = design$upper$sf,
    sfupar    = design$upper$param,
    sfl       = design$lower$sf,
    sflpar    = design$lower$param,
    tol       = design$tol,
    r         = design$r,
    n.fix     = max_info,
    timing    = timing,
    usTime    = us_time,
    lsTime    = ls_time
  )
  if (design$test.type %in% c(7, 8) && !is.null(design$harm)) {
    gs_args$sfharm     <- design$harm$sf
    gs_args$sfharmparam <- design$harm$param
  }
  gs_formals <- names(formals(gsDesign::gsDesign))
  if ("testUpper" %in% gs_formals && !is.null(design$testUpper)) {
    gs_args$testUpper <- design$testUpper
  }
  if ("testLower" %in% gs_formals && !is.null(design$testLower)) {
    gs_args$testLower <- design$testLower
  }
  if ("testHarm" %in% gs_formals && !is.null(design$testHarm)) {
    gs_args$testHarm <- design$testHarm
  }

  updated <- do.call(gsDesign::gsDesign, gs_args)

  # --- summary table ----------------------------------------------------------
  bounds <- gsDesign::gsBoundSummary(updated)

  # --- compact info table -----------------------------------------------------
  eff_spend <- cumsum(updated$upper$spend)
  fut_spend <- cumsum(updated$lower$spend)
  st <- if (!is.null(us_time)) us_time else timing

  info_df <- data.frame(
    Analysis      = seq_len(design$k),
    IF            = round(timing, 4),
    spending_time = round(st, 4),
    upper_bound   = round(updated$upper$bound, 4),
    lower_bound   = round(updated$lower$bound, 4),
    cum_upper_spend = round(eff_spend, 6),
    cum_lower_spend = round(fut_spend, 6)
  )

  list(design = updated, bounds = bounds, info = info_df)
}
