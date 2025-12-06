#' Group Sequential Design for Negative Binomial Outcomes
#'
#' Creates a group sequential design for negative binomial outcomes based on
#' sample size calculations from \code{\link{sample_size_nbinom}}.
#'
#' @param x An object of class \code{sample_size_nbinom_result} from
#'   \code{\link{sample_size_nbinom}}.
#' @param k Number of analyses (interim + final). Default is 3.
#' @param test.type Test type as in \code{\link[gsDesign]{gsDesign}}:
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
#' @param timing Timing of interim analyses. May be a vector of length k-1
#'   with values between 0 and 1 representing information fractions.
#'   Default is 1 (equally spaced).
#' @param sfu Spending function for upper bound. Default is \code{gsDesign::sfHSD}.
#' @param sfupar Parameter for upper spending function. Default is -4.
#' @param sfl Spending function for lower bound. Default is \code{gsDesign::sfHSD}.
#' @param sflpar Parameter for lower spending function. Default is -2.
#' @param tol Tolerance for convergence. Default is 1e-06.
#' @param r Integer controlling grid size for numerical integration.
#'   Default is 18.
#' @param usTime Spending time for upper bound (optional).
#' @param lsTime Spending time for lower bound (optional).
#' @param analysis_times Optional vector of calendar times for each analysis.
#'   If provided, must have length k. These times are stored in the \code{T}
#'   element and displayed by \code{\link[gsDesign]{gsBoundSummary}}.
#'
#' @return An object of class \code{gsNB} which inherits from \code{gsDesign}
#'   and \code{sample_size_nbinom_result}. Contains all elements from
#'   \code{gsDesign::gsDesign()} plus:
#'   \describe{
#'     \item{nb_design}{The original \code{sample_size_nbinom_result} object}
#'     \item{n1}{Sample size per analysis for group 1}
#'     \item{n2}{Sample size per analysis for group 2}
#'     \item{T}{Calendar time at each analysis (if \code{analysis_times} provided)}
#'   }
#'
#' @references
#' Jennison, C. and Turnbull, B.W. (2000), \emph{Group Sequential Methods with
#' Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
#'
#' @export
#'
#' @examples
#' # First create a sample size calculation
#' nb_ss <- sample_size_nbinom(
#'   lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
#'   accrual_rate = 10, accrual_duration = 20, trial_duration = 24
#' )
#'
#' # Then create a group sequential design with analysis times
#' gs_design <- gsNBCalendar(nb_ss, k = 3, test.type = 4,
#'                          analysis_times = c(10, 18, 24))
#'
#' @importFrom gsDesign gsDesign sfHSD
gsNBCalendar <- function(x,
                         k = 3,
                         test.type = 4,
                         alpha = 0.025,
                         beta = 0.1,
                         astar = 0,
                         delta = 0,
                         timing = 1,
                         sfu = gsDesign::sfHSD,
                         sfupar = -4,
                         sfl = gsDesign::sfHSD,
                         sflpar = -2,
                         tol = 1e-06,
                         r = 18,
                         usTime = NULL,
                         lsTime = NULL,
                         analysis_times = NULL) {
  # Validate input

  if (!inherits(x, "sample_size_nbinom_result")) {
    stop("x must be an object of class 'sample_size_nbinom_result'")
  }

  # Calculate delta1 based on the log rate ratio
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
    timing = timing,
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
  n_total_fixed <- x$n_total
  ratio <- x$inputs$ratio

  # Calculate cumulative sample sizes at each analysis
  # n.I / n.fix gives the information fraction, multiply by fixed sample size
  n_cumulative <- (gs$n.I / info_fixed) * n_total_fixed

  # Per-group sample sizes (cumulative)
  n1_cumulative <- n_cumulative / (1 + ratio)
  n2_cumulative <- n_cumulative * ratio / (1 + ratio)

  # Build result object

  # Start with the gsDesign object
  result <- gs

  # Add negative binomial specific components
  result$nb_design <- x
  result$n1 <- n1_cumulative
  result$n2 <- n2_cumulative
  result$n_total <- n_cumulative
  result$ratio <- ratio

 # Add calendar times if provided (for gsBoundSummary display)
  if (!is.null(analysis_times)) {
    if (length(analysis_times) != k) {
      stop("analysis_times must have length k (", k, ")")
    }
    result$T <- analysis_times
  }

  # Set the class to inherit from both gsDesign and sample_size_nbinom_result
  class(result) <- c("gsNB", "gsDesign", "sample_size_nbinom_result")

  return(result)
}


#' Compute Statistical Information at Analysis Time
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
#'
#' @return The statistical information (inverse of variance) at the analysis time.
#'
#' @details
#' For subjects enrolled at time \code{t}, their exposure at analysis time \code{T}
#' is \code{T - t}. The function integrates over all enrolled subjects to compute
#' the total expected events and variance of the log rate ratio.
#'
#' The variance formula for the log rate ratio is:
#' \deqn{Var = \sum_i \frac{1/\mu_{1i} + k}{n_{1i}} + \frac{1/\mu_{2i} + k}{n_{2i}}}
#'
#' where \eqn{\mu_{ji} = \lambda_j \times exposure_i} and \eqn{k} is the dispersion.
#'
#' @export
#'
#' @examples
#' # Compute information at month 12 for a trial with:
#' # - 10 subjects/month enrollment rate
#' # - 20 month enrollment period
#' # - Control event rate 0.5, treatment rate 0.3
#' # - Dispersion 0.1
#' compute_info_at_time(
#'   analysis_time = 12,
#'   accrual_rate = 10,
#'   accrual_duration = 20,
#'   lambda1 = 0.5,
#'   lambda2 = 0.3,
#'   dispersion = 0.1
#' )
compute_info_at_time <- function(analysis_time, accrual_rate, accrual_duration,
                                  lambda1, lambda2, dispersion, ratio = 1) {
  # Number of subjects enrolled by analysis_time
  enrollment_time <- min(analysis_time, accrual_duration)
  n_total <- accrual_rate * enrollment_time
  n1 <- n_total / (1 + ratio)
  n2 <- n_total * ratio / (1 + ratio)


  # For uniform enrollment, the average exposure is:

  # If analysis_time <= accrual_duration: avg_exposure = analysis_time / 2
  # If analysis_time > accrual_duration: 
  #   Subjects enrolled at time t have exposure = analysis_time - t
  #   For t in [0, accrual_duration], avg_exposure = analysis_time - accrual_duration/2

  if (analysis_time <= accrual_duration) {
    # Still enrolling - average exposure is half the analysis time
    avg_exposure <- analysis_time / 2
  } else {
    # Enrollment complete - average exposure is analysis_time minus midpoint of enrollment
    avg_exposure <- analysis_time - accrual_duration / 2
  }

  # Expected events per subject
  mu1 <- lambda1 * avg_exposure
  mu2 <- lambda2 * avg_exposure

  # Variance of log rate ratio
  # Var(log(lambda2/lambda1)) = (1/mu1 + k)/n1 + (1/mu2 + k)/n2
  k <- dispersion
  variance <- (1 / mu1 + k) / n1 + (1 / mu2 + k) / n2


  # Information is inverse of variance
  info <- 1 / variance

  return(info)
}


#' Summary for gsNB Objects
#'
#' Provides a textual summary of a group sequential design for negative binomial
#' outcomes, similar to the summary provided by \code{\link[gsDesign]{gsDesign}}.
#' For tabular output, use \code{\link[gsDesign]{gsBoundSummary}} directly on
#' the gsNB object.
#'
#' @param object An object of class \code{gsNB}.
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
#' # For tabular bounds summary, use gsBoundSummary directly:
#' # gsDesign::gsBoundSummary(gs_design)
summary.gsNB <- function(object, ...) {
  # Extract negative binomial design information
  nb <- object$nb_design
  inputs <- nb$inputs

  # Calculate risk ratio (lambda2/lambda1)
  risk_ratio <- inputs$lambda2 / inputs$lambda1

  # Determine test type description
  test_type_desc <- switch(
    as.character(object$test.type),
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
    exposure_text <- sprintf("average exposure (calendar) %.2f, (at-risk n1=%.2f, n2=%.2f)", 
                             nb$exposure, nb$exposure_at_risk_n1, nb$exposure_at_risk_n2)
  }

  # Build the summary text
  max_followup_str <- if (!is.null(inputs$max_followup)) sprintf(", max follow-up %.1f", inputs$max_followup) else ""
  
  summary_text <- sprintf(
    paste0(
      "%s group sequential design for negative binomial outcomes, ",
      "%d analyses, ",
      "total sample size %.1f, ",
      "%.0f percent power, ",
      "%.1f percent (1-sided) Type I error. ",
      "Control rate %.4f, treatment rate %.4f, ",
      "risk ratio %.4f, dispersion %.4f. ",
      "Accrual duration %.1f, trial duration %.1f%s, ",
      "%s. ",
      "Randomization ratio %.0f:1."
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
    exposure_text,
    inputs$ratio
  )

  class(summary_text) <- "gsNBsummary"
  summary_text
}


#' Print Method for gsNBsummary Objects
#'
#' @param x An object of class \code{gsNBsummary}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.gsNBsummary <- function(x, ...) {
  cat(strwrap(x, width = 80), sep = "\n")
  cat("\n")
  invisible(x)
}


#' Convert Group Sequential Design to Integer Sample Sizes
#'
#' Generic function to round sample sizes in a group sequential design to integers.
#' This extends the \code{\link[gsDesign]{toInteger}} function from the gsDesign
#' package to work with \code{gsNB} objects.
#'
#' @param x An object of class \code{gsNB} or \code{gsDesign}.
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


#' @describeIn toInteger Method for gsDesign objects (calls gsDesign::toInteger)
#' @param ratio Randomization ratio. See \code{\link[gsDesign]{toInteger}}.
#' @param roundUpFinal Logical. See \code{\link[gsDesign]{toInteger}}.
#' @export
toInteger.gsDesign <- function(x, ratio = x$ratio, roundUpFinal = TRUE, ...) {
 gsDesign::toInteger(x, ratio = ratio, roundUpFinal = roundUpFinal)
}


#' @describeIn toInteger Method for gsNB objects
#'
#' Rounds sample sizes in a group sequential negative binomial design to integers,
#' respecting the randomization ratio.
#'
#' @param ratio Randomization ratio (n2/n1). If an integer is provided, rounding
#'   is done to a multiple of \code{ratio + 1}. Default uses the ratio from the
#'   original design.
#' @param roundUpFinal If \code{TRUE} (default), the final sample size is rounded
#'   up to ensure the target is met. If \code{FALSE}, rounding is to the nearest
#'   integer.
#'
#' @details
#' This function rounds sample sizes at each analysis to integers while maintaining
#' the randomization ratio and ensuring monotonically increasing sample sizes across
#' analyses. Only the final analysis sample size is rounded to an integer;
#' interim sample sizes remain as expected (non-integer) values based on
#' the information fraction.
#'
#' When \code{analysis_times} were provided to \code{\link{gsNBCalendar}},
#' the statistical information (\code{n.I}) is recomputed at each analysis
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
  if (roundUpFinal) {
    if (is.numeric(ratio) && ratio == floor(ratio) && ratio >= 0) {
      # Round up to nearest multiple of (ratio + 1)
      group_size <- ratio + 1
      n_total_new[k] <- ceiling(x$n_total[k] / group_size) * group_size
    } else {
      n_total_new[k] <- ceiling(x$n_total[k])
    }
  } else {
    if (is.numeric(ratio) && ratio == floor(ratio) && ratio >= 0) {
      group_size <- ratio + 1
      n_total_new[k] <- round(x$n_total[k] / group_size) * group_size
    } else {
      n_total_new[k] <- round(x$n_total[k])
    }
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
        ratio = ratio
      )
    }
    result$n.I <- info_at_analyses

    # Update n.fix to be the final analysis information
    result$n.fix <- info_at_analyses[k]

    # Update timing based on new information fractions
    result$timing <- info_at_analyses / info_at_analyses[k]
  } else {
    # No analysis times - use simple scaling
    result$n.I <- result$timing * result$n.fix
  }

  result
}