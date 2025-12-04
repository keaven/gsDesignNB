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
#' @param n.fix Sample size inflation factor. Default is 1.
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
#'
#' @return An object of class \code{gsNB} which inherits from \code{gsDesign}
#'   and \code{sample_size_nbinom_result}. Contains all elements from
#'   \code{gsDesign::gsDesign()} plus:
#'   \describe{
#'     \item{nb_design}{The original \code{sample_size_nbinom_result} object}
#'     \item{n1}{Sample size per analysis for group 1}
#'     \item{n2}{Sample size per analysis for group 2}
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
#' # Then create a group sequential design
#' gs_design <- gsNBCalendar(nb_ss, k = 3, test.type = 4)
#'
#' @importFrom gsDesign gsDesign sfHSD
gsNBCalendar <- function(x,
                         k = 3,
                         test.type = 4,
                         alpha = 0.025,
                         beta = 0.1,
                         astar = 0,
                         delta = 0,
                         n.fix = 1,
                         timing = 1,
                         sfu = gsDesign::sfHSD,
                         sfupar = -4,
                         sfl = gsDesign::sfHSD,
                         sflpar = -2,
                         tol = 1e-06,
                         r = 18,
                         usTime = NULL,
                         lsTime = NULL) {
  # Validate input

  if (!inherits(x, "sample_size_nbinom_result")) {
    stop("x must be an object of class 'sample_size_nbinom_result'")
  }

  # Call gsDesign with the provided parameters

  gs <- gsDesign::gsDesign(
    k = k,
    test.type = test.type,
    alpha = alpha,
    beta = beta,
    astar = astar,
    delta = delta,
    n.fix = n.fix,
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
  # gs$n.I contains the cumulative sample size at each analysis
  # We need to scale this by the fixed sample size from nb design
  n_total_fixed <- x$n_total
  ratio <- x$inputs$ratio

  # gs$n.I is the information (sample size) at each analysis relative to n.fix


  # Calculate cumulative sample sizes at each analysis
  n_cumulative <- gs$n.I * n_total_fixed

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

  # Set the class to inherit from both gsDesign and sample_size_nbinom_result
  class(result) <- c("gsNB", "gsDesign", "sample_size_nbinom_result")

  return(result)
}

