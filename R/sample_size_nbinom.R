#' Sample Size Calculation for Negative Binomial Distribution
#'
#' Computes the sample size for comparing two treatment groups assuming a negative
#' binomial distribution for the outcome.
#'
#' @param lambda1 Rate in group 1 (control).
#' @param lambda2 Rate in group 2 (treatment).
#' @param dispersion Dispersion parameter \code{k} such that \eqn{Var(Y) = \mu + k \mu^2}.
#'   Note that this is equivalent to \code{1/size} in R's \code{rnbinom} parameterization.
#' @param power Power of the test (1 - beta). Default is 0.9.
#' @param alpha Significance level. Default is 0.025.
#' @param sided One-sided or two-sided test. 1 for one-sided, 2 for two-sided. Default is 1.
#' @param ratio Allocation ratio n2/n1. Default is 1.
#' @param accrual_rate Vector of accrual rates (patients per unit time).
#' @param accrual_duration Vector of durations for each accrual rate. Must be same length
#'   as \code{accrual_rate}.
#' @param trial_duration Total planned duration of the trial.
#' @param dropout_rate Dropout rate (hazard rate). Default is 0.
#' @param max_followup Maximum follow-up time for any patient. Default is NULL (infinite).
#' @param event_gap Gap duration after each event during which no new events are counted.
#'   Default is NULL (no gap). If provided, the effective event rate is reduced.
#' @param method Method for sample size calculation. "zhu" for Zhu and Lakkis (2014),
#'   "friede" for Friede and Schmidli (2010) / Mütze et al. (2018), or
#'   "AZG" for the Anderson-Zhang-Gemini method which adjusts exposure for event gaps
#'   rather than rates.
#'
#' @return An object of class \code{sample_size_nbinom_result}, which is a list containing:
#' \describe{
#'   \item{inputs}{Named list of the original function arguments.}
#'   \item{n1}{Sample size for group 1}
#'   \item{n2}{Sample size for group 2}
#'   \item{n_total}{Total sample size}
#'   \item{exposure}{Average exposure time used in calculation (calendar time)}
#'   \item{exposure1}{Average at-risk exposure time for group 1 (if AZG method)}
#'   \item{exposure2}{Average at-risk exposure time for group 2 (if AZG method)}
#' }
#'
#' @references
#' Zhu, H., & Lakkis, H. (2014). Sample size calculation for comparing two negative
#' binomial rates in clinical trials. \emph{Statistics in Biopharmaceutical Research},
#' 6(1), 107-115.
#'
#' Friede, T., & Schmidli, H. (2010). Sample size estimation for clinical trials
#' with negative binomial rates. \emph{Methods of Information in Medicine},
#' 49(6), 623-631.
#'
#' Mütze, T., Glimm, E., Schmidli, H., & Friede, T. (2018). Group sequential designs
#' for negative binomial outcomes. \emph{Statistical Methods in Medical Research},
#' 27(10), 2978-2993.
#'
#' @seealso
#' \code{vignette("sample_size_nbinom", package = "gsDesignNB")} for a detailed explanation of the methodology and the difference between "zhu" and "AZG" methods.
#'
#' @export
#'
#' @examples
#' # Calculate sample size for lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1
#' # with fixed recruitment of 10/month for 20 months, 24 month trial duration
#' x <- sample_size_nbinom(
#'   lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
#'   accrual_rate = 10, accrual_duration = 20, trial_duration = 24
#' )
#' class(x)
#' summary(x)
#'
#' # With piecewise accrual
#' # 5 patients/month for 3 months, then 10 patients/month for 3 months
#' # Trial ends at month 12.
#' x2 <- sample_size_nbinom(
#'   lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
#'   accrual_rate = c(5, 10), accrual_duration = c(3, 3),
#'   trial_duration = 12
#' )
#' summary(x2)
#'
#' @importFrom stats pnorm qnorm
sample_size_nbinom <- function(lambda1, lambda2, dispersion, power = NULL,
                               alpha = 0.025, sided = 1, ratio = 1,
                               accrual_rate, accrual_duration,
                               trial_duration, dropout_rate = 0,
                               max_followup = NULL, event_gap = NULL, method = "zhu") {
  if (lambda1 <= 0 || lambda2 <= 0) {
    stop("Rates lambda1 and lambda2 must be positive.")
  }
  if (dispersion < 0) {
    stop("Dispersion parameter must be non-negative.")
  }
  if (!is.null(power) && (power <= 0 || power >= 1)) {
    stop("Power must be between 0 and 1.")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("Alpha must be between 0 and 1.")
  }
  if (dropout_rate < 0) {
    stop("Dropout rate must be non-negative.")
  }
  if (!is.null(max_followup) && max_followup <= 0) {
    stop("max_followup must be positive.")
  }

  power_input <- power
  inputs <- list(
    lambda1 = lambda1,
    lambda2 = lambda2,
    dispersion = dispersion,
    power = power_input,
    alpha = alpha,
    sided = sided,
    ratio = ratio,
    accrual_rate = accrual_rate,
    accrual_duration = accrual_duration,
    trial_duration = trial_duration,
    dropout_rate = dropout_rate,
    max_followup = max_followup,
    event_gap = event_gap,
    method = method
  )

  # Determine mode: Calculate N or Calculate Power
  # We don't have 'mode' logic based on missing accrual_rate anymore since it's required.
  # If power is NULL, solve for N. If power is provided, solve for Power (implied).
  # BUT: If accrual_rate is fixed (provided), we effectively have a fixed N.
  # In typical 'solve_n' mode, we want to SCALE the accrual rate to achieve the target power.
  # In 'solve_power' mode, we keep accrual rate fixed and compute power.

  mode <- "solve_n"
  if (!is.null(power)) {
    # We are solving for N (scaling accrual) to meet this power
    mode <- "solve_n"
  } else {
    # Power is null, so we calculate power for the fixed accrual
    mode <- "solve_power"
    # Default power? No, we are calculating it.
  }

  # Helper function for average exposure over [u_min, u_max] given dropout_rate
  avg_exp_func <- function(u_min, u_max, dropout_rate) {
    if (u_min >= u_max) {
      return(0)
    } # Should not happen in valid segment
    if (dropout_rate == 0) {
      return((u_min + u_max) / 2)
    } else {
      term1 <- 1 / dropout_rate
      term2 <- (exp(-dropout_rate * u_min) - exp(-dropout_rate * u_max)) / (dropout_rate^2 * (u_max - u_min))
      return(term1 - term2)
    }
  }

  # Calculate average exposure
  current_time <- 0
  total_n_accrual <- 0
  total_exposure_mass <- 0

  # Handle variable accrual to calculate average exposure
  # Accrual arguments are now mandatory, so no checks for is.null needed (R will error if missing)
  
  if (length(accrual_rate) != length(accrual_duration)) {
    stop("accrual_rate and accrual_duration must have the same length.")
  }

  total_accrual_time <- sum(accrual_duration)
  if (total_accrual_time > trial_duration) {
    stop("Total accrual duration cannot exceed trial duration.")
  }

  for (i in seq_along(accrual_rate)) {
    r <- accrual_rate[i]
    d <- accrual_duration[i]

    n_seg <- r * d
    if (n_seg > 0) {
      # Potential follow-up range (administrative censoring only)
      u_max <- trial_duration - current_time
      u_min <- trial_duration - (current_time + d)

      avg_followup <- 0

      if (is.null(max_followup) || is.infinite(max_followup) || u_max <= max_followup) {
        # Case 1: No truncation by max_followup (or negligible)
        avg_followup <- avg_exp_func(u_min, u_max, dropout_rate)
      } else if (u_min >= max_followup) {
        # Case 2: All truncated by max_followup
        # Effectively fixed exposure of max_followup
        if (dropout_rate == 0) {
          avg_followup <- max_followup
        } else {
          avg_followup <- (1 - exp(-dropout_rate * max_followup)) / dropout_rate
        }
      } else {
        # Case 3: Split
        tau_star <- trial_duration - max_followup
        d1 <- tau_star - current_time
        d2 <- (current_time + d) - tau_star

        if (dropout_rate == 0) {
          avg_1 <- max_followup
        } else {
          avg_1 <- (1 - exp(-dropout_rate * max_followup)) / dropout_rate
        }

        avg_2 <- avg_exp_func(u_min, max_followup, dropout_rate)
        avg_followup <- (d1 * avg_1 + d2 * avg_2) / d
      }

      total_n_accrual <- total_n_accrual + n_seg
      total_exposure_mass <- total_exposure_mass + n_seg * avg_followup
    }
    current_time <- current_time + d
  }

  if (total_n_accrual == 0) {
    stop("Accrual results in 0 patients.")
  }

  exposure_calendar <- total_exposure_mass / total_n_accrual

  # Setup effective rates and exposures based on method and event_gap
  exposure1 <- exposure_calendar
  exposure2 <- exposure_calendar
  lambda1_eff <- lambda1
  lambda2_eff <- lambda2

  if (method == "AZG") {
    # Anderson-Zhang-Gemini method: adjust exposure, keep rates raw
    if (!is.null(event_gap) && !is.na(event_gap) && event_gap > 0) {
      exposure1 <- exposure_calendar / (1 + lambda1 * event_gap)
      exposure2 <- exposure_calendar / (1 + lambda2 * event_gap)
    }
    # rates remain unadjusted
  } else {
    # zhu or friede: adjust rates, keep exposure calendar
    if (!is.null(event_gap) && !is.na(event_gap) && event_gap > 0) {
      lambda1_eff <- lambda1 / (1 + lambda1 * event_gap)
      lambda2_eff <- lambda2 / (1 + lambda2 * event_gap)
    }
  }

  mu1 <- lambda1_eff * exposure1
  mu2 <- lambda2_eff * exposure2
  k <- dispersion

  z_alpha <- qnorm(1 - alpha / sided)

  n1_c <- 0
  n2_c <- 0
  n_total_c <- 0
  computed_accrual_rate <- NULL

  if (mode == "solve_n") {
    z_beta <- qnorm(power)

    if (method == "zhu" || method == "AZG") {
      num <- (z_alpha + z_beta)^2 * ((1 / mu1 + k) + (1 / ratio) * (1 / mu2 + k))
      den <- (log(lambda1 / lambda2))^2
      n1 <- num / den
      n2 <- n1 * ratio
    } else if (method == "friede") {
      # Variance term for Friede
      p1 <- 1 / (1 + ratio)
      p2 <- ratio / (1 + ratio)
      V <- (1 / mu1 + k) / p1 + (1 / mu2 + k) / p2
      num <- (z_alpha + z_beta)^2 * V
      den <- (log(lambda1 / lambda2))^2
      n_total <- num / den
      n1 <- n_total * p1
      n2 <- n_total * p2
    } else {
      stop("Unknown method. Choose 'zhu', 'friede', or 'AZG'.")
    }

    n1_c <- ceiling(n1)
    n2_c <- ceiling(n2)
    n_total_c <- n1_c + n2_c

    # Scaling accrual logic
    if (!is.null(accrual_rate)) {
      scaling_factor <- n_total_c / total_n_accrual
      computed_accrual_rate <- accrual_rate * scaling_factor
    }
  } else {
    # solve_power
    computed_accrual_rate <- accrual_rate
    n_total_c <- total_n_accrual
    n1_c <- n_total_c / (1 + ratio)
    n2_c <- n_total_c * ratio / (1 + ratio)

    if (method == "zhu" || method == "AZG") {
      # z_beta = sqrt( n1 * (log(mu1/mu2))^2 / V ) - z_alpha
      V <- (1 / mu1 + k) + (1 / ratio) * (1 / mu2 + k)
      z_beta <- sqrt(n1_c * (log(lambda1 / lambda2))^2 / V) - z_alpha
    } else if (method == "friede") {
      # z_beta = sqrt( n_total * (log(lambda1/lambda2))^2 / V_bar ) - z_alpha
      p1 <- 1 / (1 + ratio)
      p2 <- ratio / (1 + ratio)
      V <- (1 / mu1 + k) / p1 + (1 / mu2 + k) / p2
      z_beta <- sqrt(n_total_c * (log(lambda1 / lambda2))^2 / V) - z_alpha
    } else {
      stop("Unknown method. Choose 'zhu', 'friede', or 'AZG'.")
    }

    power <- pnorm(z_beta)
  }

  # Calculate variance of log rate ratio log(lambda1/lambda2) with calculated sample size
  # Note: For solve_power, we use the non-rounded n1_c/n2_c derived from accrual if they are floats?
  # Accrual rates imply sample size. If sample size is not integer, what is the variance?
  # Usually variance depends on N.
  # I'll use the n1_c and n2_c as calculated (which might be floats in solve_power mode).

  variance <- (1 / mu1 + k) / n1_c + (1 / mu2 + k) / n2_c

  # Calculate expected events
  # Expected events = n * mu = n * lambda * exposure
  # Note: mu1 and mu2 are already lambda * exposure
  events_n1 <- n1_c * mu1
  events_n2 <- n2_c * mu2
  total_events <- events_n1 + events_n2

  result <- c(
    list(inputs = inputs),
    list(
      n1 = n1_c,
      n2 = n2_c,
      n_total = n_total_c,
      alpha = alpha,
      sided = sided,
      power = power,
      exposure = exposure_calendar,
      exposure1 = exposure1,
      exposure2 = exposure2,
      events_n1 = events_n1,
      events_n2 = events_n2,
      total_events = total_events,
      variance = variance,
      accrual_rate = computed_accrual_rate,
      accrual_duration = accrual_duration
    )
  )
  class(result) <- c("sample_size_nbinom_result", "list")
  result
}


#' Print Method for sample_size_nbinom_result Objects
#'
#' Prints a concise summary of the sample size calculation results.
#'
#' @param x An object of class \code{sample_size_nbinom_result}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.sample_size_nbinom_result <- function(x, ...) {
  cat("Sample Size for Negative Binomial Outcome\n")
  cat("==========================================\n\n")
  
  cat(sprintf("Method:          %s\n", x$inputs$method))
  cat(sprintf("Sample size:     n1 = %d, n2 = %d, total = %d\n", 
              x$n1, x$n2, x$n_total))
  cat(sprintf("Expected events: %.1f (n1: %.1f, n2: %.1f)\n",
              x$total_events, x$events_n1, x$events_n2))
  cat(sprintf("Power: %.0f%%, Alpha: %.3f (%d-sided)\n",
              x$power * 100, x$alpha, x$sided))
  cat(sprintf("Rates: control = %.4f, treatment = %.4f (RR = %.4f)\n",
              x$inputs$lambda1, x$inputs$lambda2, 
              x$inputs$lambda2 / x$inputs$lambda1))
  
  if (x$inputs$method == "AZG" && !is.null(x$inputs$event_gap) && x$inputs$event_gap > 0) {
    cat(sprintf("Dispersion: %.4f, Avg exposure (calendar): %.2f\n",
                x$inputs$dispersion, x$exposure))
    cat(sprintf("Avg exposure (at-risk): n1 = %.2f, n2 = %.2f\n",
                x$exposure1, x$exposure2))
  } else {
    cat(sprintf("Dispersion: %.4f, Avg exposure: %.2f\n",
                x$inputs$dispersion, x$exposure))
  }
  
  cat(sprintf("Accrual: %.1f, Trial duration: %.1f\n",
              sum(x$inputs$accrual_duration), x$inputs$trial_duration))
  
  invisible(x)
}


#' Summary for sample_size_nbinom_result Objects
#'
#' Provides a textual summary of the sample size calculation for negative binomial
#' outcomes, similar to the summary for gsNB objects.
#'
#' @param object An object of class \code{sample_size_nbinom_result}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A character string summarizing the design (invisibly). The summary
#'   is also printed to the console.
#'
#' @export
#'
#' @examples
#' x <- sample_size_nbinom(
#'   lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
#'   accrual_rate = 10, accrual_duration = 20, trial_duration = 24
#' )
#' class(x)
#' summary(x)
summary.sample_size_nbinom_result <- function(object, ...) {
  inputs <- object$inputs
  risk_ratio <- inputs$lambda2 / inputs$lambda1
  
  # Build the summary text
  summary_text <- sprintf(
    paste0(
      "Fixed sample size design for negative binomial outcome (%s method), ",
      "total sample size %d (n1=%d, n2=%d), ",
      "%.0f percent power, ",
      "%.1f percent (%d-sided) Type I error. ",
      "Control rate %.4f, treatment rate %.4f, ",
      "risk ratio %.4f, dispersion %.4f. ",
      "Accrual duration %.1f, trial duration %.1f, ",
      "average exposure %.2f. ",
      "Expected events %.1f. ",
      "Randomization ratio %.0f:1."
    ),
    inputs$method,
    object$n_total,
    object$n1,
    object$n2,
    object$power * 100,
    object$alpha * 100,
    inputs$sided,
    inputs$lambda1,
    inputs$lambda2,
    risk_ratio,
    inputs$dispersion,
    sum(inputs$accrual_duration),
    inputs$trial_duration,
    object$exposure,
    object$total_events,
    inputs$ratio
  )
  
  class(summary_text) <- "sample_size_nbinom_summary"
  summary_text
}


#' Print Method for sample_size_nbinom_summary Objects
#'
#' @param x An object of class \code{sample_size_nbinom_summary}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.sample_size_nbinom_summary <- function(x, ...) {
  cat(strwrap(x, width = 80), sep = "\n")
  cat("\n")
  invisible(x)
}
