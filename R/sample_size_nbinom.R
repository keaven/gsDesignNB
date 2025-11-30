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
#' @param exposure Duration of exposure. Default is 1. Ignored if \code{accrual_rate}
#'   and related parameters are provided.
#' @param ratio Allocation ratio n2/n1. Default is 1.
#' @param accrual_rate Vector of accrual rates (patients per unit time). If provided,
#'   \code{accrual_duration} and \code{trial_duration} must also be provided.
#' @param accrual_duration Vector of durations for each accrual rate. Must be same length
#'   as \code{accrual_rate}.
#' @param trial_duration Total planned duration of the trial.
#' @param dropout_rate Dropout rate (hazard rate). Default is 0.
#' @param max_followup Maximum follow-up time for any patient. Default is NULL (infinite).
#' @param method Method for sample size calculation. "zhu" for Zhu and Lakkis (2014),
#'   "friede" for Friede and Schmidli (2010) / Mütze et al. (2018).
#'
#' @return A list containing:
#' \describe{
#'   \item{n1}{Sample size for group 1}
#'   \item{n2}{Sample size for group 2}
#'   \item{n_total}{Total sample size}
#'   \item{exposure}{Average exposure time used in calculation}
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
#' @export
#'
#' @examples
#' # Calculate sample size for lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1
#' sample_size_nbinom(lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8)
#'
#' # With piecewise accrual
#' # 5 patients/month for 3 months, then 10 patients/month for 3 months
#' # Trial ends at month 12.
#' sample_size_nbinom(
#'   lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
#'   accrual_rate = c(5, 10), accrual_duration = c(3, 3),
#'   trial_duration = 12
#' )
#'
#' @importFrom stats pnorm qnorm
sample_size_nbinom <- function(lambda1, lambda2, dispersion, power = NULL,
                               alpha = 0.025, sided = 1, exposure = NULL, ratio = 1,
                               accrual_rate = NULL, accrual_duration = NULL,
                               trial_duration = NULL, dropout_rate = 0,
                               max_followup = NULL, method = "zhu") {
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

  # Determine mode: Calculate N or Calculate Power
  mode <- "solve_n"
  if (is.null(power)) {
    if (!is.null(accrual_rate)) {
      mode <- "solve_power"
    } else {
      power <- 0.9 # Default if accrual not provided
    }
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

  if (!is.null(accrual_rate)) {
    # Handle variable accrual to calculate average exposure
    if (is.null(accrual_duration) || is.null(trial_duration)) {
      stop("If accrual_rate is provided, accrual_duration and trial_duration must also be provided.")
    }
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

    exposure <- total_exposure_mass / total_n_accrual
  } else {
    # No accrual info provided
    if (is.null(exposure)) {
      exposure <- 1 # Default if not provided
    }

    # Adjust exposure if max_followup is provided (act as a cap)
    if (!is.null(max_followup) && !is.infinite(max_followup)) {
      exposure <- min(exposure, max_followup)
    }

    # If dropout_rate > 0, adjust fixed exposure
    if (dropout_rate > 0) {
      exposure <- (1 - exp(-dropout_rate * exposure)) / dropout_rate
    }
  }

  mu1 <- lambda1 * exposure
  mu2 <- lambda2 * exposure
  k <- dispersion

  z_alpha <- qnorm(1 - alpha / sided)

  n1_c <- 0
  n2_c <- 0
  n_total_c <- 0
  computed_accrual_rate <- NULL

  if (mode == "solve_n") {
    z_beta <- qnorm(power)

    if (method == "zhu") {
      num <- (z_alpha + z_beta)^2 * ((1 / mu1 + k) + (1 / ratio) * (1 / mu2 + k))
      den <- (log(mu1 / mu2))^2
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
      stop("Unknown method. Choose 'zhu' or 'friede'.")
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

    if (method == "zhu") {
      # z_beta = sqrt( n1 * (log(mu1/mu2))^2 / V ) - z_alpha
      V <- (1 / mu1 + k) + (1 / ratio) * (1 / mu2 + k)
      z_beta <- sqrt(n1_c * (log(mu1 / mu2))^2 / V) - z_alpha
    } else if (method == "friede") {
      # z_beta = sqrt( n_total * (log(lambda1/lambda2))^2 / V_bar ) - z_alpha
      p1 <- 1 / (1 + ratio)
      p2 <- ratio / (1 + ratio)
      V <- (1 / mu1 + k) / p1 + (1 / mu2 + k) / p2
      z_beta <- sqrt(n_total_c * (log(lambda1 / lambda2))^2 / V) - z_alpha
    } else {
      stop("Unknown method. Choose 'zhu' or 'friede'.")
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

  return(list(
    n1 = n1_c,
    n2 = n2_c,
    n_total = n_total_c,
    alpha = alpha,
    sided = sided,
    power = power,
    exposure = exposure,
    events_n1 = events_n1,
    events_n2 = events_n2,
    total_events = total_events,
    variance = variance,
    accrual_rate = computed_accrual_rate,
    accrual_duration = accrual_duration
  ))
}
