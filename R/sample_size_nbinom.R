#' Sample size calculation for negative binomial outcomes
#'
#' Computes the sample size (or power) for comparing two treatment groups
#' assuming negative binomial distributed event counts. When
#' `test_type = "wald"` (default), the formula uses a single variance
#' evaluated under the alternative, corresponding to Method 3 of
#' Zhu & Lakkis (2014) and the formulas of Friede & Schmidli (2010)
#' and Mutze et al. (2019). When `test_type = "score"`, separate null
#' and alternative variances are used (Farrington & Manning style),
#' which gives better Type I error control for the score test.
#'
#' @param lambda1 Event rate for group 1 (control), in events per unit time.
#' @param lambda2 Event rate for group 2 (treatment), in events per unit time.
#' @param rr0 Rate ratio under the null hypothesis
#'   (\eqn{\lambda_2 / \lambda_1}). Default is 1 (superiority).
#'   For non-inferiority, use a value > 1 (e.g., 1.1).
#'   For super-superiority, use a value < 1 (e.g., 0.8).
#' @param dispersion Dispersion parameter \eqn{k} such that
#'   \eqn{\mathrm{Var}(Y) = \mu + k\mu^2}.
#'   Equivalent to `1/size` in [stats::rnbinom()]. Can be a scalar (common
#'   dispersion) or a vector of length 2 (group-specific: control, treatment).
#' @param power Target power (\eqn{1 - \beta}). If `NULL`, power is computed
#'   for the given accrual rates (no sample size scaling). Default is 0.9.
#' @param alpha Significance level. Default is 0.025.
#' @param sided Number of sides for the test: 1 (one-sided) or 2 (two-sided).
#'   Default is 1.
#' @param ratio Allocation ratio \eqn{r = n_2/n_1}. Default is 1 (equal
#'   allocation).
#' @param accrual_rate Vector of accrual rates (patients per unit time) for
#'   each recruitment segment.
#' @param accrual_duration Vector of durations for each accrual segment.
#'   Must be the same length as `accrual_rate`.
#' @param trial_duration Total planned duration of the trial. If
#'   `trial_duration` is less than the sum of `accrual_duration`, accrual
#'   is truncated at `trial_duration`.
#' @param dropout_rate Dropout hazard rate. Can be:
#'   \itemize{
#'     \item A scalar (common constant rate for both groups). Default is 0.
#'     \item A vector of length 2 (group-specific constant rates: control,
#'       treatment).
#'     \item A data frame with columns `rate` and `duration` (and optionally
#'       `treatment`) defining piecewise constant dropout hazards. When a
#'       `treatment` column is present, use 1 for control and 2 for treatment.
#'       Without a `treatment` column, the same piecewise schedule applies to
#'       both groups. The last `duration` may be `Inf` to extend the final rate
#'       indefinitely.
#'   }
#' @param max_followup Maximum follow-up time for any patient. Default is
#'   `NULL` (infinite). Can be a vector of length 2 for group-specific caps.
#' @param test_type Type of test for which to size the study:
#'   `"wald"` (default) uses a single variance under the alternative;
#'   `"score"` uses separate null and alternative variances
#'   (\eqn{z_\alpha \sqrt{V_0} + z_\beta \sqrt{V_1}}).
#' @param event_gap Gap duration after each event during which no new events
#'   are counted (e.g., a recovery period). Default is `NULL` (no gap).
#'   When specified, the effective rate is reduced to
#'   \eqn{\lambda_{\mathrm{eff}} = \lambda / (1 + \lambda \cdot \mathrm{gap})}.
#'
#' @details
#' ## Sample size formula
#'
#' **Wald test** (`test_type = "wald"`):
#' \deqn{n_1 = \frac{(z_{\alpha/s} + z_\beta)^2 V_1}{(\theta - \theta_0)^2}}
#'
#' **Score test** (`test_type = "score"`):
#' \deqn{n_1 = \frac{(z_{\alpha/s} \sqrt{V_0} + z_\beta \sqrt{V_1})^2}{(\theta - \theta_0)^2}}
#'
#' where \eqn{\theta = \log(\lambda_2/\lambda_1)},
#' \eqn{\theta_0 = \log(\mathrm{rr}_0)}, and:
#' \deqn{V_1 = \left(\frac{1}{\mu_1} + k_1\right) + \frac{1}{r}\left(\frac{1}{\mu_2} + k_2\right)}
#' is the variance under \eqn{H_1}. Under \eqn{H_0} (pooled rate
#' \eqn{\lambda_0 = (\lambda_1 + r \lambda_2 \mathrm{rr}_0) / (1 + r)}):
#' \deqn{V_0 = \left(\frac{1}{\mu_0} + k_0\right)\left(1 + \frac{1}{r}\right)}
#' with \eqn{\mu_g = \lambda_g \bar{t}_g} the expected event count and
#' \eqn{\bar{t}_g} the average exposure for group \eqn{g}.
#'
#' ## Average exposure
#'
#' The average exposure \eqn{\bar{t}_g} accounts for piecewise accrual,
#' piecewise exponential dropout, and maximum follow-up truncation. With
#' piecewise constant dropout hazards \eqn{\delta_1, \delta_2, \ldots}
#' over successive intervals, the survival function is
#' \eqn{S(t) = \exp(-\sum_j \delta_j \ell_j)} where \eqn{\ell_j} is the
#' time spent in interval \eqn{j}. The expected exposure for a patient with
#' potential follow-up \eqn{u} is \eqn{m(u) = \int_0^u S(t)\,dt}, computed
#' as a sum of exponential integrals over each piece. For a single constant
#' rate \eqn{\delta > 0} this simplifies to
#' \eqn{m(u) = (1 - e^{-\delta u})/\delta}.
#' The overall average is a weighted mean across accrual segments.
#'
#' ## Variance inflation
#'
#' When follow-up times are variable, the dispersion is inflated by a factor
#' \eqn{Q_g = \mathrm{E}[t_g^2] / (\mathrm{E}[t_g])^2 \ge 1} (Zhu & Lakkis,
#' 2014) to account for the non-linear dependence of the NB variance on
#' exposure.
#'
#' ## Event gap correction (Jensen's inequality)
#'
#' When `event_gap` > 0, the naive effective rate
#' \eqn{\lambda / (1 + \lambda g)} overestimates the true population-level
#' effective rate because of subject-level heterogeneity (frailty).
#' In the Gamma-Poisson mixture, each subject's rate
#' \eqn{\Lambda_i \sim \mathrm{Gamma}(1/k, k\lambda)} is random.
#' Since \eqn{f(x) = x/(1+xg)} is concave, Jensen's inequality gives
#' \eqn{\mathrm{E}[f(\Lambda)] < f(\mathrm{E}[\Lambda])}.
#'
#' A second-order Taylor correction is applied:
#' \deqn{\lambda_{\mathrm{eff}} \approx \frac{\lambda}{1+\lambda g}
#'   \left(1 - \frac{k \lambda g}{(1+\lambda g)^2}\right)}
#' This uses \eqn{f''(\lambda) = -2g/(1+\lambda g)^3} and
#' \eqn{\mathrm{Var}(\Lambda) = k\lambda^2}.
#'
#' @return An object of class `sample_size_nbinom_result`, which is a list
#'   containing:
#' \describe{
#'   \item{inputs}{Named list of the original function arguments.}
#'   \item{n1}{Sample size for group 1 (control).}
#'   \item{n2}{Sample size for group 2 (treatment).}
#'   \item{n_total}{Total sample size (\eqn{n_1 + n_2}).}
#'   \item{alpha}{Significance level used.}
#'   \item{sided}{One-sided or two-sided test.}
#'   \item{power}{Power of the test.}
#'   \item{exposure}{Average calendar exposure \eqn{\bar{t}_g} (vector of
#'     length 2 for control and treatment).}
#'   \item{exposure_at_risk_n1}{Average at-risk exposure for group 1
#'     (adjusted for event gap).}
#'   \item{exposure_at_risk_n2}{Average at-risk exposure for group 2
#'     (adjusted for event gap).}
#'   \item{events_n1}{Expected number of events in group 1.}
#'   \item{events_n2}{Expected number of events in group 2.}
#'   \item{total_events}{Total expected number of events.}
#'   \item{variance}{Variance of the log rate ratio
#'     \eqn{\mathrm{Var}(\hat\theta)}.}
#'   \item{accrual_rate}{Accrual rate(s) used (possibly scaled to achieve
#'     target power).}
#'   \item{accrual_duration}{Accrual duration(s) used.}
#' }
#'
#' @references
#' Zhu, H., & Lakkis, H. (2014).
#' Sample size calculation for comparing two negative binomial rates.
#' _Statistics in Medicine_,
#' 33(3), 376--387. \doi{10.1002/sim.5947}
#'
#' Friede, T., & Schmidli, H. (2010).
#' Blinded sample size reestimation with negative binomial counts in
#' superiority and non-inferiority trials.
#' _Methods of Information in Medicine_,
#' 49(06), 618--624. \doi{10.3414/ME09-02-0060}
#'
#' Mutze, T., Glimm, E., Schmidli, H., & Friede, T. (2019).
#' Group sequential designs for negative binomial outcomes.
#' _Statistical Methods in Medical Research_,
#' 28(8), 2326--2347. \doi{10.1177/0962280218773115}
#'
#' @seealso
#' [compute_info_at_time()] for computing statistical information at a given
#' analysis time; [blinded_ssr()] for blinded sample size reestimation;
#' [gsNBCalendar()] for group sequential designs;
#' `vignette("sample-size-nbinom", package = "gsDesignNB")` for detailed
#' methodology.
#'
#' @importFrom stats pnorm qnorm
#'
#' @export
#'
#' @examples
#' # Basic sample size calculation
#' x <- sample_size_nbinom(
#'   lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
#'   accrual_rate = 10, accrual_duration = 20, trial_duration = 24
#' )
#' class(x)
#' summary(x)
#'
#' # With piecewise accrual
#' x2 <- sample_size_nbinom(
#'   lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
#'   accrual_rate = c(5, 10), accrual_duration = c(3, 3),
#'   trial_duration = 12
#' )
#' summary(x2)
#'
#' # Compute power for a fixed design (power = NULL)
#' sample_size_nbinom(
#'   lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = NULL,
#'   accrual_rate = 10, accrual_duration = 20, trial_duration = 24
#' )
sample_size_nbinom <- function(
  lambda1, lambda2, dispersion, power = NULL,
  alpha = 0.025, sided = 1, ratio = 1, rr0 = 1,
  accrual_rate, accrual_duration,
  trial_duration, dropout_rate = 0,
  max_followup = NULL, test_type = c("wald", "score"),
  event_gap = NULL
) {
  test_type <- match.arg(test_type)
  if (lambda1 <= 0 || lambda2 <= 0) {
    stop("Rates lambda1 and lambda2 must be positive.")
  }
  if (rr0 <= 0) {
    stop("Null hypothesis rate ratio rr0 must be positive.")
  }
  if (any(dispersion < 0)) {
    stop("Dispersion parameter must be non-negative.")
  }
  if (length(dispersion) == 1) {
    dispersion <- rep(dispersion, 2)
  } else if (length(dispersion) != 2) {
    stop("Dispersion must be a scalar or a vector of length 2.")
  }

  if (!is.null(power) && (power <= 0 || power >= 1)) {
    stop("Power must be between 0 and 1.")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("Alpha must be between 0 and 1.")
  }
  # Normalize dropout_rate to internal piecewise specification:
  # a list of 2 data.frames (one per group), each with 'rate' and 'duration'.
  if (is.data.frame(dropout_rate)) {
    if (!all(c("rate", "duration") %in% names(dropout_rate))) {
      stop("dropout_rate data frame must have 'rate' and 'duration' columns.")
    }
    if (any(dropout_rate$rate < 0)) {
      stop("Dropout rates must be non-negative.")
    }
    dropout_spec <- list()
    for (g in 1:2) {
      if ("treatment" %in% names(dropout_rate)) {
        sub <- dropout_rate[dropout_rate$treatment == g, , drop = FALSE]
        if (nrow(sub) == 0) {
          sub <- dropout_rate[is.na(dropout_rate$treatment), , drop = FALSE]
        }
      } else {
        sub <- dropout_rate
      }
      if (nrow(sub) == 0) {
        dropout_spec[[g]] <- data.frame(rate = 0, duration = Inf)
      } else {
        dropout_spec[[g]] <- data.frame(rate = sub$rate, duration = sub$duration)
      }
    }
  } else {
    if (any(dropout_rate < 0)) {
      stop("Dropout rate must be non-negative.")
    }
    if (length(dropout_rate) == 1) {
      dropout_rate <- rep(dropout_rate, 2)
    } else if (length(dropout_rate) != 2) {
      stop("Dropout rate must be a scalar, a vector of length 2, or a data frame.")
    }
    dropout_spec <- list(
      data.frame(rate = dropout_rate[1], duration = Inf),
      data.frame(rate = dropout_rate[2], duration = Inf)
    )
  }

  if (!is.null(max_followup)) {
    if (any(max_followup <= 0)) {
      stop("max_followup must be positive.")
    }
    if (length(max_followup) == 1) {
      max_followup <- rep(max_followup, 2)
    } else if (length(max_followup) != 2) {
      stop("max_followup must be a scalar or a vector of length 2.")
    }
  } else {
    max_followup <- c(Inf, Inf)
  }

  # Truncate accrual if trial_duration < sum(accrual_duration)
  cum_accrual_duration <- cumsum(accrual_duration)
  total_accrual_duration <- utils::tail(cum_accrual_duration, 1)
  
  if (trial_duration < total_accrual_duration) {
    # Find which period contains trial_duration
    idx <- which(cum_accrual_duration >= trial_duration)[1]
    
    # Truncate vectors
    accrual_rate <- accrual_rate[1:idx]
    accrual_duration <- accrual_duration[1:idx]
    
    # Adjust the last duration
    prev_dur <- if (idx == 1) 0 else cum_accrual_duration[idx - 1]
    accrual_duration[idx] <- trial_duration - prev_dur
  }

  power_input <- power
  inputs <- list(
    lambda1 = lambda1,
    lambda2 = lambda2,
    rr0 = rr0,
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
    test_type = test_type,
    event_gap = event_gap
  )

  # Determine mode: Calculate N or Calculate Power
  mode <- "solve_n"
  if (!is.null(power)) {
    # We are solving for N (scaling accrual) to meet this power
    mode <- "solve_n"
  } else {
    # Power is null, so we calculate power for the fixed accrual
    mode <- "solve_power"
  }

  # --- Piecewise exponential dropout helpers ---
  # Building-block integrals.
  # I1(delta, L) = integral_0^L exp(-delta*s) ds
  .I1 <- function(delta, L) {
    if (L <= 0) return(0)
    if (delta == 0) return(L)
    (1 - exp(-delta * L)) / delta
  }
  # I2(delta, L) = integral_0^L s * exp(-delta*s) ds
  .I2 <- function(delta, L) {
    if (L <= 0) return(0)
    if (delta == 0) return(L^2 / 2)
    (1 - exp(-delta * L) * (1 + delta * L)) / delta^2
  }
  # J1(delta, sa, sb) = integral_{sa}^{sb} I1(delta, s) ds
  .J1 <- function(delta, sa, sb) {
    if (sb <= sa) return(0)
    if (delta == 0) return((sb^2 - sa^2) / 2)
    (sb - sa) / delta - (exp(-delta * sa) - exp(-delta * sb)) / delta^2
  }
  # J2(delta, sa, sb) = integral_{sa}^{sb} I2(delta, s) ds
  .J2 <- function(delta, sa, sb) {
    if (sb <= sa) return(0)
    if (delta == 0) return((sb^3 - sa^3) / 6)
    ((sb - sa) - exp(-delta * sa) * (2 / delta + sa) +
      exp(-delta * sb) * (2 / delta + sb)) / delta^2
  }

  # Expected exposure m1(u) = integral_0^u S(t) dt
  pw_m1 <- function(u, dr_spec) {
    if (u <= 0) return(0)
    rates <- dr_spec$rate
    durations <- dr_spec$duration
    cum_start <- 0
    S_curr <- 1
    m1 <- 0
    for (j in seq_along(rates)) {
      piece_end <- cum_start + durations[j]
      L <- min(u, piece_end) - cum_start
      if (L <= 0) break
      m1 <- m1 + S_curr * .I1(rates[j], L)
      if (u <= piece_end) break
      S_curr <- S_curr * exp(-rates[j] * durations[j])
      cum_start <- piece_end
    }
    m1
  }

  # E[min(T,u)^2] = 2 * integral_0^u t * S(t) dt
  pw_m2 <- function(u, dr_spec) {
    if (u <= 0) return(0)
    rates <- dr_spec$rate
    durations <- dr_spec$duration
    cum_start <- 0
    S_curr <- 1
    m2 <- 0
    for (j in seq_along(rates)) {
      piece_end <- cum_start + durations[j]
      L <- min(u, piece_end) - cum_start
      if (L <= 0) break
      m2 <- m2 + 2 * S_curr * (cum_start * .I1(rates[j], L) + .I2(rates[j], L))
      if (u <= piece_end) break
      S_curr <- S_curr * exp(-rates[j] * durations[j])
      cum_start <- piece_end
    }
    m2
  }

  # Average of m1(u) over u uniformly distributed in [u_min, u_max]
  pw_avg_m1 <- function(u_min, u_max, dr_spec) {
    if (u_min >= u_max) return(0)
    rates <- dr_spec$rate
    durations <- dr_spec$duration
    cum_start <- 0
    S_curr <- 1
    m1_at_bp <- 0
    total <- 0
    for (j in seq_along(rates)) {
      piece_end <- cum_start + durations[j]
      a <- max(u_min, cum_start)
      b <- min(u_max, piece_end)
      if (a < b) {
        sa <- a - cum_start
        sb <- b - cum_start
        total <- total + m1_at_bp * (b - a) + S_curr * .J1(rates[j], sa, sb)
      }
      if (u_max <= piece_end) break
      m1_at_bp <- m1_at_bp + S_curr * .I1(rates[j], durations[j])
      S_curr <- S_curr * exp(-rates[j] * durations[j])
      cum_start <- piece_end
    }
    total / (u_max - u_min)
  }

  # Average of m2(u) over u uniformly distributed in [u_min, u_max]
  pw_avg_m2 <- function(u_min, u_max, dr_spec) {
    if (u_min >= u_max) return(0)
    rates <- dr_spec$rate
    durations <- dr_spec$duration
    cum_start <- 0
    S_curr <- 1
    m1_at_bp <- 0
    m2_at_bp <- 0
    total <- 0
    for (j in seq_along(rates)) {
      piece_end <- cum_start + durations[j]
      a <- max(u_min, cum_start)
      b <- min(u_max, piece_end)
      if (a < b) {
        sa <- a - cum_start
        sb <- b - cum_start
        total <- total + m2_at_bp * (b - a) +
          2 * S_curr * (cum_start * .J1(rates[j], sa, sb) + .J2(rates[j], sa, sb))
      }
      if (u_max <= piece_end) break
      L_j <- durations[j]
      m1_at_bp <- m1_at_bp + S_curr * .I1(rates[j], L_j)
      m2_at_bp <- m2_at_bp + 2 * S_curr * (cum_start * .I1(rates[j], L_j) + .I2(rates[j], L_j))
      S_curr <- S_curr * exp(-rates[j] * L_j)
      cum_start <- piece_end
    }
    total / (u_max - u_min)
  }

  # Calculate average exposure
  current_time <- 0
  total_n_accrual <- 0
  total_exposure_mass <- c(0, 0)
  total_exposure_sq_mass <- c(0, 0)

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

      # Calculate for each group (1 and 2)
      for (g in 1:2) {
        dr_spec <- dropout_spec[[g]]
        mf <- max_followup[g]

        avg_followup <- 0
        avg_followup_sq <- 0

        if (is.infinite(mf) || u_max <= mf) {
          # Case 1: No truncation by max_followup
          avg_followup <- pw_avg_m1(u_min, u_max, dr_spec)
          avg_followup_sq <- pw_avg_m2(u_min, u_max, dr_spec)
        } else if (u_min >= mf) {
          # Case 2: All truncated by max_followup
          avg_followup <- pw_m1(mf, dr_spec)
          avg_followup_sq <- pw_m2(mf, dr_spec)
        } else {
          # Case 3: Split at max_followup
          len_truncated <- u_max - mf
          len_not_truncated <- mf - u_min

          avg_1 <- pw_m1(mf, dr_spec)
          avg_sq_1 <- pw_m2(mf, dr_spec)
          avg_2 <- pw_avg_m1(u_min, mf, dr_spec)
          avg_sq_2 <- pw_avg_m2(u_min, mf, dr_spec)

          avg_followup <- (len_truncated * avg_1 + len_not_truncated * avg_2) / d
          avg_followup_sq <- (len_truncated * avg_sq_1 + len_not_truncated * avg_sq_2) / d
        }

        total_exposure_mass[g] <- total_exposure_mass[g] + n_seg * avg_followup
        total_exposure_sq_mass[g] <- total_exposure_sq_mass[g] + n_seg * avg_followup_sq
      }

      total_n_accrual <- total_n_accrual + n_seg
    }
    current_time <- current_time + d
  }

  if (total_n_accrual == 0) {
    stop("Accrual results in 0 patients.")
  }

  exposure_calendar <- total_exposure_mass / total_n_accrual
  exposure_sq_avg <- total_exposure_sq_mass / total_n_accrual

  # Calculate inflation factor Q for variance due to variable follow-up
  # Q = E[t^2] / (E[t])^2
  # If exposure is constant, Q = 1.
  Q_inflation <- rep(1, 2)
  for (g in 1:2) {
    if (exposure_calendar[g] > 0) {
      Q_inflation[g] <- exposure_sq_avg[g] / (exposure_calendar[g]^2)
    }
  }

  # Setup effective rates and exposures based on event_gap
  if (!is.null(event_gap) && !is.na(event_gap) && event_gap > 0) {
    # Second-order Taylor correction for Jensen's inequality bias.
    # With subject-level frailty Lambda ~ Gamma(1/k, k*lambda), the
    # population-level effective rate is E[Lambda / (1 + Lambda*gap)].
    # Since f(x) = x/(1+x*g) is concave, E[f(Lambda)] < f(E[Lambda])
    # by Jensen's inequality.
    #
    # Taylor expansion around E[Lambda] = lambda:
    #   E[f(Lambda)] ~ f(lambda) + f''(lambda) * Var(Lambda) / 2
    # where f''(x) = -2g / (1+xg)^3 and Var(Lambda) = k * lambda^2.
    #
    # This gives:
    #   lambda_eff ~ lambda/(1+lambda*g) * (1 - k*lambda*g / (1+lambda*g)^2)
    #
    # Note: dispersion has already been expanded to length 2. Use the
    # *base* dispersion (before Q inflation) for this correction since Q
    # addresses variable follow-up, not frailty.
    lambda1_eff <- lambda1 / (1 + lambda1 * event_gap) *
      (1 - dispersion[1] * lambda1 * event_gap / (1 + lambda1 * event_gap)^2)
    lambda2_eff <- lambda2 / (1 + lambda2 * event_gap) *
      (1 - dispersion[2] * lambda2 * event_gap / (1 + lambda2 * event_gap)^2)

    # Ensure effective rates remain positive (correction is small for
    # reasonable parameter combinations but could go negative in extreme cases)
    lambda1_eff <- max(lambda1_eff, 0)
    lambda2_eff <- max(lambda2_eff, 0)

    # Adjusted exposures for reporting (at-risk)
    exposure1_at_risk <- exposure_calendar[1] / (1 + lambda1 * event_gap)
    exposure2_at_risk <- exposure_calendar[2] / (1 + lambda2 * event_gap)
  } else {
    lambda1_eff <- lambda1
    lambda2_eff <- lambda2

    exposure1_at_risk <- exposure_calendar[1]
    exposure2_at_risk <- exposure_calendar[2]
  }

  mu1 <- lambda1_eff * exposure_calendar[1]
  mu2 <- lambda2_eff * exposure_calendar[2]

  # Apply inflation factor to dispersion
  k1 <- dispersion[1] * Q_inflation[1]
  k2 <- dispersion[2] * Q_inflation[2]

  z_alpha <- qnorm(1 - alpha / sided)

  n1_c <- 0
  n2_c <- 0
  n_total_c <- 0
  computed_accrual_rate <- NULL

  # Variance under H1 (alternative)
  V1 <- (1 / mu1 + k1) + (1 / ratio) * (1 / mu2 + k2)

  # Variance under H0 (null) — used only for score test
  # Pooled rate under H0: weighted average assuming equal exposure
  lambda0_eff <- (lambda1_eff + ratio * lambda2_eff * rr0) / (1 + ratio)
  mu0 <- lambda0_eff * exposure_calendar[1]  # use control exposure as reference
  k0 <- mean(dispersion) * mean(Q_inflation)  # pooled dispersion
  V0 <- (1 / mu0 + k0) * (1 + 1 / ratio)

  if (mode == "solve_n") {
    z_beta <- qnorm(power)

    if (test_type == "score") {
      num <- (z_alpha * sqrt(V0) + z_beta * sqrt(V1))^2
    } else {
      num <- (z_alpha + z_beta)^2 * V1
    }
    den <- (log(lambda1 * rr0 / lambda2))^2
    n1 <- num / den
    n2 <- n1 * ratio

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

    theta <- log(lambda1 * rr0 / lambda2)
    if (test_type == "score") {
      # Solve: z_alpha * sqrt(V0) + z_beta * sqrt(V1) = |theta| * sqrt(n1)
      z_beta <- (abs(theta) * sqrt(n1_c) - z_alpha * sqrt(V0)) / sqrt(V1)
    } else {
      z_beta <- sqrt(n1_c * theta^2 / V1) - z_alpha
    }

    power <- pnorm(z_beta)
  }

  variance <- (1 / mu1 + k1) / n1_c + (1 / mu2 + k2) / n2_c

  # Calculate expected events
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
      test_type = test_type,
      exposure = exposure_calendar,
      exposure_at_risk_n1 = exposure1_at_risk,
      exposure_at_risk_n2 = exposure2_at_risk,
      events_n1 = events_n1,
      events_n2 = events_n2,
      total_events = total_events,
      variance = variance,
      variance_null = V0 / n1_c * (1 + ratio) / ratio,
      accrual_rate = computed_accrual_rate,
      accrual_duration = accrual_duration
    )
  )
  class(result) <- c("sample_size_nbinom_result", "list")
  result
}


#' Print method for sample_size_nbinom_result objects
#'
#' Prints a concise summary of the sample size calculation results.
#'
#' @param x An object of class `sample_size_nbinom_result`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @examples
#' x <- sample_size_nbinom(
#'   lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
#'   accrual_rate = 10, accrual_duration = 20, trial_duration = 24
#' )
#' print(x)
#'
#' @export
print.sample_size_nbinom_result <- function(x, ...) {
  cat("Sample size for negative binomial outcome\n")
  cat("==========================================\n\n")

  cat(sprintf(
    "Sample size:     n1 = %d, n2 = %d, total = %d\n",
    x$n1, x$n2, x$n_total
  ))
  cat(sprintf(
    "Expected events: %.1f (n1: %.1f, n2: %.1f)\n",
    x$total_events, x$events_n1, x$events_n2
  ))
  cat(sprintf(
    "Power: %.0f%%, Alpha: %.3f (%d-sided)\n",
    x$power * 100, x$alpha, x$sided
  ))
  cat(sprintf(
    "Rates: control = %.4f, treatment = %.4f (RR = %.4f)\n",
    x$inputs$lambda1, x$inputs$lambda2,
    x$inputs$lambda2 / x$inputs$lambda1
  ))

  if (!is.null(x$inputs$rr0) && x$inputs$rr0 != 1) {
    cat(sprintf("Null hypothesis RR: %.4f\n", x$inputs$rr0))
  }

  # Handle dispersion display
  if (x$inputs$dispersion[1] == x$inputs$dispersion[2]) {
    disp_str <- sprintf("%.4f", x$inputs$dispersion[1])
  } else {
    disp_str <- sprintf("%.4f (n1), %.4f (n2)", x$inputs$dispersion[1], x$inputs$dispersion[2])
  }

  # Handle exposure display
  if (abs(x$exposure[1] - x$exposure[2]) < 1e-6) {
    exp_str <- sprintf("%.2f", x$exposure[1])
  } else {
    exp_str <- sprintf("%.2f (n1), %.2f (n2)", x$exposure[1], x$exposure[2])
  }

  cat(sprintf(
    "Dispersion: %s, Avg exposure (calendar): %s\n",
    disp_str, exp_str
  ))

  if (!is.null(x$inputs$event_gap) && x$inputs$event_gap > 0) {
    cat(sprintf(
      "Avg exposure (at-risk): n1 = %.2f, n2 = %.2f\n",
      x$exposure_at_risk_n1, x$exposure_at_risk_n2
    ))
    cat(sprintf("Event gap: %.2f\n", x$inputs$event_gap))
  }

  if (!is.null(x$inputs$dropout_rate)) {
    dr <- x$inputs$dropout_rate
    if (is.data.frame(dr)) {
      cat("Dropout rate: piecewise\n")
      if ("treatment" %in% names(dr)) {
        for (g in unique(dr$treatment)) {
          sub <- dr[dr$treatment == g, ]
          pieces <- paste0(sprintf("%.4f (%.1f)", sub$rate, sub$duration), collapse = ", ")
          cat(sprintf("  Group %s: %s\n", g, pieces))
        }
      } else {
        pieces <- paste0(sprintf("%.4f (%.1f)", dr$rate, dr$duration), collapse = ", ")
        cat(sprintf("  Both groups: %s\n", pieces))
      }
    } else if (any(dr > 0)) {
      if (dr[1] == dr[2]) {
        cat(sprintf("Dropout rate: %.4f\n", dr[1]))
      } else {
        cat(sprintf("Dropout rate: %.4f (n1), %.4f (n2)\n", dr[1], dr[2]))
      }
    }
  }

  cat(sprintf(
    "Accrual: %.1f, Trial duration: %.1f\n",
    sum(x$inputs$accrual_duration), x$inputs$trial_duration
  ))

  if (!is.null(x$inputs$max_followup)) {
    if (all(is.infinite(x$inputs$max_followup))) {
      # Do nothing if both are infinite (default)
    } else {
      cat(sprintf("Max follow-up: %.1f\n", x$inputs$max_followup[1]))
    }
  }

  invisible(x)
}


#' Summary for sample_size_nbinom_result objects
#'
#' Provides a textual summary of the sample size calculation for negative binomial
#' outcomes, similar to the summary for gsNB objects.
#'
#' @param object An object of class `sample_size_nbinom_result`.
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

  rr0_text <- ""
  if (!is.null(inputs$rr0) && inputs$rr0 != 1) {
    rr0_text <- sprintf("null hypothesis RR %.4f, ", inputs$rr0)
  }

  # Handle dispersion
  if (inputs$dispersion[1] == inputs$dispersion[2]) {
    disp_text <- sprintf("dispersion %.4f", inputs$dispersion[1])
  } else {
    disp_text <- sprintf("dispersion %.4f (n1) / %.4f (n2)", inputs$dispersion[1], inputs$dispersion[2])
  }

  # Handle exposure
  if (abs(object$exposure[1] - object$exposure[2]) < 1e-6) {
    exp_text <- sprintf("average exposure %.2f", object$exposure[1])
  } else {
    exp_text <- sprintf("average exposure %.2f (n1) / %.2f (n2)", object$exposure[1], object$exposure[2])
  }

  # Handle event gap in summary
  gap_text <- ""
  if (!is.null(inputs$event_gap) && inputs$event_gap > 0) {
    gap_text <- sprintf(
      " Event gap %.2f implies average at-risk exposure %.2f (n1) / %.2f (n2).",
      inputs$event_gap, object$exposure_at_risk_n1, object$exposure_at_risk_n2
    )
  }

  # Build the summary text
  summary_text <- sprintf(
    paste0(
      "Fixed sample size design for negative binomial outcome, ",
      "total sample size %d (n1=%d, n2=%d), ",
      "%.0f percent power, ",
      "%.1f percent (%d-sided) Type I error. ",
      "Control rate %.4f, treatment rate %.4f, ",
      "risk ratio %.4f, %s%s. ",
      "Accrual duration %.1f, trial duration %.1f, ",
      "%s.%s ",
      "Expected events %.1f. ",
      "Randomization ratio %.0f:1."
    ),
    object$n_total,
    object$n1,
    object$n2,
    object$power * 100,
    object$alpha * 100,
    inputs$sided,
    inputs$lambda1,
    inputs$lambda2,
    risk_ratio,
    rr0_text,
    disp_text,
    sum(inputs$accrual_duration),
    inputs$trial_duration,
    exp_text,
    gap_text,
    object$total_events,
    inputs$ratio
  )

  class(summary_text) <- "sample_size_nbinom_summary"
  summary_text
}


#' Print method for sample_size_nbinom_summary objects
#'
#' @param x An object of class `sample_size_nbinom_summary`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @examples
#' x <- sample_size_nbinom(
#'   lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
#'   accrual_rate = 10, accrual_duration = 20, trial_duration = 24
#' )
#' s <- summary(x)
#' print(s)
#'
#' @export
print.sample_size_nbinom_summary <- function(x, ...) {
  cat(strwrap(x, width = 80), sep = "\n")
  cat("\n")
  invisible(x)
}
