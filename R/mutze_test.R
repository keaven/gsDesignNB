#' Wald or score test for treatment effect using negative binomial model
#'
#' Fits a negative binomial (or Poisson) log-rate model to the aggregated
#' subject-level data produced by [cut_data_by_date()]. With
#' `test_type = "wald"` (default), the method matches the
#' Wald test described by Mutze et al. (2019). With `test_type = "score"`,
#' the function fits only the null (no treatment effect) model and computes
#' the score statistic, which evaluates all quantities under \eqn{H_0} and
#' avoids the finite-sample anti-conservatism of the Wald test.
#'
#' When the maximum likelihood negative binomial fit
#' is unreliable, the test automatically switches to one of two statistically
#' sensible fallbacks: a Poisson test when the data are essentially
#' Poisson, or a method-of-moments (MoM) variance estimate plugged into the
#' same negative binomial information formula when the data are extremely
#' overdispersed or the ML fit fails to converge.
#'
#' @param data A data frame with at least the columns `treatment`, `events`, and
#'   `tte` (follow-up time). Typically output from [cut_data_by_date()].
#' @param method Type of model to fit: "nb" (default) uses a negative binomial
#'   GLM via [MASS::glm.nb()], "poisson" fits a Poisson GLM.
#' @param test_type Type of test statistic: `"wald"` (default) or `"score"`.
#'   The Wald test estimates the treatment effect under the alternative and
#'   divides by its standard error. The score test fits only the null model
#'   and evaluates the derivative of the log-likelihood at \eqn{\theta = 0},
#'   avoiding estimation under the alternative. The score test typically has
#'   better finite-sample Type I error control and is faster because it only
#'   fits a one-parameter null model.
#' @param conf_level Confidence level for the rate ratio interval. Default 0.95.
#' @param sided Number of sides for the test: 1 (default) or 2.
#' @param poisson_threshold Upper threshold (in units of `fit$theta`, the
#'   `MASS::glm.nb()` shape parameter \eqn{\theta_{\text{NB}} = 1/k}) above
#'   which the data are treated as essentially Poisson. Default is 50,
#'   corresponding to \eqn{\hat{k} < 0.02}.
#' @param mom_threshold Lower threshold on `fit$theta` below which the NB ML
#'   fit is considered unreliable (extreme overdispersion). Default is 20,
#'   corresponding to \eqn{\hat{k} > 20}.
#'
#' @return An object of class `mutze_test` containing:
#'   * `method`: A string indicating the test method used.
#'   * `estimate`: log rate ratio (experimental vs control). For
#'     `test_type = "score"`, this is a plug-in estimate.
#'   * `se`: standard error for the log rate ratio.
#'   * `z`: test statistic (Wald or score).
#'   * `p_value`: one-sided or two-sided p-value.
#'   * `rate_ratio`: estimated rate ratio and its confidence interval.
#'   * `dispersion`: estimated dispersion on the \eqn{\theta = 1/k} scale.
#'   * `group_summary`: observed subjects/events/exposure per treatment.
#'   * `fallback`: character label (`"ml"`, `"poisson"`, or `"mom"`).
#'   * `test_type`: character label (`"wald"` or `"score"`).
#'
#' @importFrom stats pnorm poisson qnorm
#' @importFrom utils tail
#'
#' @export
#'
#' @examples
#' enroll_rate <- data.frame(rate = 20 / (5 / 12), duration = 5 / 12)
#' fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3))
#' dropout_rate <- data.frame(
#'   treatment = c("Control", "Experimental"),
#'   rate = c(0.1, 0.05), duration = c(100, 100)
#' )
#' sim <- nb_sim(enroll_rate, fail_rate, dropout_rate, max_followup = 2, n = 40)
#' cut <- cut_data_by_date(sim, cut_date = 1.5)
#' mutze_test(cut)
#' mutze_test(cut, test_type = "score")
mutze_test <- function(data, method = c("nb", "poisson"), test_type = c("wald", "score"),
                       conf_level = 0.95, sided = 1,
                       poisson_threshold = 50, mom_threshold = 20) {
  method <- match.arg(method)
  test_type <- match.arg(test_type)
  df <- as.data.frame(data)
  required <- c("treatment", "events", "tte")
  if (!all(required %in% names(df))) {
    stop("Data must contain columns: ", paste(required, collapse = ", "), call. = FALSE)
  }
  if (any(df$tte <= 0)) {
    df <- df[df$tte > 0, , drop = FALSE]
  }
  if (nrow(df) == 0) {
    stop("No rows with positive follow-up time available.", call. = FALSE)
  }
  df$treatment <- droplevels(factor(df$treatment))
  if (nlevels(df$treatment) != 2) {
    stop("mutze_test currently supports exactly two treatment groups.", call. = FALSE)
  }
  if (poisson_threshold <= 0 || mom_threshold <= 0) {
    stop("poisson_threshold and mom_threshold must be positive.", call. = FALSE)
  }

  # Dispatch to score test if requested
  if (test_type == "score") {
    return(.mutze_score_test(df, method = method, conf_level = conf_level,
                            sided = sided, poisson_threshold = poisson_threshold,
                            mom_threshold = mom_threshold))
  }

  fallback <- NA_character_
  est <- NA_real_
  se <- NA_real_
  dispersion <- NA_real_
  method_label <- NA_character_
  fit <- NULL

  if (method == "nb") {
    fit <- tryCatch(
      suppressWarnings(MASS::glm.nb(events ~ treatment + offset(log(tte)), data = df)),
      error = function(e) NULL
    )
    ml_ok <- !is.null(fit) && isTRUE(fit$converged) && !is.na(fit$theta) &&
      is.finite(fit$theta) && fit$theta > 0

    near_poisson <- ml_ok && fit$theta > poisson_threshold
    extreme_overdisp <- ml_ok && fit$theta < 1 / mom_threshold

    if (ml_ok && !near_poisson && !extreme_overdisp) {
      model_summary <- summary(fit)
      coef_name <- tail(rownames(model_summary$coefficients), 1)
      est <- model_summary$coefficients[coef_name, "Estimate"]
      se <- model_summary$coefficients[coef_name, "Std. Error"]
      dispersion <- fit$theta
      method_label <- "Negative binomial Wald"
      fallback <- "ml"
    } else if (near_poisson || (!ml_ok && !.mom_prefers_overdispersion(df))) {
      # Poisson fallback: either ML produced a near-Poisson theta, or ML failed
      # and the data do not look overdispersed by MoM (in which case Poisson is
      # statistically correct).
      fit <- stats::glm(events ~ treatment + offset(log(tte)), data = df, family = poisson())
      model_summary <- summary(fit)
      coef_name <- tail(rownames(model_summary$coefficients), 1)
      est <- model_summary$coefficients[coef_name, "Estimate"]
      se <- model_summary$coefficients[coef_name, "Std. Error"]
      dispersion <- Inf
      method_label <- if (near_poisson) {
        "Poisson Wald (fallback, near-Poisson ML)"
      } else {
        "Poisson Wald (fallback, ML non-convergent, MoM k=0)"
      }
      fallback <- "poisson"
    } else {
      # Method-of-moments fallback: extreme overdispersion or ML failure with
      # MoM-estimated k > 0.
      mom_res <- .mutze_mom_wald(df)
      est <- mom_res$estimate
      se <- mom_res$se
      dispersion <- if (mom_res$k > 0) 1 / mom_res$k else Inf
      method_label <- if (ml_ok && extreme_overdisp) {
        "Negative binomial Wald (MoM fallback, extreme overdispersion)"
      } else {
        "Negative binomial Wald (MoM fallback, ML non-convergent)"
      }
      fallback <- "mom"
      fit <- mom_res$fit_info
    }
  } else {
    fit <- stats::glm(events ~ treatment + offset(log(tte)), data = df, family = poisson())
    model_summary <- summary(fit)
    coef_name <- tail(rownames(model_summary$coefficients), 1)
    est <- model_summary$coefficients[coef_name, "Estimate"]
    se <- model_summary$coefficients[coef_name, "Std. Error"]
    dispersion <- Inf
    method_label <- "Poisson Wald"
    fallback <- "poisson"
  }

  z <- est / se
  if (sided == 1) {
    pval <- stats::pnorm(z)
  } else {
    pval <- 2 * stats::pnorm(-abs(z))
  }
  alpha <- 1 - conf_level
  crit <- stats::qnorm(1 - alpha / 2)
  ci_log <- est + c(-1, 1) * crit * se
  rr <- exp(est)
  rr_ci <- exp(ci_log)

  group_summary <- data.table::as.data.table(df)[
    , .(subjects = .N, events = sum(events), exposure = sum(tte)),
    by = treatment
  ]
  group_summary <- as.data.frame(group_summary)

  res <- list(
    method = method_label,
    estimate = est,
    se = se,
    z = z,
    p_value = pval,
    sided = sided,
    rate_ratio = rr,
    conf_int = rr_ci,
    conf_level = conf_level,
    dispersion = dispersion,
    fallback = fallback,
    test_type = "wald",
    model = fit,
    group_summary = group_summary
  )
  class(res) <- "mutze_test"
  res
}

# Compute MoM-based point estimate, SE, and dispersion for the log rate ratio.
# Returns list(estimate, se, k, fit_info) where fit_info is a small list of
# estimates (not a fitted glm object) preserved in the result.
#' @noRd
.mutze_mom_wald <- function(df) {
  mom <- estimate_nb_mom(df, group = "treatment")
  lev <- levels(df$treatment)
  lam <- mom$lambda
  # estimate_nb_mom returns a named numeric indexed by factor levels
  lam1 <- as.numeric(lam[lev[1]])
  lam2 <- as.numeric(lam[lev[2]])
  k_hat <- as.numeric(mom$dispersion)

  if (!is.finite(lam1) || !is.finite(lam2) || lam1 <= 0 || lam2 <= 0) {
    stop("Method-of-moments estimates are not positive in both arms; ",
         "the NB Wald test cannot be formed on this interim.", call. = FALSE)
  }

  t1 <- df$tte[df$treatment == lev[1]]
  t2 <- df$tte[df$treatment == lev[2]]
  mu1 <- lam1 * t1
  mu2 <- lam2 * t2
  w1 <- sum(mu1 / (1 + k_hat * mu1))
  w2 <- sum(mu2 / (1 + k_hat * mu2))
  if (!is.finite(w1) || !is.finite(w2) || w1 <= 0 || w2 <= 0) {
    stop("Method-of-moments Fisher information is non-positive; ",
         "the NB Wald test cannot be formed on this interim.", call. = FALSE)
  }
  se <- sqrt(1 / w1 + 1 / w2)
  est <- log(lam2) - log(lam1)
  list(
    estimate = est,
    se = se,
    k = k_hat,
    fit_info = list(method = "mom", lambda1 = lam1, lambda2 = lam2, dispersion_k = k_hat)
  )
}

# Score test implementation.
# Fits only the null model (no treatment effect) and computes the score
# statistic U / sqrt(I_0) where U is the score for the treatment coefficient
# evaluated at the restricted MLE and I_0 is the corresponding Fisher
# information.
#' @noRd
.mutze_score_test <- function(df, method, conf_level, sided,
                              poisson_threshold, mom_threshold) {
  lev <- levels(df$treatment)
  is_trt <- df$treatment == lev[2]  # treatment indicator (second level)

  fallback <- NA_character_
  dispersion <- NA_real_
  method_label <- NA_character_
  fit0 <- NULL
  k_hat <- 0
  lambda0 <- NA_real_

  if (method == "nb") {
    # Fit null NB model (intercept + offset only)
    fit0 <- tryCatch(
      suppressWarnings(MASS::glm.nb(events ~ 1 + offset(log(tte)), data = df)),
      error = function(e) NULL
    )
    ml_ok <- !is.null(fit0) && isTRUE(fit0$converged) && !is.na(fit0$theta) &&
      is.finite(fit0$theta) && fit0$theta > 0

    near_poisson <- ml_ok && fit0$theta > poisson_threshold
    extreme_overdisp <- ml_ok && fit0$theta < 1 / mom_threshold

    if (ml_ok && !near_poisson && !extreme_overdisp) {
      k_hat <- 1 / fit0$theta
      lambda0 <- exp(coef(fit0)[1])
      dispersion <- fit0$theta
      method_label <- "Negative binomial score"
      fallback <- "ml"
    } else if (near_poisson || (!ml_ok && !.mom_prefers_overdispersion(df))) {
      # Poisson score test
      fit0 <- stats::glm(events ~ 1 + offset(log(tte)), data = df, family = poisson())
      k_hat <- 0
      lambda0 <- exp(coef(fit0)[1])
      dispersion <- Inf
      method_label <- if (isTRUE(near_poisson)) {
        "Poisson score (fallback, near-Poisson ML)"
      } else {
        "Poisson score (fallback, ML non-convergent, MoM k=0)"
      }
      fallback <- "poisson"
    } else {
      # MoM fallback: use pooled MoM estimates under null
      mom <- tryCatch(estimate_nb_mom(df), error = function(e) NULL)
      if (is.null(mom) || is.na(mom$dispersion)) {
        stop("Method-of-moments estimation failed under null model.", call. = FALSE)
      }
      k_hat <- max(0, as.numeric(mom$dispersion))
      lambda0 <- sum(df$events) / sum(df$tte)
      dispersion <- if (k_hat > 0) 1 / k_hat else Inf
      method_label <- if (ml_ok && extreme_overdisp) {
        "Negative binomial score (MoM fallback, extreme overdispersion)"
      } else {
        "Negative binomial score (MoM fallback, ML non-convergent)"
      }
      fallback <- "mom"
    }
  } else {
    # Forced Poisson
    fit0 <- stats::glm(events ~ 1 + offset(log(tte)), data = df, family = poisson())
    k_hat <- 0
    lambda0 <- exp(coef(fit0)[1])
    dispersion <- Inf
    method_label <- "Poisson score"
    fallback <- "poisson"
  }

  # Compute score statistic under H0
  # Edge case: if total events = 0, lambda0 = 0, mu0 = 0, U = 0, I0 = 0.
  # No events means no evidence, so z = 0.
  total_events <- sum(df$events)
  if (total_events == 0 || !is.finite(lambda0) || lambda0 <= 0) {
    z <- 0
    se <- Inf
    pval <- if (sided == 1) 0.5 else 1
    est <- 0
    rr <- 1
  } else {
    mu0 <- lambda0 * df$tte
    score_vec <- (df$events - mu0) / (1 + k_hat * mu0)
    U <- sum(score_vec[is_trt])

    # Fisher information under H0 for the treatment coefficient
    w1 <- sum(mu0[!is_trt] / (1 + k_hat * mu0[!is_trt]))
    w2 <- sum(mu0[is_trt] / (1 + k_hat * mu0[is_trt]))
    if (!is.finite(w1) || !is.finite(w2) || w1 <= 0 || w2 <= 0) {
      stop("Fisher information under null is non-positive.", call. = FALSE)
    }
    I0 <- w1 * w2 / (w1 + w2)
    z <- U / sqrt(I0)
    se <- 1 / sqrt(I0)

    if (sided == 1) {
      pval <- stats::pnorm(z)
    } else {
      pval <- 2 * stats::pnorm(-abs(z))
    }

    # Plug-in point estimate (crude log rate ratio)
    # If one group has 0 events, est = +/-Inf (legitimate: infinite evidence
    # for direction, though the score z remains finite).
    y1 <- sum(df$events[!is_trt])
    t1 <- sum(df$tte[!is_trt])
    y2 <- sum(df$events[is_trt])
    t2 <- sum(df$tte[is_trt])
    rate1 <- y1 / t1
    rate2 <- y2 / t2
    if (rate1 > 0 && rate2 > 0) {
      est <- log(rate2) - log(rate1)
    } else if (rate1 == 0 && rate2 == 0) {
      est <- 0
    } else if (rate2 == 0) {
      est <- -Inf
    } else {
      est <- Inf
    }
    rr <- exp(est)
  }

  alpha <- 1 - conf_level
  crit <- stats::qnorm(1 - alpha / 2)
  ci_log <- est + c(-1, 1) * crit * se
  rr_ci <- exp(ci_log)

  group_summary <- data.table::as.data.table(df)[
    , .(subjects = .N, events = sum(events), exposure = sum(tte)),
    by = treatment
  ]
  group_summary <- as.data.frame(group_summary)

  res <- list(
    method = method_label,
    estimate = est,
    se = se,
    z = z,
    p_value = pval,
    sided = sided,
    rate_ratio = rr,
    conf_int = rr_ci,
    conf_level = conf_level,
    dispersion = dispersion,
    fallback = fallback,
    test_type = "score",
    model = fit0,
    group_summary = group_summary
  )
  class(res) <- "mutze_test"
  res
}

# Check whether the data look overdispersed under a MoM pilot.
# Used only as a tiebreaker when glm.nb() fails to converge.
#' @noRd
.mom_prefers_overdispersion <- function(df) {
  mom <- tryCatch(estimate_nb_mom(df, group = "treatment"), error = function(e) NULL)
  if (is.null(mom) || is.na(mom$dispersion)) return(FALSE)
  isTRUE(mom$dispersion > 0)
}

#' @describeIn mutze_test Print method for `mutze_test` objects.
#' @param x An object of class `mutze_test`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object.
#' @export
print.mutze_test <- function(x, ...) {
  cat("Mutze Test Results\n")
  cat("==================\n\n")
  cat(paste("Method:    ", x$method, "\n"))
  cat(sprintf("Estimate:   %.4f\n", x$estimate))
  cat(sprintf("SE:         %.4f\n", x$se))
  cat(sprintf("Z:          %.4f\n", x$z))
  cat(sprintf("p-value:    %s\n", format.pval(x$p_value, digits = 4)))
  cat(sprintf("Rate Ratio: %.4f\n", x$rate_ratio))
  cat(sprintf("CI (%.0f%%):  [%.4f, %.4f]\n", x$conf_level * 100, x$conf_int[1], x$conf_int[2]))
  cat(sprintf("Dispersion: %.4f\n\n", x$dispersion))

  cat("Group Summary:\n")
  gs <- x$group_summary
  gs$subjects <- as.integer(gs$subjects)
  gs$events <- as.integer(gs$events)
  print(gs, row.names = FALSE)
  invisible(x)
}
