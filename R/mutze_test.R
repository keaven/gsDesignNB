#' Wald test for treatment effect using negative binomial model (Mutze et al.)
#'
#' Fits a negative binomial (or Poisson) log-rate model to the aggregated
#' subject-level data produced by [cut_data_by_date()]. The method matches the
#' Wald test described by Mutze et al. (2019) for comparing treatment arms with
#' recurrent event outcomes. When the maximum likelihood negative binomial fit
#' is unreliable, the test automatically switches to one of two statistically
#' sensible fallbacks: a Poisson Wald test when the data are essentially
#' Poisson, or a method-of-moments (MoM) variance estimate plugged into the
#' same negative binomial information formula when the data are extremely
#' overdispersed or the ML fit fails to converge. The latter avoids the
#' anti-conservative behaviour of a blind Poisson fallback under genuine
#' overdispersion.
#'
#' @param data A data frame with at least the columns `treatment`, `events`, and
#'   `tte` (follow-up time). Typically output from [cut_data_by_date()].
#' @param method Type of model to fit: "nb" (default) uses a negative binomial
#'   GLM via [MASS::glm.nb()], "poisson" fits a Poisson GLM.
#' @param conf_level Confidence level for the rate ratio interval. Default 0.95.
#' @param sided Number of sides for the test: 1 (default) or 2.
#' @param poisson_threshold Upper threshold (in units of `fit$theta`, the
#'   `MASS::glm.nb()` shape parameter \eqn{\theta_{\text{NB}} = 1/k}) above
#'   which the data are treated as essentially Poisson and the function falls
#'   back to a Poisson Wald test. Default is 50, corresponding to
#'   \eqn{\hat{k} < 0.02}, by which point NB and Poisson Wald standard errors
#'   are numerically indistinguishable at typical trial sample sizes.
#' @param mom_threshold Lower threshold on `fit$theta` below which the NB ML
#'   fit is considered unreliable (extreme overdispersion). When triggered, or
#'   when `glm.nb()` fails to converge, the function falls back to a
#'   method-of-moments NB Wald test: rates and dispersion are re-estimated via
#'   [estimate_nb_mom()] and the Wald variance is computed from the Fisher
#'   information formula \eqn{\mathcal{I} = 1/(1/W_1 + 1/W_2)} with
#'   \eqn{W_g = \sum_i \mu_{g,i}/(1 + \hat{k}\mu_{g,i})}. Default is 20,
#'   corresponding to \eqn{\hat{k} > 20}. This avoids the anti-conservative
#'   variance of a Poisson fallback when the data are truly overdispersed.
#'
#' @return An object of class `mutze_test` containing the fitted model summary with elements:
#'   * `method`: A string indicating the test method used.
#'   * `estimate`: log rate ratio (experimental vs control).
#'   * `se`: standard error for the log rate ratio.
#'   * `z`: Wald statistic.
#'   * `p_value`: one-sided or two-sided p-value.
#'   * `rate_ratio`: estimated rate ratio and its confidence interval.
#'   * `dispersion`: estimated dispersion. For consistency this is reported on
#'     the `MASS::glm.nb()` scale (\eqn{\theta_{\text{NB}} = 1/k}), regardless
#'     of whether ML, Poisson fallback, or MoM fallback was used. `Inf`
#'     indicates Poisson (\eqn{k=0}).
#'   * `group_summary`: observed subjects/events/exposure per treatment.
#'   * `fallback`: character label describing which fit path was used
#'     (`"ml"`, `"poisson"`, or `"mom"`).
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
mutze_test <- function(data, method = c("nb", "poisson"), conf_level = 0.95, sided = 1,
                       poisson_threshold = 50, mom_threshold = 20) {
  method <- match.arg(method)
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
