#' Wald test for treatment effect using negative binomial model (Mutze et al.)
#'
#' Fits a negative binomial (or Poisson) log-rate model to the aggregated
#' subject-level data produced by [cut_data_by_date()]. The method matches the
#' Wald test described by Mutze et al. (2019) for comparing treatment arms with
#' recurrent event outcomes.
#'
#' @param data A data frame with at least the columns `treatment`, `events`, and
#'   `tte` (follow-up time). Typically output from [cut_data_by_date()].
#' @param method Type of model to fit: "nb" (default) uses a negative binomial
#'   GLM via [MASS::glm.nb()], "poisson" fits a Poisson GLM.
#' @param conf_level Confidence level for the rate ratio interval. Default 0.95.
#' @param sided Number of sides for the test: 1 (default) or 2.
#' @param poisson_threshold When `method = "nb"`, the model falls back to
#'   Poisson regression if theta (the NB shape parameter) is outside the range
#'   `[1/poisson_threshold, poisson_threshold]`. Very large theta indicates
#'   near-Poisson data, while very small theta indicates extreme overdispersion
#'   with unstable estimates. Default is 1000.
#'
#' @return An object of class `mutze_test` containing the fitted model summary with elements:
#'   * `method`: A string indicating the test method used.
#'   * `estimate`: log rate ratio (experimental vs control).
#'   * `se`: standard error for the log rate ratio.
#'   * `z`: Wald statistic.
#'   * `p_value`: one-sided or two-sided p-value.
#'   * `rate_ratio`: estimated rate ratio and its confidence interval.
#'   * `dispersion`: estimated dispersion (theta) when `method = "nb"`.
#'   * `group_summary`: observed subjects/events/exposure per treatment.
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
                       poisson_threshold = 1000) {
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
  if (method == "nb") {
    fit <- tryCatch(
      suppressWarnings(MASS::glm.nb(events ~ treatment + offset(log(tte)), data = df)),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      if (!isTRUE(fit$converged) || is.na(fit$theta)) {
        fit <- NULL
      }
    }
    # If theta is very large (dispersion near 0), NB is essentially Poisson
    # If theta is very small (extreme overdispersion), the model is unstable
    # Fall back to Poisson for numerical stability in both cases
    if (!is.null(fit) && (fit$theta > poisson_threshold || fit$theta < 1/poisson_threshold)) {
      fit <- NULL
    }
    if (is.null(fit)) {
      fit <- stats::glm(events ~ treatment + offset(log(tte)), data = df, family = poisson())
      dispersion <- Inf
      method_label <- "Poisson Wald (fallback)"
    } else {
      dispersion <- fit$theta
      method_label <- "Negative binomial Wald"
    }
  } else {
    fit <- stats::glm(events ~ treatment + offset(log(tte)), data = df, family = poisson())
    dispersion <- Inf
    method_label <- "Poisson Wald"
  }
  model_summary <- summary(fit)
  # Identify coefficient comparing treatment levels
  coef_name <- tail(rownames(model_summary$coefficients), 1)
  est <- model_summary$coefficients[coef_name, "Estimate"]
  se <- model_summary$coefficients[coef_name, "Std. Error"]
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
    model = fit,
    group_summary = group_summary
  )
  class(res) <- "mutze_test"
  res
}

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
