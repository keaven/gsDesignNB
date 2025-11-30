#' Wald test for treatment effect using Negative Binomial model (Mütze et al.)
#'
#' Fits a Negative Binomial (or Poisson) log-rate model to the aggregated
#' subject-level data produced by [cut_data_by_date()]. The method matches the
#' Wald test described by Mütze et al. (2018) for comparing treatment arms with
#' recurrent event outcomes.
#'
#' @param data A data frame with at least the columns `treatment`, `events`, and
#'   `tte` (follow-up time). Typically output from [cut_data_by_date()].
#' @param method Type of model to fit: "nb" (default) uses a Negative Binomial
#'   GLM via [MASS::glm.nb()], "poisson" fits a Poisson GLM.
#' @param conf_level Confidence level for the rate ratio interval. Default 0.95.
#'
#' @return A list containing the fitted model summary with elements:
#'   * `estimate`: log rate ratio (experimental vs control).
#'   * `se`: standard error for the log rate ratio.
#'   * `z`: Wald statistic.
#'   * `p_value`: two-sided p-value.
#'   * `rate_ratio`: estimated rate ratio and its confidence interval.
#'   * `dispersion`: estimated dispersion (theta) when `method = "nb"`.
#'   * `group_summary`: observed subjects/events/exposure per treatment.
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
#'
#' @export
#'
#' @importFrom stats pnorm poisson qnorm
#' @importFrom utils tail
mutze_test <- function(data, method = c("nb", "poisson"), conf_level = 0.95) {
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
    if (is.null(fit)) {
      fit <- stats::glm(events ~ treatment + offset(log(tte)), data = df, family = poisson())
      dispersion <- Inf
      method_label <- "Poisson Wald (fallback)"
    } else {
      dispersion <- fit$theta
      method_label <- "Negative Binomial Wald"
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
  pval <- 2 * stats::pnorm(-abs(z))
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

  list(
    method = method_label,
    estimate = est,
    se = se,
    z = z,
    p_value = pval,
    rate_ratio = rr,
    conf_int = rr_ci,
    conf_level = conf_level,
    dispersion = dispersion,
    model = fit,
    group_summary = group_summary
  )
}
