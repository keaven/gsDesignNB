#' Method of Moments Estimation for Negative Binomial Parameters
#'
#' Estimates the event rate(s) and common dispersion parameter (k) for
#' negative binomial count data using the method of moments.
#' This is a robust alternative to Maximum Likelihood Estimation (MLE),
#' especially when MLE fails to converge or produces boundary estimates.
#'
#' @param data A data frame containing the data. Must include columns
#'   `events` (number of events) and `tte` (total exposure/follow-up time).
#' @param group Optional character string specifying the grouping column name
#'   (e.g., "treatment"). If provided, rates are estimated separately for each
#'   group, while a common dispersion parameter is estimated across groups.
#'   If NULL (default), a single rate and dispersion are estimated (blinded case).
#'
#' @return A list containing:
#'   \item{lambda}{Estimated event rate(s). A single numeric value if `group` is NULL,
#'     or a named vector if `group` is provided.}
#'   \item{dispersion}{Estimated common dispersion parameter (k).}
#'
#' @details
#' The method of moments estimator for the dispersion parameter \eqn{k} is derived
#' by equating the theoretical variance to the observed second central moment,
#' accounting for varying exposure times.
#'
#' For a given group with rate \eqn{\lambda}, the expected count for subject \eqn{i}
#' is \eqn{\mu_i = \lambda t_i}. The variance is \eqn{V_i = \mu_i + k \mu_i^2}.
#' The estimator is calculated as:
#' \deqn{\hat{k} = \max\left(0, \frac{\sum (y_i - \hat{\mu}_i)^2 - \sum y_i}{\sum \hat{\mu}_i^2}\right)}
#' where \eqn{y_i} is the number of events, \eqn{t_i} is the exposure time,
#' and \eqn{\hat{\mu}_i = \hat{\lambda} t_i} is the estimated expected count.
#'
#' When multiple groups are present, the numerator and denominator are summed
#' across all groups to estimate a common \eqn{k}.
#'
#' @export
#'
#' @examples
#' # Blinded estimation (single group)
#' df <- data.frame(events = c(1, 2, 0, 3), tte = c(1, 1.2, 0.5, 1.5))
#' estimate_nb_mom(df)
#'
#' # Unblinded estimation (two groups)
#' df_group <- df
#' df_group$group <- c("A", "A", "B", "B")
#' estimate_nb_mom(df_group, group = "group")
estimate_nb_mom <- function(data, group = NULL) {
  # Input validation
  if (!is.data.frame(data)) stop("'data' must be a data frame")
  if (!all(c("events", "tte") %in% names(data))) {
    stop("Data must contain 'events' and 'tte' columns.")
  }

  # Filter out rows with non-positive exposure (cannot contribute to rate estimation)
  obs <- data[data$tte > 0, , drop = FALSE]
  
  if (nrow(obs) == 0) {
    warning("No data with positive exposure time.")
    return(list(lambda = NA_real_, dispersion = NA_real_))
  }

  if (is.null(group)) {
    # Blinded case: single rate
    total_events <- sum(obs$events)
    total_tte <- sum(obs$tte)
    
    if (total_tte == 0) {
      lambda_est <- NA_real_
      k_est <- NA_real_
    } else {
      lambda_est <- total_events / total_tte
      
      # Calculate fitted values and residuals
      mu_hat <- lambda_est * obs$tte
      
      # Numerator: Sum((y - mu)^2) - Sum(y)
      # Denominator: Sum(mu^2)
      num <- sum((obs$events - mu_hat)^2) - total_events
      den <- sum(mu_hat^2)
      
      if (den > 0) {
        k_est <- max(0, num / den)
      } else {
        k_est <- 0 # If estimated means are all 0 (i.e. lambda=0)
      }
    }
    
    return(list(lambda = lambda_est, dispersion = k_est))
    
  } else {
    # Unblinded case: rate per group, common dispersion
    if (!group %in% names(obs)) {
      stop(paste("Column", group, "not found in data."))
    }
    
    obs$grp_factor <- as.factor(obs[[group]])
    
    # Estimate rates per group
    # tapply returns array, convert to vector
    sums_events <- tapply(obs$events, obs$grp_factor, sum)
    sums_tte <- tapply(obs$tte, obs$grp_factor, sum)
    
    lambda_est <- sums_events / sums_tte
    lambda_est[is.na(lambda_est)] <- 0 # Handle groups with 0 tte (though filtered out, grouping might create empty levels)
    # Actually if sums_tte is 0, lambda is undefined.
    
    # Assign group lambda to each observation
    obs$lambda_grp <- lambda_est[obs$grp_factor]
    
    # Expected counts
    obs$mu_hat <- obs$lambda_grp * obs$tte
    
    # Pooled MoM for k
    # Sum over all observations (numerator and denominator additive across groups)
    total_events <- sum(obs$events)
    num <- sum((obs$events - obs$mu_hat)^2) - total_events
    den <- sum(obs$mu_hat^2)
    
    if (den > 0) {
      k_est <- max(0, num / den)
    } else {
      k_est <- 0
    }
    
    return(list(lambda = lambda_est, dispersion = k_est))
  }
}
