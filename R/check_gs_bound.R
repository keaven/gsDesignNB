#' Check group sequential bounds
#'
#' Updates the group sequential design boundaries based on observed information
#' and checks if boundaries have been crossed.
#'
#' @param sim_results Data frame of simulation results (from [sim_gs_nbinom()]).
#' @param design The planning `gsNB` object.
#' @param info_scale Character. "blinded" (default) or "unblinded" information to use for bounds.
#'
#' @return A data frame with added columns:
#'   \describe{
#'     \item{cross_upper}{Logical, true if upper bound crossed (efficacy)}
#'     \item{cross_lower}{Logical, true if lower bound crossed (futility)}
#'   }
#'
#' @export
#'
#' @examples
#' design <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 2, timing = c(0.5, 1))
#' sim_df <- data.frame(
#'   sim = c(1, 1, 2, 2),
#'   analysis = c(1, 2, 1, 2),
#'   z_stat = c(2.5, NA, -0.2, 2.2),
#'   blinded_info = c(50, 100, 50, 100),
#'   unblinded_info = c(50, 100, 50, 100)
#' )
#' check_gs_bound(sim_df, design)
check_gs_bound <- function(sim_results, design, info_scale = c("blinded", "unblinded")) {
  info_scale <- match.arg(info_scale)

  if (!inherits(design, "gsNB") && !inherits(design, "gsDesign")) {
    stop("design must be a gsDesign or gsNB object")
  }

  dt <- data.table::as.data.table(sim_results)

  # Ensure we have z_stat
  if (!"z_stat" %in% names(dt)) stop("sim_results must contain 'z_stat'")

  # Identify information column
  info_col <- if (info_scale == "blinded") "blinded_info" else "unblinded_info"
  if (!info_col %in% names(dt)) stop(paste("sim_results must contain", info_col))

  # Process each simulation separately
  process_sim <- function(sub_dt) {
    k_max <- max(sub_dt$analysis)
    max_info <- design$n.fix

    cross_upper <- rep(FALSE, k_max)
    cross_lower <- rep(FALSE, k_max)
    stopped <- FALSE

    # Current timing vector (starts with design timing)
    current_timing <- design$timing

    for (k in seq_len(k_max)) {
      if (stopped) {
        # Already stopped
        break
      }

      obs_info <- sub_dt[analysis == k][[info_col]]
      z_val <- sub_dt[analysis == k]$z_stat

      if (is.na(z_val)) next

      # Flip Z if needed?
      # In vignette: z_eff <- -test_result$z.
      # mutze_test returns Z for treatment coefficient.
      # If lambda2 < lambda1 (benefit), coeff is negative, Z is negative.
      # gsDesign assumes positive Z for efficacy.
      # So we should flip Z if the design anticipates negative effect.
      # The user vignette flips it. I will assume we flip it.
      # Ideally, we check delta1.
      if (!is.null(design$delta1) && design$delta1 < 0) {
        z_val <- -z_val
      }

      # Calculate information fraction
      frac <- obs_info / max_info

      # Cap fraction at 1 if overrunning
      if (frac > 1) frac <- 1

      # Update timing
      current_timing[k] <- frac

      # Ensure monotonicity
      if (k > 1 && current_timing[k] <= current_timing[k - 1]) {
        current_timing[k] <- current_timing[k - 1] + 0.001
      }
      if (k < k_max && current_timing[k + 1] <= current_timing[k]) {
        # Push future timings if needed (though they will be overwritten later)
        current_timing[k + 1] <- min(current_timing[k] + 0.001, 0.999)
      }

      # Handle final analysis logic
      if (k == k_max) current_timing[k] <- 1

      # Update design bounds
      # We need to recreate the design with the new timing
      # We extract parameters from the original design
      # Note: usTime is stored in gsNB objects (from my recent fix)
      # If present, we might want to respect it, BUT the vignette logic
      # explicitly overwrites timing based on information fraction.
      # Standard error spending usually adapts to information.
      # If usTime was used, gsDesign usually ignores timing?
      # Wait, if usTime is provided, gsDesign uses it for spending.
      # If we want to use Information Fraction spending (Lan-DeMets), we should NOT pass usTime?
      # The vignette recalculates bounds:
      # gsDesign(..., timing = current_timing, ...)
      # It passes usTime too!
      # If usTime is passed, gsDesign uses usTime for spending calculation.
      # So changing `timing` (information fraction) only affects the drift/B-value scaling?
      # Actually, boundaries depend on Information Fraction for covariance.
      # Spending depends on Spending Function argument `t`.
      # If `usTime` is provided, `t = usTime`. If not, `t = timing`.

      # To follow vignette logic:

      # It passes `usTime` AND `timing`.
      # Use tryCatch to handle edge cases where spending function can't
      # achieve target beta with the observed information fraction.
      temp_gs <- tryCatch(
        gsDesign::gsDesign(
          k = design$k,
          test.type = design$test.type,
          alpha = design$alpha,
          beta = design$beta,
          astar = design$astar,
          delta = design$delta,
          sfu = design$upper$sf,
          sfupar = design$upper$param,
          sfl = design$lower$sf,
          sflpar = design$lower$param,
          tol = design$tol,
          r = design$r,
          n.fix = max_info,
          timing = current_timing,
          usTime = design$usTime, # If NULL, ignored
          lsTime = design$lsTime
        ),
        error = function(e) {
          # If the spending function fails, fall back to original design bounds
          # This can happen with extreme information fractions
          design
        }
      )

      upper_bound <- temp_gs$upper$bound[k]
      lower_bound <- temp_gs$lower$bound[k]

      if (z_val > upper_bound) {
        cross_upper[k] <- TRUE
        stopped <- TRUE
      } else if (z_val < lower_bound) {
        cross_lower[k] <- TRUE
        stopped <- TRUE
      }
    }

    return(list(cross_upper = cross_upper, cross_lower = cross_lower))
  }

  # Apply to all sims
  # This might be slow for many sims if done in R loop.
  # But bounds update is fast.

  results <- dt[, process_sim(.SD), by = sim]

  # Merge back
  dt[, cross_upper := results$cross_upper]
  dt[, cross_lower := results$cross_lower]

  as.data.frame(dt)
}
