# ── Internal helpers ───────────────────────────────────────────────────────────

#' Gamma–Poisson compound draw (internal)
#'
#' Draws count observations from the Gamma–Poisson mixture that corresponds to
#' the negative binomial distribution used throughout **gsDesignNB**:
#' \deqn{\lambda_i \sim \text{Gamma}(1/k,\; \mu_i k), \quad
#'   Y_i^{(m)} \mid \lambda_i \sim \text{Poisson}(\lambda_i).}
#'
#' @param mu Numeric vector of predicted means (> 0).
#' @param k  Numeric scalar. NB dispersion: Var(Y) = mu + k * mu^2.
#' @param n_imp Integer. Number of independent draws per element of `mu`.
#' @return Integer matrix, `length(mu)` rows × `n_imp` columns.
#' @keywords internal
.impute_nb_draw <- function(mu, k, n_imp = 1L) {
  n   <- length(mu)
  out <- matrix(NA_integer_, nrow = n, ncol = n_imp)
  shape <- 1.0 / k
  for (m in seq_len(n_imp)) {
    lam     <- stats::rgamma(n, shape = shape, scale = mu * k)
    out[, m] <- stats::rpois(n, lambda = lam)
  }
  out
}

#' Bootstrap-resample subjects within strata (internal)
#'
#' Performs stratified cluster (subject-level) bootstrap resampling. When a
#' subject is drawn more than once its rows are duplicated and the
#' `subject_col` values are replaced with unique pseudo-IDs so that the GLMM
#' treats them as separate individuals (mirroring SAS `PROC SURVEYSELECT` with
#' `method=urs cluster=USUBJID`).
#'
#' @param data        Data frame to resample.
#' @param subject_col Character. Column that identifies subjects (clusters).
#' @param strata_cols Character vector. Stratification column names. May be
#'   `character(0)` for no stratification.
#' @param n_boot      Integer. Number of bootstrap replicates.
#' @return List of length `n_boot`. Each element is a copy of `data` with an
#'   added `replicate` integer column and (when `n_boot > 1`) potentially
#'   modified `subject_col` values for duplicated subjects.
#' @keywords internal
.bootstrap_by_cluster <- function(data, subject_col, strata_cols, n_boot) {
  if (n_boot == 1L) {
    data[["replicate"]] <- 1L
    return(list(data))
  }

  id_strat_cols <- c(subject_col, strata_cols)
  subjects      <- unique(data[, id_strat_cols, drop = FALSE])

  if (length(strata_cols) > 0L) {
    strata_key <- do.call(paste, c(subjects[, strata_cols, drop = FALSE], sep = "|"))
  } else {
    strata_key <- rep(".", nrow(subjects))
  }
  subjects[[".__strat__"]] <- strata_key

  result <- vector("list", n_boot)
  for (b in seq_len(n_boot)) {
    sampled_ids <- unlist(
      lapply(
        split(subjects[[subject_col]], subjects[[".__strat__"]]),
        function(ids) sample(ids, length(ids), replace = TRUE)
      ),
      use.names = FALSE
    )

    id_counts <- table(sampled_ids)
    rows_list <- lapply(names(id_counts), function(sid) {
      subj_rows <- data[data[[subject_col]] == sid, , drop = FALSE]
      cnt       <- as.integer(id_counts[[sid]])
      do.call(rbind, lapply(seq_len(cnt), function(i) {
        r2                <- subj_rows
        r2[[subject_col]] <- if (cnt == 1L) sid else paste0(sid, "_b", i)
        r2
      }))
    })

    boot_df              <- do.call(rbind, rows_list)
    boot_df[["replicate"]] <- b
    rownames(boot_df)    <- NULL
    result[[b]]          <- boot_df
  }
  result
}


# ── Exported functions ─────────────────────────────────────────────────────────

#' Fit a negative binomial GLMM for count imputation
#'
#' Fits a negative binomial generalized linear mixed model (GLMM) to
#' longitudinal count data using `glmmTMB::glmmTMB()` with
#' `family = nbinom2(link = "log")`. The fitted model and the estimated
#' dispersion parameter \eqn{k} (where Var(Y) = \eqn{\mu + k\mu^2}) are
#' returned for use in subsequent imputation steps. This mirrors PROC GLIMMIX
#' with `dist=negbin link=log` in SAS.
#'
#' @param data Data frame of **observed** (non-missing) records.
#' @param formula A two-sided formula specifying fixed and random effects, e.g.
#'   `count ~ base + trt + visit + (1 | id)`. The left-hand side should be the
#'   outcome variable; rows with a missing outcome should be excluded before
#'   calling this function (see [impute_nb()]).
#' @param replicate_col Character or `NULL`. Column in `data` identifying
#'   bootstrap replicates. The model is fitted separately within each unique
#'   value. If `NULL` (default), all rows form a single fit.
#'
#' @return A named list—one element per unique replicate—where each element
#'   contains:
#'   \describe{
#'     \item{`model`}{The fitted `glmmTMB` object.}
#'     \item{`k`}{Estimated NB dispersion \eqn{k} (= `glmmTMB::sigma(model)`).}
#'   }
#'   When `replicate_col = NULL`, the list has a single element named `"1"`.
#'
#' @details
#' The dispersion parameter follows the parameterisation used throughout
#' **gsDesignNB**: \eqn{k} such that \eqn{\text{Var}(Y) = \mu + k\mu^2}. This
#' matches `glmmTMB`'s `nbinom2` family, where `glmmTMB::sigma()` returns
#' \eqn{k} directly.
#'
#' @importFrom stats predict
#' @export
#'
#' @examples
#' \dontrun{
#' # Requires glmmTMB
#' obs_data <- long_data[!is.na(long_data$count), ]
#' fits <- fit_nb_glmm(
#'   data    = obs_data,
#'   formula = count ~ baseline + trt + visit + (1 | id)
#' )
#' fits[["1"]]$k  # estimated dispersion
#' }
fit_nb_glmm <- function(data, formula, replicate_col = NULL) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop(
      "Package 'glmmTMB' is required for model-based imputation.\n",
      "Install it with:  install.packages('glmmTMB')"
    )
  }

  use_internal_rep <- is.null(replicate_col)
  if (use_internal_rep) {
    data[[".__rep__"]] <- 1L
    replicate_col      <- ".__rep__"
  }

  reps   <- sort(unique(data[[replicate_col]]))
  result <- setNames(vector("list", length(reps)), as.character(reps))

  for (r in reps) {
    sub <- data[data[[replicate_col]] == r, , drop = FALSE]
    fit <- glmmTMB::glmmTMB(
      formula = formula,
      data    = sub,
      family  = glmmTMB::nbinom2(link = "log"),
      control = glmmTMB::glmmTMBControl(
        optCtrl = list(iter.max = 1000L, eval.max = 2000L)
      )
    )
    result[[as.character(r)]] <- list(model = fit, k = glmmTMB::sigma(fit))
  }
  result
}


#' Impute missing counts under Missing at Random (MAR)
#'
#' For observations whose missingness flag matches `mar_values`, generates
#' `n_imp` imputed counts using the GLMM predicted mean **including
#' subject-level BLUPs**. Draws use the Gamma–Poisson compound distribution:
#' \deqn{\lambda_i \sim \text{Gamma}(1/k,\; \hat\mu_i^{\text{BLUP}} k),
#'   \quad Y_i^{(m)} \mid \lambda_i \sim \text{Poisson}(\lambda_i).}
#'
#' This function handles one or more bootstrap replicates when a
#' `replicate_col` is provided and `fits` contains one model per replicate.
#' Observed rows (non-missing outcome) are passed through unchanged
#' (`imputed_value` = observed value).
#'
#' @param data         Data frame including all rows (observed and missing).
#' @param fits         Named list of fits as returned by [fit_nb_glmm()].
#' @param outcome_col  Character. Column with the count outcome (may have `NA`).
#' @param miss_flag_col Character. Column with the missingness flag.
#' @param mar_values   Character vector. Flag values treated as MAR. Default
#'   `"MAR"`. In the reference-based framework, reference-arm MNAR subjects are
#'   also imputed as MAR (pass `c("MAR", "MNAR")` if desired).
#' @param n_imp        Integer. Number of imputations per replicate. Default `5`.
#' @param replicate_col Character or `NULL`. Replicate identifier column.
#'
#' @return Data frame in long format with all original columns plus
#'   `imputation` (integer 1 to `n_imp`) and `imputed_value` (imputed count;
#'   equals the observed value for non-missing rows).
#'
#' @importFrom stats predict
#' @export
#'
#' @examples
#' \dontrun{
#' fits    <- fit_nb_glmm(obs_data, count ~ base + trt + visit + (1 | id))
#' imp_mar <- impute_nb_mar(
#'   data          = long_data,
#'   fits          = fits,
#'   outcome_col   = "count",
#'   miss_flag_col = "miss_flag",
#'   mar_values    = "MAR",
#'   n_imp         = 5L
#' )
#' }
impute_nb_mar <- function(data,
                          fits,
                          outcome_col,
                          miss_flag_col,
                          mar_values    = "MAR",
                          n_imp         = 5L,
                          replicate_col = NULL) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop(
      "Package 'glmmTMB' is required.\n",
      "Install it with:  install.packages('glmmTMB')"
    )
  }

  use_internal_rep <- is.null(replicate_col)
  if (use_internal_rep) {
    data[[".__rep__"]] <- 1L
    replicate_col      <- ".__rep__"
  }

  reps     <- sort(unique(data[[replicate_col]]))
  out_list <- vector("list", length(reps))

  for (ri in seq_along(reps)) {
    r      <- reps[ri]
    r_char <- as.character(r)
    sub    <- data[data[[replicate_col]] == r, , drop = FALSE]

    fit_obj <- fits[[r_char]]
    if (is.null(fit_obj)) stop(sprintf("No fit found for replicate '%s'.", r_char))

    model <- fit_obj$model
    k     <- fit_obj$k

    is_mar     <- !is.na(sub[[miss_flag_col]]) & sub[[miss_flag_col]] %in% mar_values
    is_missing <- is.na(sub[[outcome_col]])
    to_impute  <- is_mar & is_missing

    mu_blup <- rep(NA_real_, nrow(sub))
    if (any(to_impute)) {
      mu_blup[to_impute] <- tryCatch(
        stats::predict(model,
          newdata          = sub[to_impute, , drop = FALSE],
          type             = "response",
          re.form          = NULL,
          allow.new.levels = TRUE
        ),
        error = function(e) {
          stats::predict(model,
            newdata          = sub[to_impute, , drop = FALSE],
            type             = "response",
            re.form          = NA,
            allow.new.levels = TRUE
          )
        }
      )
    }

    imp_list <- vector("list", n_imp)
    for (m in seq_len(n_imp)) {
      row_m                    <- sub
      row_m[["imputation"]]    <- m
      row_m[["imputed_value"]] <- sub[[outcome_col]]

      if (any(to_impute)) {
        row_m[["imputed_value"]][to_impute] <-
          .impute_nb_draw(mu_blup[to_impute], k = k, n_imp = 1L)[, 1L]
      }
      imp_list[[m]] <- row_m
    }
    out_list[[ri]] <- do.call(rbind, imp_list)
  }

  result           <- do.call(rbind, out_list)
  rownames(result) <- NULL
  result
}


#' Impute missing counts under a reference-based MNAR assumption
#'
#' Implements the **copy-reference** (CR) strategy for observations whose
#' missingness flag matches `mnar_value` in non-reference treatment arms.
#' The imputation mean is the fixed-effects-only prediction under the reference
#' arm, adjusted upward (or downward) by the subject's estimated random effect:
#' \deqn{\hat\mu_i^{\text{cf}} =
#'   \hat\mu_i^{\text{FE, ref}} \times
#'   \frac{\hat\mu_i^{\text{BLUP}}}{\hat\mu_i^{\text{FE}}}.}
#' This mirrors the SAS `PROC PLM` approach that re-predicts under the
#' counterfactual treatment and then multiplies by the BLUP ratio on the
#' response scale.
#'
#' MNAR subjects already in the reference arm should be handled by
#' [impute_nb_mar()] (MAR imputation is appropriate for the reference arm
#' because there is no better arm to "copy from").
#'
#' @param data          Data frame including all rows (observed and missing).
#' @param fits          Named list of fits as returned by [fit_nb_glmm()].
#' @param outcome_col   Character. Column with the count outcome.
#' @param miss_flag_col Character. Column with the missingness flag.
#' @param mnar_value    Character. Flag value identifying MNAR rows. Default
#'   `"MNAR"`.
#' @param trt_col       Character. Column with the treatment assignment.
#' @param reference_trt Value in `trt_col` that identifies the reference arm.
#' @param n_imp         Integer. Number of imputations per replicate. Default
#'   `5`.
#' @param replicate_col Character or `NULL`. Replicate identifier column.
#'
#' @return Data frame in long format with all original columns plus
#'   `imputation` and `imputed_value`. Only MNAR non-reference rows have
#'   counterfactual imputations; all other rows pass through unchanged.
#'
#' @importFrom stats predict
#' @export
#'
#' @examples
#' \dontrun{
#' fits     <- fit_nb_glmm(obs_data, count ~ base + trt + visit + (1 | id))
#' imp_mnar <- impute_nb_mnar_ref(
#'   data          = long_data,
#'   fits          = fits,
#'   outcome_col   = "count",
#'   miss_flag_col = "miss_flag",
#'   mnar_value    = "MNAR",
#'   trt_col       = "trt",
#'   reference_trt = 0L,
#'   n_imp         = 5L
#' )
#' }
impute_nb_mnar_ref <- function(data,
                               fits,
                               outcome_col,
                               miss_flag_col,
                               mnar_value    = "MNAR",
                               trt_col,
                               reference_trt,
                               n_imp         = 5L,
                               replicate_col = NULL) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop(
      "Package 'glmmTMB' is required.\n",
      "Install it with:  install.packages('glmmTMB')"
    )
  }

  use_internal_rep <- is.null(replicate_col)
  if (use_internal_rep) {
    data[[".__rep__"]] <- 1L
    replicate_col      <- ".__rep__"
  }

  reps     <- sort(unique(data[[replicate_col]]))
  out_list <- vector("list", length(reps))

  for (ri in seq_along(reps)) {
    r      <- reps[ri]
    r_char <- as.character(r)
    sub    <- data[data[[replicate_col]] == r, , drop = FALSE]

    fit_obj <- fits[[r_char]]
    if (is.null(fit_obj)) stop(sprintf("No fit found for replicate '%s'.", r_char))

    model <- fit_obj$model
    k     <- fit_obj$k

    is_mnar    <- !is.na(sub[[miss_flag_col]]) & sub[[miss_flag_col]] == mnar_value
    is_nonref  <- sub[[trt_col]] != reference_trt
    is_missing <- is.na(sub[[outcome_col]])
    to_impute  <- is_mnar & is_nonref & is_missing

    mu_cf <- rep(NA_real_, nrow(sub))
    if (any(to_impute)) {
      mnar_sub <- sub[to_impute, , drop = FALSE]

      mu_blup <- tryCatch(
        stats::predict(model,
          newdata          = mnar_sub,
          type             = "response",
          re.form          = NULL,
          allow.new.levels = TRUE
        ),
        error = function(e)
          stats::predict(model,
            newdata          = mnar_sub,
            type             = "response",
            re.form          = NA,
            allow.new.levels = TRUE
          )
      )
      mu_fe <- stats::predict(model,
        newdata          = mnar_sub,
        type             = "response",
        re.form          = NA,
        allow.new.levels = TRUE
      )

      re_mult                         <- mu_blup / mu_fe
      re_mult[!is.finite(re_mult) | re_mult <= 0] <- 1.0

      cf_sub              <- mnar_sub
      cf_sub[[trt_col]]   <- reference_trt
      cf_mu_fe <- stats::predict(model,
        newdata          = cf_sub,
        type             = "response",
        re.form          = NA,
        allow.new.levels = TRUE
      )

      mu_cf[to_impute] <- cf_mu_fe * re_mult
    }

    imp_list <- vector("list", n_imp)
    for (m in seq_len(n_imp)) {
      row_m                    <- sub
      row_m[["imputation"]]    <- m
      row_m[["imputed_value"]] <- sub[[outcome_col]]

      if (any(to_impute)) {
        row_m[["imputed_value"]][to_impute] <-
          .impute_nb_draw(mu_cf[to_impute], k = k, n_imp = 1L)[, 1L]
      }
      imp_list[[m]] <- row_m
    }
    out_list[[ri]] <- do.call(rbind, imp_list)
  }

  result           <- do.call(rbind, out_list)
  rownames(result) <- NULL
  result
}


#' Apply the composite ICE strategy: replace post-ICE outcomes with baseline
#'
#' For subjects whose missingness flag matches `composite_value`, all missing
#' post-ICE count observations are set to the subject's baseline count. This
#' implements the composite estimand strategy for intercurrent events such as
#' death or treatment discontinuation due to disease worsening, where the
#' event itself is incorporated into the outcome (e.g., baseline carried
#' forward as a "worst case" placeholder).
#'
#' The function is intentionally simple and requires no model. It can be
#' applied to a dataset already containing `imputed_value` from a prior MAR
#' or MNAR imputation step, or directly to the original data.
#'
#' @param data              Data frame.
#' @param outcome_col       Character. Column with the original count outcome;
#'   used to identify which rows are missing (`NA`).
#' @param imputed_value_col Character. Column to update. If absent, it is
#'   created as a copy of `outcome_col` before the composite fill is applied.
#'   Default `"imputed_value"`.
#' @param miss_flag_col     Character. Column with the missingness flag.
#' @param composite_value   Character. Flag value triggering the composite
#'   strategy. Default `"Comp"`.
#' @param baseline_col      Character. Column with the baseline count used as
#'   the fill value.
#'
#' @return Data frame with `imputed_value_col` updated for composite rows.
#'
#' @export
#'
#' @examples
#' df <- data.frame(
#'   count         = c(3L, NA,   NA,    5L),
#'   imputed_value = c(3L,  7L,  NA,    5L),
#'   miss_flag     = c(NA, "MAR", "Comp", NA),
#'   baseline      = c(4L,  4L,   4L,    6L)
#' )
#' impute_nb_composite(
#'   df,
#'   outcome_col     = "count",
#'   miss_flag_col   = "miss_flag",
#'   composite_value = "Comp",
#'   baseline_col    = "baseline"
#' )
impute_nb_composite <- function(data,
                                outcome_col,
                                imputed_value_col = "imputed_value",
                                miss_flag_col,
                                composite_value   = "Comp",
                                baseline_col) {
  if (!imputed_value_col %in% names(data)) {
    data[[imputed_value_col]] <- data[[outcome_col]]
  }

  is_comp   <- !is.na(data[[miss_flag_col]]) & data[[miss_flag_col]] == composite_value
  is_missing <- is.na(data[[outcome_col]])
  to_fill   <- is_comp & is_missing

  data[[imputed_value_col]][to_fill] <- data[[baseline_col]][to_fill]
  data
}


#' Multiple imputation for longitudinal negative binomial counts
#'
#' Orchestrates the full multiple imputation (MI) pipeline for longitudinal
#' recurrent-event count data with negative binomial overdispersion:
#'
#' 1. **Bootstrap resampling** (optional): cluster-level (subject-level)
#'    stratified resampling with replacement, creating `n_boot` replicates.
#'    This propagates estimation uncertainty into the imputed values,
#'    mirroring the `PROC SURVEYSELECT method=urs cluster=USUBJID` step in the
#'    SAS macro.
#' 2. **GLMM fitting**: a negative binomial GLMM is fitted to the observed
#'    (non-missing) rows of each replicate via [fit_nb_glmm()].
#' 3. **Imputation by mechanism**:
#'    - *MAR* rows: predicted mean with subject BLUPs → Gamma–Poisson draw.
#'    - *MNAR reference-arm* rows: same as MAR (reference arm has no
#'      "better" treatment to copy from).
#'    - *MNAR non-reference-arm* rows: reference-based (copy-reference)
#'      imputation. The counterfactual mean is the fixed-effects-only
#'      prediction under the reference arm multiplied by the subject's
#'      random-effect ratio (BLUP prediction / FE prediction on the response
#'      scale). See [impute_nb_mnar_ref()].
#'    - *Composite ICE* rows: missing value set to baseline count.
#'      See [impute_nb_composite()].
#' 4. Returns a long-format data frame with one row per original observation ×
#'    bootstrap replicate × imputation.
#'
#' @param data            Data frame in long format (one row per
#'   subject × visit).
#' @param formula         Two-sided formula passed to [fit_nb_glmm()],
#'   specifying fixed and random effects. The left-hand side should be the
#'   outcome variable (with `NA` for missing observations). Example:
#'   `count ~ baseline + trt + visit + (1 | id)`.
#' @param outcome_col     Character. Column name of the count outcome.
#' @param miss_flag_col   Character. Column name of the missingness mechanism
#'   flag. Values in this column control which imputation strategy is applied:
#'   `mar_values`, `mnar_value`, or `composite_value`. Rows with `NA` in this
#'   column are treated as complete (observed).
#' @param baseline_col    Character. Column name of the baseline count used by
#'   the composite strategy.
#' @param trt_col         Character. Column name of the treatment group.
#' @param reference_trt   Value in `trt_col` identifying the reference
#'   (comparator) arm.
#' @param subject_col     Character. Column name of the subject identifier
#'   (cluster unit for bootstrap resampling).
#' @param strata_cols     Character vector of column names used to stratify the
#'   bootstrap resampling. Default `NULL` (no stratification).
#' @param mar_values      Character vector. Values of `miss_flag_col` treated as
#'   MAR. Default `"MAR"`.
#' @param mnar_value      Character. Value of `miss_flag_col` treated as MNAR
#'   (triggers reference-based imputation for non-reference arms). Default
#'   `"MNAR"`.
#' @param composite_value Character. Value of `miss_flag_col` that triggers the
#'   composite strategy (baseline carry-forward for missing rows). Default
#'   `"Comp"`.
#' @param n_imp           Integer. Number of imputations per bootstrap
#'   replicate. Default `5L`.
#' @param n_boot          Integer. Number of bootstrap replicates. Default
#'   `1L` (no resampling; a single GLMM is fitted to the original data).
#' @param seed            Integer or `NULL`. Random seed for reproducibility.
#'   Default `NULL`.
#'
#' @return A data frame with all columns from `data` plus:
#'   \describe{
#'     \item{`replicate`}{Bootstrap replicate index (1 to `n_boot`).}
#'     \item{`imputation`}{Imputation index (1 to `n_imp`).}
#'     \item{`imputed_value`}{Imputed count. Equals the observed value for
#'       non-missing rows; contains imputed draws for missing rows.}
#'   }
#'   The total number of rows is
#'   `nrow(data) * n_boot * n_imp`.
#'
#' @details
#' ## Relationship between bootstrap and MI
#' Setting `n_boot > 1` combines bootstrap and MI ("boot-MI"), which yields
#' a valid variance estimator without requiring Rubin's rules. Setting
#' `n_boot = 1` produces conventional MI; apply Rubin's rules to the `n_imp`
#' imputed datasets when pooling.
#'
#' ## Formula and GLMM specification
#' The formula is passed directly to `glmmTMB::glmmTMB()`. A typical formula
#' mirrors the PROC GLIMMIX model:
#' ```
#' outcome ~ baseline + strat1 + strat2 + trt + visit + param + (1 | id)
#' ```
#' The original SAS model also included an unstructured residual covariance
#' across visits within `id:param`:
#' ```
#' + (0 + visit | id:param)
#' ```
#' Complex random-effect structures may cause convergence issues; start with
#' a random intercept only and add complexity as needed.
#'
#' ## Composite strategy
#' The composite strategy applies only to **missing** post-ICE rows
#' (`is.na(outcome_col)` must be `TRUE`). Observed rows with `miss_flag_col ==
#' composite_value` are left unchanged.
#'
#' @importFrom stats predict rgamma rpois
#' @export
#'
#' @examples
#' \dontrun{
#' # Requires glmmTMB
#' result <- impute_nb(
#'   data          = long_data,
#'   formula       = count ~ baseline + trt + visit + (1 | id),
#'   outcome_col   = "count",
#'   miss_flag_col = "miss_flag",
#'   baseline_col  = "baseline",
#'   trt_col       = "trt",
#'   reference_trt = 0L,
#'   subject_col   = "id",
#'   strata_cols   = c("trt", "strat1"),
#'   mar_values    = "MAR",
#'   mnar_value    = "MNAR",
#'   composite_value = "Comp",
#'   n_imp         = 5L,
#'   n_boot        = 10L,
#'   seed          = 42L
#' )
#' head(result[!is.na(result$miss_flag), ])
#' }
impute_nb <- function(data,
                      formula,
                      outcome_col,
                      miss_flag_col,
                      baseline_col,
                      trt_col,
                      reference_trt,
                      subject_col,
                      strata_cols     = NULL,
                      mar_values      = "MAR",
                      mnar_value      = "MNAR",
                      composite_value = "Comp",
                      n_imp           = 5L,
                      n_boot          = 1L,
                      seed            = NULL) {
  # ---- Input validation ----
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop(
      "Package 'glmmTMB' is required for model-based imputation.\n",
      "Install it with:  install.packages('glmmTMB')"
    )
  }
  required_cols <- c(outcome_col, miss_flag_col, baseline_col,
                     trt_col, subject_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop("Columns not found in data: ", paste(missing_cols, collapse = ", "))
  }
  strata_cols <- if (is.null(strata_cols)) character(0L) else strata_cols
  if (length(strata_cols) > 0L) {
    bad_strata <- setdiff(strata_cols, names(data))
    if (length(bad_strata) > 0L) {
      stop("Strata columns not found in data: ", paste(bad_strata, collapse = ", "))
    }
  }

  n_imp  <- as.integer(n_imp)
  n_boot <- as.integer(n_boot)
  if (n_imp < 1L)  stop("'n_imp' must be >= 1.")
  if (n_boot < 1L) stop("'n_boot' must be >= 1.")

  if (!is.null(seed)) set.seed(seed)

  # ---- Bootstrap resampling ----
  replicates  <- .bootstrap_by_cluster(data, subject_col, strata_cols, n_boot)
  all_results <- vector("list", n_boot)

  for (b in seq_len(n_boot)) {
    rep_data <- replicates[[b]]

    # ---- Fit GLMM on observed rows ----
    obs_mask <- !is.na(rep_data[[outcome_col]])
    fits_b   <- fit_nb_glmm(
      data         = rep_data[obs_mask, , drop = FALSE],
      formula      = formula,
      replicate_col = NULL        # single fit for this replicate
    )
    model <- fits_b[["1"]]$model
    k     <- fits_b[["1"]]$k

    # ---- Classify each missing row ----
    miss    <- is.na(rep_data[[outcome_col]])
    flag    <- rep_data[[miss_flag_col]]
    is_mar  <- miss & !is.na(flag) & flag %in% mar_values
    is_mnar <- miss & !is.na(flag) & flag == mnar_value
    is_ref  <- rep_data[[trt_col]] == reference_trt
    is_comp <- miss & !is.na(flag) & flag == composite_value

    mar_or_mnar_ref <- is_mar | (is_mnar & is_ref)   # imputed as MAR
    mnar_nonref     <- is_mnar & !is_ref              # reference-based MNAR

    # ---- Pre-compute means (once per replicate, not per imputation) ----
    mu_mar  <- rep(NA_real_, nrow(rep_data))
    mu_mnar <- rep(NA_real_, nrow(rep_data))

    if (any(mar_or_mnar_ref)) {
      sub_mar <- rep_data[mar_or_mnar_ref, , drop = FALSE]
      mu_mar[mar_or_mnar_ref] <- tryCatch(
        stats::predict(model,
          newdata          = sub_mar,
          type             = "response",
          re.form          = NULL,
          allow.new.levels = TRUE
        ),
        error = function(e)
          stats::predict(model,
            newdata          = sub_mar,
            type             = "response",
            re.form          = NA,
            allow.new.levels = TRUE
          )
      )
    }

    if (any(mnar_nonref)) {
      sub_mnar <- rep_data[mnar_nonref, , drop = FALSE]

      mu_blup <- tryCatch(
        stats::predict(model,
          newdata          = sub_mnar,
          type             = "response",
          re.form          = NULL,
          allow.new.levels = TRUE
        ),
        error = function(e)
          stats::predict(model,
            newdata          = sub_mnar,
            type             = "response",
            re.form          = NA,
            allow.new.levels = TRUE
          )
      )
      mu_fe <- stats::predict(model,
        newdata          = sub_mnar,
        type             = "response",
        re.form          = NA,
        allow.new.levels = TRUE
      )

      re_mult <- mu_blup / mu_fe
      re_mult[!is.finite(re_mult) | re_mult <= 0] <- 1.0

      cf_sub            <- sub_mnar
      cf_sub[[trt_col]] <- reference_trt
      cf_mu_fe          <- stats::predict(model,
        newdata          = cf_sub,
        type             = "response",
        re.form          = NA,
        allow.new.levels = TRUE
      )

      mu_mnar[mnar_nonref] <- cf_mu_fe * re_mult
    }

    # ---- Draw n_imp imputed datasets ----
    imp_list <- vector("list", n_imp)
    for (m in seq_len(n_imp)) {
      row_m                    <- rep_data
      row_m[["imputation"]]    <- m
      row_m[["imputed_value"]] <- rep_data[[outcome_col]]

      if (any(mar_or_mnar_ref)) {
        row_m[["imputed_value"]][mar_or_mnar_ref] <-
          .impute_nb_draw(mu_mar[mar_or_mnar_ref], k = k, n_imp = 1L)[, 1L]
      }
      if (any(mnar_nonref)) {
        row_m[["imputed_value"]][mnar_nonref] <-
          .impute_nb_draw(mu_mnar[mnar_nonref], k = k, n_imp = 1L)[, 1L]
      }
      if (any(is_comp)) {
        row_m[["imputed_value"]][is_comp] <- row_m[[baseline_col]][is_comp]
      }

      imp_list[[m]] <- row_m
    }

    all_results[[b]] <- do.call(rbind, imp_list)
  }

  result           <- do.call(rbind, all_results)
  rownames(result) <- NULL
  result
}
