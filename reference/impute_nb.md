# Multiple imputation for longitudinal negative binomial counts

Orchestrates the full multiple imputation (MI) pipeline for longitudinal
recurrent-event count data with negative binomial overdispersion:

## Usage

``` r
impute_nb(
  data,
  formula,
  outcome_col,
  miss_flag_col,
  baseline_col,
  trt_col,
  reference_trt,
  subject_col,
  strata_cols = NULL,
  mar_values = "MAR",
  mnar_value = "MNAR",
  composite_value = "Comp",
  n_imp = 5L,
  n_boot = 1L,
  seed = NULL
)
```

## Arguments

- data:

  Data frame in long format (one row per subject × visit).

- formula:

  Two-sided formula passed to
  [`fit_nb_glmm()`](https://keaven.github.io/gsDesignNB/reference/fit_nb_glmm.md),
  specifying fixed and random effects. The left-hand side should be the
  outcome variable (with `NA` for missing observations). Example:
  `count ~ baseline + trt + visit + (1 | id)`.

- outcome_col:

  Character. Column name of the count outcome.

- miss_flag_col:

  Character. Column name of the missingness mechanism flag. Values in
  this column control which imputation strategy is applied:
  `mar_values`, `mnar_value`, or `composite_value`. Rows with `NA` in
  this column are treated as complete (observed).

- baseline_col:

  Character. Column name of the baseline count used by the composite
  strategy.

- trt_col:

  Character. Column name of the treatment group.

- reference_trt:

  Value in `trt_col` identifying the reference (comparator) arm.

- subject_col:

  Character. Column name of the subject identifier (cluster unit for
  bootstrap resampling).

- strata_cols:

  Character vector of column names used to stratify the bootstrap
  resampling. Default `NULL` (no stratification).

- mar_values:

  Character vector. Values of `miss_flag_col` treated as MAR. Default
  `"MAR"`.

- mnar_value:

  Character. Value of `miss_flag_col` treated as MNAR (triggers
  reference-based imputation for non-reference arms). Default `"MNAR"`.

- composite_value:

  Character. Value of `miss_flag_col` that triggers the composite
  strategy (baseline carry-forward for missing rows). Default `"Comp"`.

- n_imp:

  Integer. Number of imputations per bootstrap replicate. Default `5L`.

- n_boot:

  Integer. Number of bootstrap replicates. Default `1L` (no resampling;
  a single GLMM is fitted to the original data).

- seed:

  Integer or `NULL`. Random seed for reproducibility. Default `NULL`.

## Value

A data frame with all columns from `data` plus:

- `replicate`:

  Bootstrap replicate index (1 to `n_boot`).

- `imputation`:

  Imputation index (1 to `n_imp`).

- `imputed_value`:

  Imputed count. Equals the observed value for non-missing rows;
  contains imputed draws for missing rows.

The total number of rows is `nrow(data) * n_boot * n_imp`.

## Details

1.  **Bootstrap resampling** (optional): cluster-level (subject-level)
    stratified resampling with replacement, creating `n_boot`
    replicates. This propagates estimation uncertainty into the imputed
    values, mirroring the `PROC SURVEYSELECT method=urs cluster=USUBJID`
    step in the SAS macro.

2.  **GLMM fitting**: a negative binomial GLMM is fitted to the observed
    (non-missing) rows of each replicate via
    [`fit_nb_glmm()`](https://keaven.github.io/gsDesignNB/reference/fit_nb_glmm.md).

3.  **Imputation by mechanism**:

    - *MAR* rows: predicted mean with subject BLUPs → Gamma–Poisson
      draw.

    - *MNAR reference-arm* rows: same as MAR (reference arm has no
      "better" treatment to copy from).

    - *MNAR non-reference-arm* rows: reference-based (copy-reference)
      imputation. The counterfactual mean is the fixed-effects-only
      prediction under the reference arm multiplied by the subject's
      random-effect ratio (BLUP prediction / FE prediction on the
      response scale). See
      [`impute_nb_mnar_ref()`](https://keaven.github.io/gsDesignNB/reference/impute_nb_mnar_ref.md).

    - *Composite ICE* rows: missing value set to baseline count. See
      [`impute_nb_composite()`](https://keaven.github.io/gsDesignNB/reference/impute_nb_composite.md).

4.  Returns a long-format data frame with one row per original
    observation × bootstrap replicate × imputation.

### Relationship between bootstrap and MI

Setting `n_boot > 1` combines bootstrap and MI ("boot-MI"), which yields
a valid variance estimator without requiring Rubin's rules. Setting
`n_boot = 1` produces conventional MI; apply Rubin's rules to the
`n_imp` imputed datasets when pooling.

### Formula and GLMM specification

The formula is passed directly to
[`glmmTMB::glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html). A
typical formula mirrors the PROC GLIMMIX model:

    outcome ~ baseline + strat1 + strat2 + trt + visit + param + (1 | id)

The original SAS model also included an unstructured residual covariance
across visits within `id:param`:

    + (0 + visit | id:param)

Complex random-effect structures may cause convergence issues; start
with a random intercept only and add complexity as needed.

### Composite strategy

The composite strategy applies only to **missing** post-ICE rows
(`is.na(outcome_col)` must be `TRUE`). Observed rows with
`miss_flag_col == composite_value` are left unchanged.

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires glmmTMB
result <- impute_nb(
  data          = long_data,
  formula       = count ~ baseline + trt + visit + (1 | id),
  outcome_col   = "count",
  miss_flag_col = "miss_flag",
  baseline_col  = "baseline",
  trt_col       = "trt",
  reference_trt = 0L,
  subject_col   = "id",
  strata_cols   = c("trt", "strat1"),
  mar_values    = "MAR",
  mnar_value    = "MNAR",
  composite_value = "Comp",
  n_imp         = 5L,
  n_boot        = 10L,
  seed          = 42L
)
head(result[!is.na(result$miss_flag), ])
} # }
```
