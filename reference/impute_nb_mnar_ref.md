# Impute missing counts under a reference-based MNAR assumption

Implements the **copy-reference** (CR) strategy for observations whose
missingness flag matches `mnar_value` in non-reference treatment arms.
The imputation mean is the fixed-effects-only prediction under the
reference arm, adjusted upward (or downward) by the subject's estimated
random effect: \$\$\hat\mu_i^{\text{cf}} = \hat\mu_i^{\text{FE, ref}}
\times \frac{\hat\mu_i^{\text{BLUP}}}{\hat\mu_i^{\text{FE}}}.\$\$ This
mirrors the SAS `PROC PLM` approach that re-predicts under the
counterfactual treatment and then multiplies by the BLUP ratio on the
response scale.

## Usage

``` r
impute_nb_mnar_ref(
  data,
  fits,
  outcome_col,
  miss_flag_col,
  mnar_value = "MNAR",
  trt_col,
  reference_trt,
  n_imp = 5L,
  replicate_col = NULL
)
```

## Arguments

- data:

  Data frame including all rows (observed and missing).

- fits:

  Named list of fits as returned by
  [`fit_nb_glmm()`](https://keaven.github.io/gsDesignNB/reference/fit_nb_glmm.md).

- outcome_col:

  Character. Column with the count outcome.

- miss_flag_col:

  Character. Column with the missingness flag.

- mnar_value:

  Character. Flag value identifying MNAR rows. Default `"MNAR"`.

- trt_col:

  Character. Column with the treatment assignment.

- reference_trt:

  Value in `trt_col` that identifies the reference arm.

- n_imp:

  Integer. Number of imputations per replicate. Default `5`.

- replicate_col:

  Character or `NULL`. Replicate identifier column.

## Value

Data frame in long format with all original columns plus `imputation`
and `imputed_value`. Only MNAR non-reference rows have counterfactual
imputations; all other rows pass through unchanged.

## Details

MNAR subjects already in the reference arm should be handled by
[`impute_nb_mar()`](https://keaven.github.io/gsDesignNB/reference/impute_nb_mar.md)
(MAR imputation is appropriate for the reference arm because there is no
better arm to "copy from").

## Examples

``` r
if (FALSE) { # \dontrun{
fits     <- fit_nb_glmm(obs_data, count ~ base + trt + visit + (1 | id))
imp_mnar <- impute_nb_mnar_ref(
  data          = long_data,
  fits          = fits,
  outcome_col   = "count",
  miss_flag_col = "miss_flag",
  mnar_value    = "MNAR",
  trt_col       = "trt",
  reference_trt = 0L,
  n_imp         = 5L
)
} # }
```
