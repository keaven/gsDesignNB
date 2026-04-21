# Impute missing counts under Missing at Random (MAR)

For observations whose missingness flag matches `mar_values`, generates
`n_imp` imputed counts using the GLMM predicted mean **including
subject-level BLUPs**. Draws use the Gamma–Poisson compound
distribution: \$\$\lambda_i \sim \text{Gamma}(1/k,\\
\hat\mu_i^{\text{BLUP}} k), \quad Y_i^{(m)} \mid \lambda_i \sim
\text{Poisson}(\lambda_i).\$\$

## Usage

``` r
impute_nb_mar(
  data,
  fits,
  outcome_col,
  miss_flag_col,
  mar_values = "MAR",
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

  Character. Column with the count outcome (may have `NA`).

- miss_flag_col:

  Character. Column with the missingness flag.

- mar_values:

  Character vector. Flag values treated as MAR. Default `"MAR"`. In the
  reference-based framework, reference-arm MNAR subjects are also
  imputed as MAR (pass `c("MAR", "MNAR")` if desired).

- n_imp:

  Integer. Number of imputations per replicate. Default `5`.

- replicate_col:

  Character or `NULL`. Replicate identifier column.

## Value

Data frame in long format with all original columns plus `imputation`
(integer 1 to `n_imp`) and `imputed_value` (imputed count; equals the
observed value for non-missing rows).

## Details

This function handles one or more bootstrap replicates when a
`replicate_col` is provided and `fits` contains one model per replicate.
Observed rows (non-missing outcome) are passed through unchanged
(`imputed_value` = observed value).

## Examples

``` r
if (FALSE) { # \dontrun{
fits    <- fit_nb_glmm(obs_data, count ~ base + trt + visit + (1 | id))
imp_mar <- impute_nb_mar(
  data          = long_data,
  fits          = fits,
  outcome_col   = "count",
  miss_flag_col = "miss_flag",
  mar_values    = "MAR",
  n_imp         = 5L
)
} # }
```
