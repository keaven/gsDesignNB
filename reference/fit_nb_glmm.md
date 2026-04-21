# Fit a negative binomial GLMM for count imputation

Fits a negative binomial generalized linear mixed model (GLMM) to
longitudinal count data using
[`glmmTMB::glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html)
with `family = nbinom2(link = "log")`. The fitted model and the
estimated dispersion parameter \\k\\ (where Var(Y) = \\\mu + k\mu^2\\)
are returned for use in subsequent imputation steps. This mirrors PROC
GLIMMIX with `dist=negbin link=log` in SAS.

## Usage

``` r
fit_nb_glmm(data, formula, replicate_col = NULL)
```

## Arguments

- data:

  Data frame of **observed** (non-missing) records.

- formula:

  A two-sided formula specifying fixed and random effects, e.g.
  `count ~ base + trt + visit + (1 | id)`. The left-hand side should be
  the outcome variable; rows with a missing outcome should be excluded
  before calling this function (see
  [`impute_nb()`](https://keaven.github.io/gsDesignNB/reference/impute_nb.md)).

- replicate_col:

  Character or `NULL`. Column in `data` identifying bootstrap
  replicates. The model is fitted separately within each unique value.
  If `NULL` (default), all rows form a single fit.

## Value

A named list—one element per unique replicate—where each element
contains:

- `model`:

  The fitted `glmmTMB` object.

- `k`:

  Estimated NB dispersion \\k\\ (= `glmmTMB::sigma(model)`).

When `replicate_col = NULL`, the list has a single element named `"1"`.

## Details

The dispersion parameter follows the parameterisation used throughout
**gsDesignNB**: \\k\\ such that \\\text{Var}(Y) = \mu + k\mu^2\\. This
matches `glmmTMB`'s `nbinom2` family, where
[`glmmTMB::sigma()`](https://rdrr.io/r/stats/sigma.html) returns \\k\\
directly.

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires glmmTMB
obs_data <- long_data[!is.na(long_data$count), ]
fits <- fit_nb_glmm(
  data    = obs_data,
  formula = count ~ baseline + trt + visit + (1 | id)
)
fits[["1"]]$k  # estimated dispersion
} # }
```
