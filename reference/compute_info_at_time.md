# Compute statistical information at analysis time

Computes the statistical information \\\mathcal{I}\\ for the log rate
ratio \\\theta = \log(\lambda_2/\lambda_1)\\ at a given calendar
analysis time, accounting for staggered enrollment, dropout, maximum
follow-up, and event gaps.

## Usage

``` r
compute_info_at_time(
  analysis_time,
  accrual_rate,
  accrual_duration,
  lambda1,
  lambda2,
  dispersion,
  ratio = 1,
  dropout_rate = 0,
  event_gap = 0,
  max_followup = Inf
)
```

## Arguments

- analysis_time:

  Calendar time of the analysis.

- accrual_rate:

  Enrollment rate (subjects per time unit).

- accrual_duration:

  Duration of the enrollment period.

- lambda1:

  Event rate \\\lambda_1\\ for group 1 (control).

- lambda2:

  Event rate \\\lambda_2\\ for group 2 (treatment).

- dispersion:

  Dispersion parameter \\k\\ such that \\\mathrm{Var}(Y) = \mu +
  k\mu^2\\. Can be a vector of length 2.

- ratio:

  Allocation ratio \\r = n_2/n_1\\. Default is 1.

- dropout_rate:

  Dropout hazard rate. Default is 0. Can be a vector of length 2 for
  group-specific rates (control, treatment).

- event_gap:

  Gap duration after each event. Default is 0.

- max_followup:

  Maximum follow-up time per subject. Default is `Inf`. Can be a vector
  of length 2.

## Value

The statistical information \\\mathcal{I}\\ (inverse of variance) at the
analysis time.

## Details

This function delegates to
[`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md)
with `power = NULL` and returns \\\mathcal{I} =
1/\mathrm{Var}(\hat\theta)\\ from the resulting variance. This ensures
full consistency with package design calculations, including piecewise
accrual, dropout, max follow-up truncation, event-gap correction, and
follow-up variability inflation (\\Q_g\\).

## Examples

``` r
compute_info_at_time(
  analysis_time = 12,
  accrual_rate = 10,
  accrual_duration = 10,
  lambda1 = 0.5,
  lambda2 = 0.3,
  dispersion = 0.1
)
#> [1] 50.20492
```
