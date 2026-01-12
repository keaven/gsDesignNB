# Compute statistical information at analysis time

Computes the statistical information for the log rate ratio at a given
analysis time, accounting for staggered enrollment and varying exposure
times.

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

  The calendar time of the analysis.

- accrual_rate:

  The enrollment rate (subjects per time unit).

- accrual_duration:

  The duration of the enrollment period.

- lambda1:

  Event rate for group 1 (control).

- lambda2:

  Event rate for group 2 (treatment).

- dispersion:

  The negative binomial dispersion parameter.

- ratio:

  Allocation ratio (n2/n1). Default is 1.

- dropout_rate:

  Dropout rate (hazard rate). Default is 0.

- event_gap:

  Gap duration after each event during which no new events are counted.
  Default is 0.

- max_followup:

  Maximum follow-up time per subject. Exposure time is truncated at this
  value. Default is `Inf` (no truncation).

## Value

The statistical information (inverse of variance) at the analysis time.

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
#> [1] 51.9802
```
