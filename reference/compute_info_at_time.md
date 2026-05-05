# Compute Statistical Information at Analysis Time

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
  event_gap = 0
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

## Value

The statistical information (inverse of variance) at the analysis time.
