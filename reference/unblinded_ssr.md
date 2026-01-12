# Unblinded sample size re-estimation for recurrent events

Estimates the event rates and dispersion from unblinded interim data and
calculates the required sample size to maintain power, assuming the
planned treatment effect holds (or using the observed control rate).

## Usage

``` r
unblinded_ssr(
  data,
  ratio = 1,
  lambda1_planning,
  lambda2_planning,
  rr0 = 1,
  power = 0.8,
  alpha = 0.025,
  accrual_rate,
  accrual_duration,
  trial_duration,
  dropout_rate = 0,
  max_followup = NULL,
  event_gap = NULL
)
```

## Arguments

- data:

  A data frame containing the unblinded interim data. Must include
  columns `events` (number of events), `tte` (total exposure/follow-up
  time), and `treatment` (treatment group identifier, e.g., 1 for
  control, 2 for experimental). This is typically the output of
  [`cut_data_by_date()`](https://keaven.github.io/gsDesignNB/reference/cut_data_by_date.md).

- ratio:

  Planned allocation ratio (experimental / control). Default is 1.

- lambda1_planning:

  Planned event rate for the control group used in original calculation.

- lambda2_planning:

  Planned event rate for the experimental group used in original
  calculation.

- rr0:

  Rate ratio under the null hypothesis (lambda2/lambda1). Default is 1.

- power:

  Target power (1 - beta). Default is 0.8.

- alpha:

  One-sided significance level. Default is 0.025.

- accrual_rate:

  Vector of accrual rates (patients per unit time).

- accrual_duration:

  Vector of durations for each accrual rate. Must be same length as
  `accrual_rate`.

- trial_duration:

  Total planned duration of the trial.

- dropout_rate:

  Dropout rate (hazard rate). Default is 0.

- max_followup:

  Maximum follow-up time for any patient. Default is NULL (infinite).

- event_gap:

  Gap duration after each event during which no new events are counted.
  Default is NULL (no gap).

## Value

A list containing:

- n_total_unblinded:

  Re-estimated total sample size using unblinded estimates.

- dispersion_unblinded:

  Estimated dispersion parameter (k) from unblinded data.

- lambda1_unblinded:

  Estimated control event rate from unblinded data.

- lambda2_unblinded:

  Estimated experimental event rate from unblinded data.

- info_fraction:

  Estimated information fraction at interim (unblinded information /
  target information).

- unblinded_info:

  Estimated statistical information from the unblinded interim data.

- target_info:

  Target statistical information required for the planned power.
