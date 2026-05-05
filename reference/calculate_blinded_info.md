# Calculate blinded statistical information

Estimates the blinded dispersion and event rate from aggregated interim
data and calculates the observed statistical information for the log
rate ratio, assuming the planned allocation ratio and treatment effect.

## Usage

``` r
calculate_blinded_info(
  data,
  ratio = 1,
  lambda1_planning,
  lambda2_planning,
  event_gap = NULL
)
```

## Arguments

- data:

  A data frame containing the blinded interim data. Must include columns
  `events` (number of events) and `tte` (total exposure/follow-up time).

- ratio:

  Planned allocation ratio (experimental / control). Default is 1.

- lambda1_planning:

  Planned event rate for the control group.

- lambda2_planning:

  Planned event rate for the experimental group.

- event_gap:

  Optional. Gap duration (numeric) to adjust planning rates if provided.
  If provided, planning rates are adjusted as lambda / (1 + lambda \*
  gap).

## Value

A list containing:

- blinded_info:

  Estimated statistical information.

- dispersion_blinded:

  Estimated dispersion parameter (k).

- lambda_blinded:

  Estimated overall event rate.

- lambda1_adjusted:

  Re-estimated control rate.

- lambda2_adjusted:

  Re-estimated experimental rate.
