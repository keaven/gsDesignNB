# Blinded sample size re-estimation for recurrent events

Estimates the blinded dispersion and event rate from aggregated interim
data and calculates the required sample size to maintain power, assuming
the planned treatment effect holds. This function supports constant
rates (Friede & Schmidli 2010) and accommodates future extensions for
time-varying rates (Schneider et al. 2013) by using the
exposure-adjusted rate.

## Usage

``` r
blinded_ssr(
  data,
  ratio = 1,
  lambda1_planning,
  lambda2_planning,
  rr0 = 1,
  power = 0.8,
  alpha = 0.025,
  method = "friede",
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

  A data frame containing the blinded interim data. Must include columns
  `events` (number of events) and `tte` (total exposure/follow-up time).
  This is typically the output of
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

- method:

  Method for sample size recalculation. Currently "friede" (Friede &
  Schmidli 2010) is implemented, which uses the blinded nuisance
  parameter estimates.

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

- n_total_blinded:

  Re-estimated total sample size using blinded estimates.

- dispersion_blinded:

  Estimated dispersion parameter (k) from blinded data.

- lambda_blinded:

  Estimated overall event rate from blinded data.

- info_fraction:

  Estimated information fraction at interim (blinded information /
  target information).

- blinded_info:

  Estimated statistical information from the blinded interim data.

- target_info:

  Target statistical information required for the planned power.

## References

Friede, T., & Schmidli, H. (2010). Blinded sample size reestimation with
count data: methods and applications in multiple sclerosis. *Statistics
in Medicine*, 29(10), 1145–1156.
[doi:10.1002/sim.3861](https://doi.org/10.1002/sim.3861)

Schneider, S., Schmidli, H., & Friede, T. (2013). Blinded sample size
re-estimation for recurrent event data with time trends. *Statistics in
Medicine*, 32(30), 5448–5457.
[doi:10.1002/sim.5977](https://doi.org/10.1002/sim.5977)

## Examples

``` r
interim <- data.frame(events = c(1, 2, 1, 3), tte = c(0.8, 1.0, 1.2, 0.9))
blinded_ssr(
  interim,
  ratio = 1,
  lambda1_planning = 0.5,
  lambda2_planning = 0.3,
  power = 0.8,
  alpha = 0.025,
  accrual_rate = 10,
  accrual_duration = 12,
  trial_duration = 18
)
#> $n_total_blinded
#> (Intercept) 
#>           6 
#> 
#> $dispersion_blinded
#> [1] 1.617394e-05
#> 
#> $lambda_blinded
#> (Intercept) 
#>    1.794874 
#> 
#> $lambda1_adjusted
#> (Intercept) 
#>    2.243592 
#> 
#> $lambda2_adjusted
#> (Intercept) 
#>    1.346155 
#> 
#> $info_fraction
#> [1] 0.05454259
#> 
#> $blinded_info
#> [1] 1.640582
#> 
#> $target_info
#> [1] 30.07893
#> 
```
