# Simulate group sequential clinical trial for negative binomial outcomes

Simulates multiple replicates of a group sequential clinical trial with
negative binomial outcomes, performing interim analyses at specified
calendar times.

## Usage

``` r
sim_gs_nbinom(
  n_sims,
  enroll_rate,
  fail_rate,
  dropout_rate = NULL,
  max_followup,
  event_gap = 0,
  analysis_times = NULL,
  n_target = NULL,
  design = NULL,
  data_cut = cut_data_by_date,
  cuts = NULL
)
```

## Arguments

- n_sims:

  Number of simulations to run.

- enroll_rate:

  Enrollment rates (data frame with `rate` and `duration`).

- fail_rate:

  Failure rates (data frame with `treatment`, `rate`, `dispersion`).

- dropout_rate:

  Dropout rates (data frame with `treatment`, `rate`, `duration`).

- max_followup:

  Maximum follow-up time.

- event_gap:

  Event gap duration.

- analysis_times:

  Vector of calendar times for interim and final analyses. Optional if
  `cuts` is provided.

- n_target:

  Total sample size to enroll (optional, if not defined by
  `enroll_rate`).

- design:

  An object of class `gsNB` or `sample_size_nbinom_result`. Used to
  extract planning parameters (`lambda1`, `lambda2`, `ratio`) for
  blinded information estimation.

- data_cut:

  Function to cut data for analysis. Defaults to
  [`cut_data_by_date()`](https://keaven.github.io/gsDesignNB/reference/cut_data_by_date.md).
  The function must accept `sim_data`, `cut_date`, and `event_gap` as
  arguments.

- cuts:

  A list of cutting criteria for each analysis. Each element of the list
  should be a list of arguments for
  [`get_cut_date()`](https://keaven.github.io/gsDesignNB/reference/get_cut_date.md)
  (e.g., `planned_calendar`, `target_events`, `target_info`). If
  provided, `analysis_times` is ignored (or used as a fallback if
  `planned_calendar` is missing in a cut).

## Value

A data frame containing simulation results for each analysis of each
trial. Columns include:

- sim:

  Simulation ID

- analysis:

  Analysis index

- analysis_time:

  Calendar time of analysis

- n_enrolled:

  Number of subjects enrolled

- events_total:

  Total events observed

- events_ctrl:

  Events in control group

- events_exp:

  Events in experimental group

- exposure_ctrl:

  Total exposure in control group

- exposure_exp:

  Total exposure in experimental group

- z_stat:

  Z-statistic from the Wald test (positive favors experimental if rate
  ratio \< 1)

- blinded_info:

  Estimated blinded statistical information

- unblinded_info:

  Observed unblinded statistical information

## Examples

``` r
set.seed(123)
enroll_rate <- data.frame(rate = 10, duration = 3)
fail_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0.6, 0.4),
  dispersion = 0.2
)
dropout_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0.05, 0.05),
  duration = c(6, 6)
)
design <- sample_size_nbinom(
  lambda1 = 0.6, lambda2 = 0.4, dispersion = 0.2, power = 0.8,
  accrual_rate = enroll_rate$rate, accrual_duration = enroll_rate$duration,
  trial_duration = 6
)
cuts <- list(
  list(planned_calendar = 2),
  list(planned_calendar = 4)
)
sim_results <- sim_gs_nbinom(
  n_sims = 2,
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  dropout_rate = dropout_rate,
  max_followup = 4,
  n_target = 30,
  design = design,
  cuts = cuts
)
head(sim_results)
#>   sim analysis analysis_time n_enrolled n_ctrl n_exp events_total events_ctrl
#> 1   1        1             2         23     11    12           10           8
#> 2   1        2             4         30     15    15           34          25
#> 3   2        1             2         15      8     7            4           3
#> 4   2        2             4         30     15    15           24          18
#>   events_exp exposure_at_risk_ctrl exposure_at_risk_exp exposure_total_ctrl
#> 1          2              13.14534            12.328930            13.14534
#> 2          9              37.48061            38.937872            37.48061
#> 3          1               7.44672             5.387439             7.44672
#> 4          6              31.08699            26.047756            31.08699
#>   exposure_total_exp     z_stat   estimate        se             method_used
#> 1          12.328930 -1.6724385 -1.3221756 0.7905675 Poisson Wald (fallback)
#> 2          38.937872 -2.7264734 -1.0597950 0.3887054 Poisson Wald (fallback)
#> 3           5.387439 -0.6710972 -0.7749088 1.1546894 Poisson Wald (fallback)
#> 4          26.047756 -1.7504986 -0.8915788 0.5093285  Negative binomial Wald
#>   dispersion blinded_info unblinded_info
#> 1        Inf     2.391287      1.6000076
#> 2        Inf     8.157183      6.6184888
#> 3        Inf     0.959958      0.7500145
#> 4   4.406507     3.842925      3.8548194
```
