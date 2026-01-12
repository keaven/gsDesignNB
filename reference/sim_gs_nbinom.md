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
#>   sim analysis analysis_time n_enrolled events_total events_ctrl events_exp
#> 1   1        1             2         23           10           8          2
#> 2   1        2             4         30           34          25          9
#> 3   2        1             2         15            4           3          1
#> 4   2        2             4         30           24          18          6
#>   exposure_ctrl exposure_exp     z_stat blinded_info unblinded_info
#> 1      13.14534    12.328930 -1.6724385     2.391287      1.6000076
#> 2      37.48061    38.937872 -2.7264734     8.157183      6.6184888
#> 3       7.44672     5.387439 -0.6710972     0.959958      0.7500145
#> 4      31.08699    26.047756 -1.7504986     3.842925      3.8548194
```
