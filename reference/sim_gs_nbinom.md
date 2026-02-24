# Simulate group sequential clinical trial for negative binomial outcomes

Simulates multiple replicates of a group sequential clinical trial with
negative binomial outcomes, performing interim analyses at specified
calendar times. Supports parallel execution via the future framework for
faster simulation with reproducible random number generation.

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
  cuts = NULL,
  seed = TRUE
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

- seed:

  Random seed for reproducible simulations. Controls the `future.seed`
  argument of
  [`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html):

  - `TRUE` (default): Automatically generates parallel-safe
    L'Ecuyer-CMRG random number streams. Results are reproducible when
    preceded by [`set.seed()`](https://rdrr.io/r/base/Random.html)
    regardless of the number of workers.

  - An integer: Used as the seed for L'Ecuyer-CMRG streams directly
    (equivalent to calling
    [`set.seed()`](https://rdrr.io/r/base/Random.html) with this value
    before the run).

  - `FALSE` or `NULL`: No special RNG handling (not recommended; results
    may not be reproducible in parallel).

  When future.apply is not installed, `seed` is used with
  [`set.seed()`](https://rdrr.io/r/base/Random.html) for sequential
  execution. See **Details** for parallel usage.

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

- n_ctrl:

  Number of subjects in control group

- n_exp:

  Number of subjects in experimental group

- events_total:

  Total events observed

- events_ctrl:

  Events in control group

- events_exp:

  Events in experimental group

- exposure_at_risk_ctrl:

  Exposure at risk in control group (adjusted for event gaps)

- exposure_at_risk_exp:

  Exposure at risk in experimental group (adjusted for event gaps)

- exposure_total_ctrl:

  Total exposure in control group (calendar follow-up)

- exposure_total_exp:

  Total exposure in experimental group (calendar follow-up)

- z_stat:

  Z-statistic from the Wald test (positive favors experimental if rate
  ratio \< 1)

- estimate:

  Estimated log rate ratio from the model

- se:

  Standard error of the estimate

- method_used:

  Method used for inference ("nb" or "poisson")

- dispersion:

  Estimated dispersion parameter from the model

- blinded_info:

  Estimated blinded statistical information (ML)

- unblinded_info:

  Observed unblinded statistical information (ML)

- info_unblinded_ml:

  Observed unblinded statistical information (ML)

- info_blinded_ml:

  Estimated blinded statistical information (ML)

- info_unblinded_mom:

  Observed unblinded statistical information (Method of Moments)

- info_blinded_mom:

  Estimated blinded statistical information (Method of Moments)

## Details

### Parallel execution

This function uses
[`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html)
to distribute simulation replicates across workers. By default,
simulations run sequentially (equivalent to
[`lapply()`](https://rdrr.io/r/base/lapply.html)). To enable parallel
execution, set a future plan before calling this function:

    library(future)
    plan(multisession, workers = 4)   # use 4 parallel workers
    results <- sim_gs_nbinom(...)
    plan(sequential)                  # restore default

### Reproducibility

The default `seed = TRUE` ensures that results are fully reproducible
regardless of the future plan (sequential or parallel) and regardless of
the number of workers. This is achieved via the L'Ecuyer-CMRG algorithm
which generates statistically independent random number streams for each
simulation replicate. To obtain the same results across runs:

    set.seed(42)
    res1 <- sim_gs_nbinom(n_sims = 100, ..., seed = TRUE)

    set.seed(42)
    res2 <- sim_gs_nbinom(n_sims = 100, ..., seed = TRUE)

    identical(res1, res2)  # TRUE, even with different plan()

## Examples

``` r
# Basic sequential usage with reproducible seed
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
  cuts = cuts,
  seed = TRUE
)
head(sim_results)
#>   sim analysis analysis_time n_enrolled n_ctrl n_exp events_total events_ctrl
#> 1   1        1             2         22     12    10           10           4
#> 2   1        2             4         30     15    15           32          17
#> 3   2        1             2         17      9     8            2           1
#> 4   2        2             4         30     15    15           20          11
#>   events_exp exposure_at_risk_ctrl exposure_at_risk_exp exposure_total_ctrl
#> 1          6              7.903375             9.141767            7.903375
#> 2         15             30.785474            35.601758           30.785474
#> 3          1              7.567636             7.918122            7.567636
#> 4          9             31.001526            33.869952           31.001526
#>   exposure_total_exp      z_stat    estimate        se             method_used
#> 1           9.141767  0.40264747  0.25990122 0.6454808 Poisson Wald (fallback)
#> 2          35.601758 -0.68724537 -0.29289609 0.4261885  Negative binomial Wald
#> 3           7.918122 -0.03201309 -0.04527334 1.4142134 Poisson Wald (fallback)
#> 4          33.869952 -0.62972725 -0.29464890 0.4678992  Negative binomial Wald
#>   dispersion blinded_info unblinded_info info_unblinded_ml info_blinded_ml
#> 1        Inf    2.3997486      2.4001220         2.4001220       2.3997486
#> 2   2.689392    5.2833148      5.5054966         5.5054966       5.2833148
#> 3        Inf    0.4799814      0.5000001         0.5000001       0.4799814
#> 4   7.281096    4.3857137      4.5676764         4.5676764       4.3857137
#>   info_unblinded_mom info_blinded_mom
#> 1           2.400000         2.400000
#> 2           5.801325         5.606178
#> 3           0.500000         0.480000
#> 4           4.692655         4.513762

if (FALSE) { # \dontrun{
# Parallel execution (requires future and future.apply)
library(future)
plan(multisession, workers = 4)
set.seed(42)
sim_results <- sim_gs_nbinom(
  n_sims = 1000,
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  dropout_rate = dropout_rate,
  max_followup = 4,
  n_target = 30,
  design = design,
  cuts = cuts,
  seed = TRUE
)
plan(sequential)
} # }
```
