# Summarize adaptive SSR simulation results

Produces trial-level and analysis-level summaries from the output of
[`sim_ssr_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_ssr_nbinom.md).
The summary includes expected sample size, stopping probabilities,
expected participants with events, and expected events observed.

## Usage

``` r
summarize_ssr_sim(x, by = "strategy")
```

## Arguments

- x:

  A `sim_ssr_nbinom` object or a data frame containing the
  `trial_results` component returned by
  [`sim_ssr_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_ssr_nbinom.md).

- by:

  Character vector of grouping columns for the trial-level summary.
  Defaults to `"strategy"`.

## Value

A list with:

- trial_summary:

  Trial-level grouped summary.

- analysis_summary:

  Analysis-level grouped summary.

## Examples

``` r
set.seed(123)
enroll_rate <- data.frame(rate = 10, duration = 4)
fail_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0.5, 0.35),
  dispersion = 0.3
)
fixed_design <- sample_size_nbinom(
  lambda1 = 0.5,
  lambda2 = 0.35,
  dispersion = 0.3,
  power = 0.8,
  alpha = 0.025,
  accrual_rate = 10,
  accrual_duration = 4,
  trial_duration = 8,
  max_followup = 4
)
gs_design <- gsNBCalendar(
  fixed_design,
  k = 3,
  test.type = 4,
  alpha = 0.025,
  analysis_times = c(3, 5, 8)
)
sim_res <- sim_ssr_nbinom(
  n_sims = 2,
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  max_followup = 4,
  design = gs_design,
  strategies = "No adaptation",
  seed = 321
)
summarize_ssr_sim(sim_res)$trial_summary
#>        strategy n_sims rejection_rate futility_rate mean_n sd_n pct_adapted
#> 1 No adaptation      2              1             0    240    0           0
#>   expected_participants_with_events expected_events_observed pct_reach_ia1
#> 1                             153.5                    349.5           100
#>   mean_if_ia1 mean_ia1_time mean_participants_with_events_ia1
#> 1    0.348474      11.05364                                65
#>   mean_events_observed_ia1 n_fallback_ia1 cross_ia1 pct_futility_ia1
#> 1                    161.5              0         0                0
#>   pct_reach_ia2 mean_if_ia2 mean_ia2_time mean_participants_with_events_ia2
#> 1           100   0.7963431      22.37044                             153.5
#>   mean_events_observed_ia2 n_fallback_ia2 cross_ia2 pct_futility_ia2
#> 1                    349.5              0         1                0
#>   pct_reach_final mean_if_final mean_final_time
#> 1               0           NaN             NaN
#>   mean_participants_with_events_final mean_events_observed_final
#> 1                                 NaN                        NaN
#>   n_fallback_final cross_final pct_futility_final mean_adapt_cut_time
#> 1                0           0                  0            22.37044
#>   mean_adapt_enroll_pct mean_adapt_months_to_close pct_adapt_allowed
#> 1              96.04167                  0.8668991               100
#>   pct_adapt_applied
#> 1                 0
```
