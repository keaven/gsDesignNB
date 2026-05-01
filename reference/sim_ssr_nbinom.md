# Simulate adaptive group sequential trials with sample size re-estimation

Simulates recurrent-event group sequential trials with information-based
interim analyses and optional sample size re-estimation (SSR). Interim
timing follows the blinded-information targeting used in the SSR study
vignette:
[`get_cut_date()`](https://keaven.github.io/gsDesignNB/reference/get_cut_date.md)
searches for the earliest cut date where the blinded information reaches
the planned information fraction, with an unblinded fallback if the
blinded fit is unavailable.

## Usage

``` r
sim_ssr_nbinom(
  n_sims,
  enroll_rate,
  fail_rate,
  dropout_rate = NULL,
  max_followup,
  design,
  n_max = NULL,
  strategies = c("No adaptation", "Blinded SSR", "Unblinded SSR"),
  adapt_analysis = NULL,
  min_if_futility = 0,
  max_enrollment_frac_for_adapt = 1,
  min_months_to_close_for_adapt = 0,
  analysis_lag_months = 0,
  event_gap = NULL,
  bound_info = c("unblinded_ml", "blinded_ml", "unblinded_mom", "blinded_mom"),
  first_min_time = 1,
  min_analysis_gap = 0.5,
  ignore_futility = FALSE,
  metadata = NULL,
  test_type = c("wald", "score"),
  seed = TRUE
)
```

## Arguments

- n_sims:

  Number of simulated trials.

- enroll_rate:

  Enrollment-rate data frame passed to
  [`nb_sim()`](https://keaven.github.io/gsDesignNB/reference/nb_sim.md).

- fail_rate:

  Failure-rate data frame passed to
  [`nb_sim()`](https://keaven.github.io/gsDesignNB/reference/nb_sim.md).
  This defines the data-generating truth.

- dropout_rate:

  Optional dropout-rate data frame passed to
  [`nb_sim()`](https://keaven.github.io/gsDesignNB/reference/nb_sim.md).

- max_followup:

  Maximum follow-up per patient in the simulated trial.

- design:

  A planning object of class `gsNB`, typically returned by
  [`gsNBCalendar()`](https://keaven.github.io/gsDesignNB/reference/gsNBCalendar.md)
  (optionally after
  [`toInteger()`](https://keaven.github.io/gsDesignNB/reference/toInteger.md)).

- n_max:

  Maximum total enrollment allowed after SSR. Defaults to 150% of the
  planned final enrollment, rounded up and constrained to be at least
  the planned sample size.

- strategies:

  Character vector of strategies to simulate. Must be chosen from
  `"No adaptation"`, `"Blinded SSR"`, and `"Unblinded SSR"`.

- adapt_analysis:

  Interim analysis index where SSR may be applied. Defaults to the last
  interim analysis (`design$k - 1`).

- min_if_futility:

  Minimum observed information fraction required before allowing a
  futility stop. Default is `0`.

- max_enrollment_frac_for_adapt:

  Maximum fraction of the planned enrollment already accrued at the
  adaptation cut for SSR to be allowed. Default is `1`.

- min_months_to_close_for_adapt:

  Minimum predicted months remaining to planned enrollment close
  required to allow SSR. Default is `0`.

- analysis_lag_months:

  Additional months of enrollment counted after a futility stop to
  approximate operational lag. Default is `0`.

- event_gap:

  Optional event-gap duration. If `NULL`, inherits
  `design$nb_design$inputs$event_gap` when available; otherwise defaults
  to `0`.

- bound_info:

  Information measure used for efficacy/futility bounds. Choices are
  `"unblinded_ml"` (default), `"blinded_ml"`, `"unblinded_mom"`, and
  `"blinded_mom"`.

- first_min_time:

  Minimum calendar time allowed for the first information-based interim
  search. Default is `1`.

- min_analysis_gap:

  Minimum gap between successive information-based interim searches.
  Default is `0.5`.

- ignore_futility:

  Logical. If `TRUE`, lower-bound crossings do not stop the trial.

- metadata:

  Optional named list or one-row data frame of scenario labels to repeat
  across the returned rows.

- test_type:

  Type of test statistic passed to
  [`mutze_test()`](https://keaven.github.io/gsDesignNB/reference/mutze_test.md):
  `"wald"` (default) or `"score"`. See
  [`mutze_test()`](https://keaven.github.io/gsDesignNB/reference/mutze_test.md)
  for details.

- seed:

  Random-seed control passed to
  [`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html).
  Follows the same conventions as
  [`sim_gs_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_gs_nbinom.md).

## Value

An object of class `sim_ssr_nbinom` with components:

- trial_results:

  One row per simulation and strategy, containing trial-level outcomes
  and stage-specific summaries in wide format.

- analysis_results:

  One row per simulation, strategy, and analysis in long format.

- settings:

  A list of key design and simulation settings.

## Details

The function compares one or more strategies:

- "No adaptation":

  The trial keeps the planned sample size.

- "Blinded SSR":

  A blinded nuisance-parameter update is applied at `adapt_analysis`
  using
  [`blinded_ssr()`](https://keaven.github.io/gsDesignNB/reference/blinded_ssr.md).

- "Unblinded SSR":

  An unblinded nuisance-parameter update is applied at `adapt_analysis`
  using
  [`unblinded_ssr()`](https://keaven.github.io/gsDesignNB/reference/unblinded_ssr.md).

The returned `trial_results` include stage-specific columns such as
`z_ia1`, `if_ia1`, `ia1_time`, `participants_with_events_ia1`,
`events_observed_ia1`, and analogous columns for later analyses. For the
actual stopping stage, `participants_with_events_stop` and
`events_observed_stop` summarize the observed burden carried into the
final decision for each simulated trial.

## Examples

``` r
set.seed(123)
enroll_rate <- data.frame(rate = 12, duration = 6)
fail_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0.6, 0.42),
  dispersion = 0.4
)
dropout_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0.05, 0.05),
  duration = c(12, 12)
)
fixed_design <- sample_size_nbinom(
  lambda1 = 0.6,
  lambda2 = 0.42,
  dispersion = 0.4,
  power = 0.8,
  alpha = 0.025,
  accrual_rate = 12,
  accrual_duration = 6,
  trial_duration = 12,
  max_followup = 6
)
gs_design <- gsNBCalendar(
  fixed_design,
  k = 3,
  test.type = 4,
  alpha = 0.025,
  sfu = sfHSD,
  sfupar = -2,
  sfl = sfHSD,
  sflpar = 1,
  analysis_times = c(4, 8, 12)
)
sim_res <- sim_ssr_nbinom(
  n_sims = 2,
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  dropout_rate = dropout_rate,
  max_followup = 6,
  design = gs_design,
  seed = 123
)
names(sim_res)
#> [1] "trial_results"    "analysis_results" "settings"        
head(sim_res$trial_results)
#>                sim      strategy reject futility reject_stage futility_stage
#> No adaptation    1 No adaptation   TRUE    FALSE          IA2    No futility
#> Blinded SSR      1   Blinded SSR   TRUE    FALSE          IA2    No futility
#> Unblinded SSR    1 Unblinded SSR   TRUE    FALSE          IA2    No futility
#> No adaptation1   2 No adaptation   TRUE    FALSE          IA1    No futility
#> Blinded SSR1     2   Blinded SSR   TRUE    FALSE          IA1    No futility
#> Unblinded SSR1   2 Unblinded SSR   TRUE    FALSE          IA1    No futility
#>                stop_stage n_adapted adapted adapt_stage adapt_cut_time
#> No adaptation         IA2       229   FALSE         IA2       36.24644
#> Blinded SSR           IA2       229   FALSE         IA2       36.24644
#> Unblinded SSR         IA2       229   FALSE         IA2       36.24644
#> No adaptation1        IA1       229   FALSE         IA2             NA
#> Blinded SSR1          IA1       229   FALSE         IA2             NA
#> Unblinded SSR1        IA1       229   FALSE         IA2             NA
#>                adapt_enroll_frac adapt_months_to_close_pred adapt_allowed
#> No adaptation                  1                          0          TRUE
#> Blinded SSR                    1                          0          TRUE
#> Unblinded SSR                  1                          0          TRUE
#> No adaptation1                NA                         NA         FALSE
#> Blinded SSR1                  NA                         NA         FALSE
#> Unblinded SSR1                NA                         NA         FALSE
#>                adapt_applied participants_with_events_stop events_observed_stop
#> No adaptation          FALSE                           178                  645
#> Blinded SSR            FALSE                           178                  645
#> Unblinded SSR          FALSE                           178                  645
#> No adaptation1         FALSE                            64                  174
#> Blinded SSR1           FALSE                            64                  174
#> Unblinded SSR1         FALSE                            64                  174
#>                ia2_adapt_cut_time ia2_enroll_frac ia2_months_to_close_pred
#> No adaptation            36.24644               1                        0
#> Blinded SSR              36.24644               1                        0
#> Unblinded SSR            36.24644               1                        0
#> No adaptation1                 NA              NA                       NA
#> Blinded SSR1                   NA              NA                       NA
#> Unblinded SSR1                 NA              NA                       NA
#>                ia2_adapt_allowed ia2_adapt_applied    z_ia1    if_ia1 info_ia1
#> No adaptation               TRUE             FALSE 1.519561 0.3288325 25.45254
#> Blinded SSR                 TRUE             FALSE 1.519561 0.3288325 25.45254
#> Unblinded SSR               TRUE             FALSE 1.519561 0.3288325 25.45254
#> No adaptation1             FALSE             FALSE 4.065476 0.3467588 26.84008
#> Blinded SSR1               FALSE             FALSE 4.065476 0.3467588 26.84008
#> Unblinded SSR1             FALSE             FALSE 4.065476 0.3467588 26.84008
#>                ia1_time n_enrolled_ia1 participants_with_events_ia1
#> No adaptation  9.521823            124                           78
#> Blinded SSR    9.521823            124                           78
#> Unblinded SSR  9.521823            124                           78
#> No adaptation1 7.874376             87                           64
#> Blinded SSR1   7.874376             87                           64
#> Unblinded SSR1 7.874376             87                           64
#>                events_observed_ia1 ia1_info_method ia1_fallback    z_ia2
#> No adaptation                  245         blinded        FALSE 3.063249
#> Blinded SSR                    245         blinded        FALSE 3.063249
#> Unblinded SSR                  245         blinded        FALSE 3.063249
#> No adaptation1                 174         blinded        FALSE       NA
#> Blinded SSR1                   174         blinded        FALSE       NA
#> Unblinded SSR1                 174         blinded        FALSE       NA
#>                   if_ia2 info_ia2 ia2_time n_enrolled_ia2
#> No adaptation  0.8536041 66.07131 36.24644            229
#> Blinded SSR    0.8536041 66.07131 36.24644            229
#> Unblinded SSR  0.8536041 66.07131 36.24644            229
#> No adaptation1        NA       NA       NA             NA
#> Blinded SSR1          NA       NA       NA             NA
#> Unblinded SSR1        NA       NA       NA             NA
#>                participants_with_events_ia2 events_observed_ia2 ia2_info_method
#> No adaptation                           178                 645         blinded
#> Blinded SSR                             178                 645         blinded
#> Unblinded SSR                           178                 645         blinded
#> No adaptation1                           NA                  NA         blinded
#> Blinded SSR1                             NA                  NA         blinded
#> Unblinded SSR1                           NA                  NA         blinded
#>                ia2_fallback z_final if_final info_final final_time
#> No adaptation         FALSE      NA       NA         NA         NA
#> Blinded SSR           FALSE      NA       NA         NA         NA
#> Unblinded SSR         FALSE      NA       NA         NA         NA
#> No adaptation1        FALSE      NA       NA         NA         NA
#> Blinded SSR1          FALSE      NA       NA         NA         NA
#> Unblinded SSR1        FALSE      NA       NA         NA         NA
#>                n_enrolled_final participants_with_events_final
#> No adaptation                NA                             NA
#> Blinded SSR                  NA                             NA
#> Unblinded SSR                NA                             NA
#> No adaptation1               NA                             NA
#> Blinded SSR1                 NA                             NA
#> Unblinded SSR1               NA                             NA
#>                events_observed_final final_info_method final_fallback
#> No adaptation                     NA              <NA>          FALSE
#> Blinded SSR                       NA              <NA>          FALSE
#> Unblinded SSR                     NA              <NA>          FALSE
#> No adaptation1                    NA              <NA>          FALSE
#> Blinded SSR1                      NA              <NA>          FALSE
#> Unblinded SSR1                    NA              <NA>          FALSE
```
