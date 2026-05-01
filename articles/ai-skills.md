# Using the gsDesignNB AI skill

``` r

library(gsDesignNB)
```

## Purpose

`gsDesignNB` includes two AI-facing documentation files:

- [`SKILL.md`](https://keaven.github.io/gsDesignNB/SKILL.md), a concise
  workflow guide for assistants or humans using the package for negative
  binomial recurrent-event trial design.
- [`llms.txt`](https://keaven.github.io/gsDesignNB/llms.txt), a
  generated pkgdown index that helps language models find the package
  documentation, articles, and reference pages.

This vignette demonstrates how the skill is intended to be used. The
skill does not replace the package documentation, the manuscript, or
statistical review. Instead, it keeps recurring workflows on track: use
package-native functions, align time units, carry event-gap assumptions
consistently, and match sample-size calculations to the planned final
test statistic.

## Example task

Suppose the task is:

> Plan a recurrent-event superiority trial with monthly rates, a 28-day
> inter-event gap, staggered enrollment, dropout, and a final score
> test.

The skill points to the following package-native workflow:

1.  Compute fixed-design sample size with
    [`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md).
2.  Use the score test when Type I error calibration is the priority,
    especially for adaptive or group sequential designs.
3.  Carry `event_gap` through planning, simulation, and data cutting.
4.  Use
    [`gsNBCalendar()`](https://keaven.github.io/gsDesignNB/reference/gsNBCalendar.md)
    for calendar-time group sequential monitoring.
5.  Use `mutze_test(test_type = "score")` for the planned final test.

## Time-scale setup

The most common preventable error in this package is mixing time units.
Here all rates and durations use months. The event gap is 28 days
converted to months.

``` r

event_gap_months <- 28 / 30.4375

design_args <- list(
  lambda1 = 0.08,
  lambda2 = 0.056,
  dispersion = 0.6,
  power = 0.80,
  alpha = 0.025,
  sided = 1,
  accrual_rate = 10,
  accrual_duration = 18,
  trial_duration = 30,
  dropout_rate = 0.01,
  max_followup = 12,
  event_gap = event_gap_months
)
```

## Wald versus score sizing

The skill’s current recommendation is to compare Wald and score sizing,
then choose the final sample size using simulation evidence. The two
calculations use different variance references:

- Wald sizing uses the alternative variance for both the Type I and
  power components.
- Score sizing uses a null variance for the Type I component and an
  alternative variance for the power component.

``` r

wald_design <- do.call(
  sample_size_nbinom,
  c(design_args, list(test_type = "wald"))
)

score_design <- do.call(
  sample_size_nbinom,
  c(design_args, list(test_type = "score"))
)

design_comparison <- data.frame(
  test_type = c(wald_design$test_type, score_design$test_type),
  n_total = c(wald_design$n_total, score_design$n_total),
  n1 = c(wald_design$n1, score_design$n1),
  n2 = c(wald_design$n2, score_design$n2),
  total_events = round(c(wald_design$total_events, score_design$total_events), 1),
  variance_alt = round(c(wald_design$variance, score_design$variance), 4),
  variance_null = round(c(wald_design$variance_null, score_design$variance_null), 4)
)

design_comparison
#>   test_type n_total  n1  n2 total_events variance_alt variance_null
#> 1      wald     518 259 259        361.5       0.0162        0.0159
#> 2     score     512 256 256        357.3       0.0164        0.0161
```

In this scenario, score sizing is slightly smaller than Wald sizing.
That is not a general rule, but it illustrates why the sizing rule and
the analysis test should not be conflated. In the package simulation
grid, the traditional Wald/Zhu–Lakkis sample size paired with the score
test preserved Type I error and provided a small practical power margin.
The skill therefore reminds the analyst to compare sizing rules, choose
the final test deliberately, and verify operating characteristics by
simulation for the actual design setting.

## Calendar-time group sequential design

The same fixed-design result can be passed to
[`gsNBCalendar()`](https://keaven.github.io/gsDesignNB/reference/gsNBCalendar.md)
to construct a calendar-time group sequential design. Here the
Wald-sized fixed design is used as a practical baseline sample size,
while the planned analysis and simulation use the score test for Type I
error control.

``` r

gs_design <- gsNBCalendar(
  wald_design,
  k = 3,
  test.type = 4,
  analysis_times = c(18, 24, 30)
)

data.frame(
  analysis = seq_along(gs_design$n.I),
  calendar_month = c(18, 24, 30),
  planned_information = round(gs_design$n.I, 2),
  information_fraction = round(gs_design$timing, 3)
)
#>   analysis calendar_month planned_information information_fraction
#> 1        1             18               46.88                0.709
#> 2        2             24               61.96                0.937
#> 3        3             30               66.13                1.000
```

## Simulate, cut, and test a small data set

For a quick executable demonstration, simulate a small trial, cut it at
12 months, and run the score test. This is intentionally tiny;
production operating-characteristic work should use
[`sim_gs_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_gs_nbinom.md)
or
[`sim_ssr_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_ssr_nbinom.md)
with many replicates and saved seeds.

``` r

set.seed(2026)

enroll_rate <- data.frame(rate = 30 / 6, duration = 6)
fail_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(design_args$lambda1, design_args$lambda2),
  dispersion = c(design_args$dispersion, design_args$dispersion)
)
dropout_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(design_args$dropout_rate, design_args$dropout_rate),
  duration = c(100, 100)
)

sim_data <- nb_sim(
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  dropout_rate = dropout_rate,
  max_followup = design_args$max_followup,
  n = 60,
  event_gap = design_args$event_gap
)

cut_data <- cut_data_by_date(
  sim_data,
  cut_date = 12,
  event_gap = design_args$event_gap
)

head(cut_data)
#>   id    treatment enroll_time       tte tte_total events
#> 1  1      Control  0.07946939 11.000613 11.920531      1
#> 2  2      Control  0.10208159 11.897918 11.897918      0
#> 3  3 Experimental  0.40356440 11.596436 11.596436      0
#> 4  4 Experimental  0.57077248  9.589392 11.429228      2
#> 5  5 Experimental  0.59292008 11.407080 11.407080      0
#> 6  6      Control  1.40774599  8.621072  9.147637      1
```

``` r

score_test <- mutze_test(cut_data, test_type = "score", sided = 1)
score_test
#> Mutze Test Results
#> ==================
#> 
#> Method:     Negative binomial score 
#> Estimate:   -0.3207
#> SE:         0.6422
#> Z:          -0.5003
#> p-value:    0.3084
#> Rate Ratio: 0.7257
#> CI (95%):  [0.2061, 2.5548]
#> Dispersion: 0.9245
#> 
#> Group Summary:
#>     treatment subjects events exposure
#>       Control       21      8 105.1446
#>  Experimental       20      6 108.6704
```

## Production workflow reminder

The small example above is useful for checking assumptions and object
shapes. For design claims, the skill directs the user to the simulation
routines and summary helpers:

``` r

sim_results <- sim_gs_nbinom(
  n_sims = 10000,
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  dropout_rate = dropout_rate,
  max_followup = design_args$max_followup,
  event_gap = design_args$event_gap,
  design = gs_design,
  analysis_times = c(18, 24, 30),
  test_type = "score",
  seed = TRUE
)

bounded <- check_gs_bound(sim_results, gs_design)
summarize_gs_sim(bounded)
```

For sample size re-estimation studies, use
[`sim_ssr_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_ssr_nbinom.md)
and
[`summarize_ssr_sim()`](https://keaven.github.io/gsDesignNB/reference/summarize_ssr_sim.md).
The score final test is especially important in SSR because adaptation
can increase information under nuisance misspecification; the score test
helps preserve Type I error where a Wald analysis may be mildly
anti-conservative. The adapted sample size itself should still be
checked by simulation rather than assumed from the formula alone.

## What this skill is and is not

The skill is a workflow aid. It is useful for:

- choosing package-native functions rather than reimplementing logic;
- preserving time-scale, event-gap, and test-statistic consistency;
- finding the right vignette or reference page quickly;
- reminding users when simulations are needed to support
  recommendations.

It is not a substitute for protocol-level statistical judgment. Clinical
trial designs still require review of assumptions, estimands,
missing-data handling, operating characteristics, and regulatory
context.
