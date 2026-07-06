# Update group sequential bounds with observed information

Given a planned `gsNB` (or `gsDesign`) object and observed statistical
information at one or more analyses, recompute the group sequential
boundaries and return an updated design object together with a
[`gsDesign::gsBoundSummary()`](https://keaven.github.io/gsDesign/reference/gsBoundSummary.html)-style
table.

## Usage

``` r
update_gsNB(design, observed_info, spending_time = NULL)
```

## Arguments

- design:

  A `gsNB` or `gsDesign` object produced by
  [`gsNBCalendar()`](https://keaven.github.io/gsDesignNB/reference/gsNBCalendar.md)
  (or
  [`gsDesign::gsDesign()`](https://keaven.github.io/gsDesign/reference/gsDesign.html)).

- observed_info:

  Numeric vector of observed statistical information at each analysis
  conducted so far. Its length must be between 1 and `design$k`. If
  shorter than `design$k`, information at future analyses is projected
  from the planned design.

- spending_time:

  Optional numeric vector the same length as `observed_info` giving the
  spending time (between 0 and 1) for each analysis. When `NULL`,
  spending time equals the information fraction. If shorter than
  `design$k`, future spending times are taken from the planned design.

## Value

A list with components:

- design:

  The updated `gsDesign` object with recalculated boundaries.

- bounds:

  A data frame from
  [`gsDesign::gsBoundSummary()`](https://keaven.github.io/gsDesign/reference/gsBoundSummary.html)
  showing Z-boundaries, nominal p-values, approximate treatment effects
  at the boundary, and cumulative crossing probabilities at each
  analysis.

- info:

  A data frame with one row per analysis containing the information
  fraction (`IF`), spending time (`spending_time`), upper and lower
  Z-boundaries, and cumulative upper and lower spending.

## Details

The observed information determines the covariance structure of the test
statistics (via the information fraction `timing`), while
`spending_time` controls how much of the error-spending budget has been
used. When `spending_time` is `NULL` (the default), spending is driven
by the observed information fraction. Supplying an explicit
`spending_time` is useful when the monitoring charter specifies
calendar-driven spending that differs from the observed information
fraction.

## Examples

``` r
library(gsDesign)
#> 
#> Attaching package: ‘gsDesign’
#> The following object is masked from ‘package:gsDesignNB’:
#> 
#>     toInteger
nb_ss <- sample_size_nbinom(
  lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
  accrual_rate = 10, accrual_duration = 20, trial_duration = 24
)
gs <- gsNBCalendar(nb_ss, k = 3, analysis_times = c(12, 18, 24))

# After observing information at the first interim
upd <- update_gsNB(gs, observed_info = gs$n.I[1] * 0.95)
upd$bounds
#>   Analysis               Value Efficacy Futility
#>  IA 1: 32%                   Z   3.0253  -0.2814
#>      N: 15         p (1-sided)   0.0012   0.6108
#>                ~delta at bound   1.5860  -0.1475
#>            P(Cross) if delta=0   0.0012   0.3892
#>            P(Cross) if delta=1   0.1318   0.0143
#>  IA 2: 64%                   Z   2.5838   0.8587
#>      N: 28         p (1-sided)   0.0049   0.1953
#>                ~delta at bound   0.9615   0.3195
#>            P(Cross) if delta=0   0.0057   0.8138
#>            P(Cross) if delta=1   0.5475   0.0411
#>      Final                   Z   1.9967   1.9967
#>      N: 44         p (1-sided)   0.0229   0.0229
#>                ~delta at bound   0.5963   0.5963
#>            P(Cross) if delta=0   0.0234   0.9766
#>            P(Cross) if delta=1   0.9000   0.1000
upd$info
#>   Analysis     IF spending_time upper_bound lower_bound cum_upper_spend
#> 1        1 0.3245        0.3245      3.0253     -0.2814        0.001242
#> 2        2 0.6441        0.6441      2.5838      0.8587        0.005668
#> 3        3 1.0000        1.0000      1.9967      1.9967        0.025000
#>   cum_lower_spend
#> 1        0.014302
#> 2        0.041108
#> 3        0.100000
```
