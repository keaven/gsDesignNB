# Summarize group sequential simulation results

Provides a summary of the operating characteristics of the group
sequential design based on simulation results.

## Usage

``` r
summarize_gs_sim(x)
```

## Arguments

- x:

  A data frame returned by
  [`check_gs_bound()`](https://keaven.github.io/gsDesignNB/reference/check_gs_bound.md)
  (or
  [`sim_gs_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_gs_nbinom.md)
  if bounds are manually checked). Must contain columns `cross_upper`,
  `cross_lower`.

## Value

A list containing:

- n_sim:

  Number of simulations

- power:

  Overall power (probability of crossing upper bound)

- futility:

  Overall futility rate (probability of crossing lower bound and not
  upper)

- analysis_summary:

  Data frame with per-analysis statistics (events, crossings)

## Examples

``` r
design <- gsDesign::gsDesign(k = 2, n.fix = 80, test.type = 2, timing = c(0.5, 1))
sim_df <- data.frame(
  sim = c(1, 1, 2, 2),
  analysis = c(1, 2, 1, 2),
  z_stat = c(2.4, NA, -0.5, 1.9),
  blinded_info = c(40, 80, 40, 80),
  unblinded_info = c(40, 80, 40, 80),
  n_enrolled = c(30, 60, 30, 60),
  events_total = c(12, 25, 10, 22)
)
bounds_checked <- check_gs_bound(sim_df, design)
summarize_gs_sim(bounds_checked)
#> $n_sim
#> [1] 2
#> 
#> $power
#> [1] 0
#> 
#> $futility
#> [1] 0
#> 
#> $analysis_summary
#>   analysis n_enrolled events info_blinded info_unblinded n_cross_upper
#> 1        1         30   11.0           40             40             0
#> 2        2         60   23.5           80             80             0
#>   n_cross_lower prob_cross_upper prob_cross_lower cum_prob_upper
#> 1             0                0                0              0
#> 2             0                0                0              0
#> 
```
