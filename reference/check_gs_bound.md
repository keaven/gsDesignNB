# Check group sequential bounds

Updates the group sequential design boundaries based on observed
information and checks if boundaries have been crossed.

## Usage

``` r
check_gs_bound(sim_results, design, info_scale = c("blinded", "unblinded"))
```

## Arguments

- sim_results:

  Data frame of simulation results (from
  [`sim_gs_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_gs_nbinom.md)).

- design:

  The planning `gsNB` object.

- info_scale:

  Character. "blinded" (default) or "unblinded" information to use for
  bounds.

## Value

A data frame with added columns:

- cross_upper:

  Logical, true if upper bound crossed (efficacy)

- cross_lower:

  Logical, true if lower bound crossed (futility)

## Examples

``` r
design <- gsDesign::gsDesign(k = 2, n.fix = 100, test.type = 2, timing = c(0.5, 1))
sim_df <- data.frame(
  sim = c(1, 1, 2, 2),
  analysis = c(1, 2, 1, 2),
  z_stat = c(2.5, NA, -0.2, 2.2),
  blinded_info = c(50, 100, 50, 100),
  unblinded_info = c(50, 100, 50, 100)
)
check_gs_bound(sim_df, design)
#>   sim analysis z_stat blinded_info unblinded_info cross_upper cross_lower
#> 1   1        1    2.5           50             50       FALSE       FALSE
#> 2   1        2     NA          100            100       FALSE       FALSE
#> 3   2        1   -0.2           50             50       FALSE       FALSE
#> 4   2        2    2.2          100            100        TRUE       FALSE
```
