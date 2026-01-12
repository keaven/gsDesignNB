# Wald test for treatment effect using negative binomial model (Mütze et al.)

Fits a negative binomial (or Poisson) log-rate model to the aggregated
subject-level data produced by
[`cut_data_by_date()`](https://keaven.github.io/gsDesignNB/reference/cut_data_by_date.md).
The method matches the Wald test described by Mütze et al. (2019) for
comparing treatment arms with recurrent event outcomes.

## Usage

``` r
mutze_test(
  data,
  method = c("nb", "poisson"),
  conf_level = 0.95,
  sided = 1,
  poisson_threshold = 1000
)
```

## Arguments

- data:

  A data frame with at least the columns `treatment`, `events`, and
  `tte` (follow-up time). Typically output from
  [`cut_data_by_date()`](https://keaven.github.io/gsDesignNB/reference/cut_data_by_date.md).

- method:

  Type of model to fit: "nb" (default) uses a negative binomial GLM via
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html),
  "poisson" fits a Poisson GLM.

- conf_level:

  Confidence level for the rate ratio interval. Default 0.95.

- sided:

  Number of sides for the test: 1 (default) or 2.

- poisson_threshold:

  When `method = "nb"`, the model falls back to Poisson regression if
  theta (the NB shape parameter) is outside the range
  `[1/poisson_threshold, poisson_threshold]`. Very large theta indicates
  near-Poisson data, while very small theta indicates extreme
  overdispersion with unstable estimates. Default is 1000.

## Value

An object of class `mutze_test` containing the fitted model summary with
elements:

- `method`: A string indicating the test method used.

- `estimate`: log rate ratio (experimental vs control).

- `se`: standard error for the log rate ratio.

- `z`: Wald statistic.

- `p_value`: one-sided or two-sided p-value.

- `rate_ratio`: estimated rate ratio and its confidence interval.

- `dispersion`: estimated dispersion (theta) when `method = "nb"`.

- `group_summary`: observed subjects/events/exposure per treatment.

## Examples

``` r
enroll_rate <- data.frame(rate = 20 / (5 / 12), duration = 5 / 12)
fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3))
dropout_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0.1, 0.05), duration = c(100, 100)
)
sim <- nb_sim(enroll_rate, fail_rate, dropout_rate, max_followup = 2, n = 40)
cut <- cut_data_by_date(sim, cut_date = 1.5)
mutze_test(cut)
#> Mütze Test Results
#> ==================
#> 
#> Method:     Poisson Wald (fallback) 
#> Estimate:   -0.0933
#> SE:         0.4859
#> Z:          -0.1921
#> p-value:    0.4238
#> Rate Ratio: 0.9109
#> CI (95%):  [0.3515, 2.3608]
#> Dispersion: Inf
#> 
#> Group Summary:
#>     treatment subjects events exposure
#>  Experimental       20      8 21.02991
#>       Control       20      9 21.55013
```
