# Wald test for treatment effect using negative binomial model (Mütze et al.)

Fits a negative binomial (or Poisson) log-rate model to the aggregated
subject-level data produced by
[`cut_data_by_date()`](https://keaven.github.io/gsDesignNB/reference/cut_data_by_date.md).
The method matches the Wald test described by Mütze et al. (2019) for
comparing treatment arms with recurrent event outcomes.

## Usage

``` r
mutze_test(data, method = c("nb", "poisson"), conf_level = 0.95)
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

## Value

A list containing the fitted model summary with elements:

- `estimate`: log rate ratio (experimental vs control).

- `se`: standard error for the log rate ratio.

- `z`: Wald statistic.

- `p_value`: two-sided p-value.

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
#> $method
#> [1] "Negative binomial Wald"
#> 
#> $estimate
#> [1] -0.4016237
#> 
#> $se
#> [1] 0.4928251
#> 
#> $z
#> [1] -0.8149418
#> 
#> $p_value
#> [1] 0.4151056
#> 
#> $rate_ratio
#> [1] 0.6692325
#> 
#> $conf_int
#> [1] 0.254732 1.758209
#> 
#> $conf_level
#> [1] 0.95
#> 
#> $dispersion
#> [1] 5580.973
#> 
#> $model
#> 
#> Call:  MASS::glm.nb(formula = events ~ treatment + offset(log(tte)), 
#>     data = df, init.theta = 5580.973267, link = log)
#> 
#> Coefficients:
#>           (Intercept)  treatmentExperimental  
#>               -0.7111                -0.4016  
#> 
#> Degrees of Freedom: 39 Total (i.e. Null);  38 Residual
#> Null Deviance:       35.59 
#> Residual Deviance: 34.92     AIC: 72.15
#> 
#> $group_summary
#>      treatment subjects events exposure
#> 1      Control       20     10 20.36196
#> 2 Experimental       20      7 21.29804
#> 
```
