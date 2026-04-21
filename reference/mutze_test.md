# Wald test for treatment effect using negative binomial model (Mutze et al.)

Fits a negative binomial (or Poisson) log-rate model to the aggregated
subject-level data produced by
[`cut_data_by_date()`](https://keaven.github.io/gsDesignNB/reference/cut_data_by_date.md).
The method matches the Wald test described by Mutze et al. (2019) for
comparing treatment arms with recurrent event outcomes. When the maximum
likelihood negative binomial fit is unreliable, the test automatically
switches to one of two statistically sensible fallbacks: a Poisson Wald
test when the data are essentially Poisson, or a method-of-moments (MoM)
variance estimate plugged into the same negative binomial information
formula when the data are extremely overdispersed or the ML fit fails to
converge. The latter avoids the anti-conservative behaviour of a blind
Poisson fallback under genuine overdispersion.

## Usage

``` r
mutze_test(
  data,
  method = c("nb", "poisson"),
  conf_level = 0.95,
  sided = 1,
  poisson_threshold = 50,
  mom_threshold = 20
)

# S3 method for class 'mutze_test'
print(x, ...)
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

  Upper threshold (in units of `fit$theta`, the
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) shape
  parameter \\\theta\_{\text{NB}} = 1/k\\) above which the data are
  treated as essentially Poisson and the function falls back to a
  Poisson Wald test. Default is 50, corresponding to \\\hat{k} \<
  0.02\\, by which point NB and Poisson Wald standard errors are
  numerically indistinguishable at typical trial sample sizes.

- mom_threshold:

  Lower threshold on `fit$theta` below which the NB ML fit is considered
  unreliable (extreme overdispersion). When triggered, or when
  `glm.nb()` fails to converge, the function falls back to a
  method-of-moments NB Wald test: rates and dispersion are re-estimated
  via
  [`estimate_nb_mom()`](https://keaven.github.io/gsDesignNB/reference/estimate_nb_mom.md)
  and the Wald variance is computed from the Fisher information formula
  \\\mathcal{I} = 1/(1/W_1 + 1/W_2)\\ with \\W_g = \sum_i
  \mu\_{g,i}/(1 + \hat{k}\mu\_{g,i})\\. Default is 20, corresponding to
  \\\hat{k} \> 20\\. This avoids the anti-conservative variance of a
  Poisson fallback when the data are truly overdispersed.

- x:

  An object of class `mutze_test`.

- ...:

  Additional arguments (currently ignored).

## Value

An object of class `mutze_test` containing the fitted model summary with
elements:

- `method`: A string indicating the test method used.

- `estimate`: log rate ratio (experimental vs control).

- `se`: standard error for the log rate ratio.

- `z`: Wald statistic.

- `p_value`: one-sided or two-sided p-value.

- `rate_ratio`: estimated rate ratio and its confidence interval.

- `dispersion`: estimated dispersion. For consistency this is reported
  on the [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html)
  scale (\\\theta\_{\text{NB}} = 1/k\\), regardless of whether ML,
  Poisson fallback, or MoM fallback was used. `Inf` indicates Poisson
  (\\k=0\\).

- `group_summary`: observed subjects/events/exposure per treatment.

- `fallback`: character label describing which fit path was used
  (`"ml"`, `"poisson"`, or `"mom"`).

Invisibly returns the input object.

## Methods (by generic)

- `print(mutze_test)`: Print method for `mutze_test` objects.

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
#> Mutze Test Results
#> ==================
#> 
#> Method:     Poisson Wald (fallback, near-Poisson ML) 
#> Estimate:   -0.1709
#> SE:         0.5175
#> Z:          -0.3302
#> p-value:    0.3706
#> Rate Ratio: 0.8429
#> CI (95%):  [0.3057, 2.3245]
#> Dispersion: Inf
#> 
#> Group Summary:
#>     treatment subjects events exposure
#>  Experimental       20      7 21.78820
#>       Control       20      8 20.98945
```
