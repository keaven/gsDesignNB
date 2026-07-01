# Wald or score test for treatment effect using negative binomial model

Fits a negative binomial (or Poisson) log-rate model to the aggregated
subject-level data produced by
[`cut_data_by_date()`](https://keaven.github.io/gsDesignNB/reference/cut_data_by_date.md).
With `test_type = "wald"` (default), the method matches the Wald test
described by Mutze et al. (2019). With `test_type = "score"`, the
function fits only the null (no treatment effect) model and computes the
score statistic, which evaluates all quantities under \\H_0\\ and avoids
the finite-sample anti-conservatism of the Wald test.

## Usage

``` r
mutze_test(
  data,
  method = c("nb", "poisson"),
  test_type = c("wald", "score"),
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

- test_type:

  Type of test statistic: `"wald"` (default) or `"score"`. The Wald test
  estimates the treatment effect under the alternative and divides by
  its standard error. The score test fits only the null model and
  evaluates the derivative of the log-likelihood at \\\theta = 0\\,
  avoiding estimation under the alternative. The score test typically
  has better finite-sample Type I error control and is faster because it
  only fits a one-parameter null model.

- conf_level:

  Confidence level for the rate ratio interval. Default 0.95.

- sided:

  Number of sides for the test: 1 (default) or 2.

- poisson_threshold:

  Upper threshold (in units of `fit$theta`, the
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) shape
  parameter \\\theta\_{\text{NB}} = 1/k\\) above which the data are
  treated as essentially Poisson. Default is 50, corresponding to
  \\\hat{k} \< 0.02\\.

- mom_threshold:

  Lower threshold on `fit$theta` below which the NB ML fit is considered
  unreliable (extreme overdispersion). Default is 20, corresponding to
  \\\hat{k} \> 20\\.

- x:

  An object of class `mutze_test`.

- ...:

  Additional arguments (currently ignored).

## Value

An object of class `mutze_test` containing:

- `method`: A string indicating the test method used.

- `estimate`: log rate ratio (experimental vs control). For
  `test_type = "score"`, this is a plug-in estimate.

- `se`: standard error for the log rate ratio.

- `z`: test statistic (Wald or score).

- `p_value`: one-sided or two-sided p-value.

- `rate_ratio`: estimated rate ratio and its confidence interval.

- `dispersion`: estimated dispersion on the \\\theta = 1/k\\ scale.

- `group_summary`: observed subjects/events/exposure per treatment.

- `fallback`: character label (`"ml"`, `"poisson"`, or `"mom"`).

- `test_type`: character label (`"wald"` or `"score"`).

Invisibly returns the input object.

## Details

When the maximum likelihood negative binomial fit is unreliable, the
test automatically switches to one of two statistically sensible
fallbacks: a Poisson test when the data are essentially Poisson, or a
method-of-moments (MoM) variance estimate plugged into the same negative
binomial information formula when the data are extremely overdispersed
or the ML fit fails to converge.

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
mutze_test(cut, test_type = "score")
#> Mutze Test Results
#> ==================
#> 
#> Method:     Poisson score (fallback, near-Poisson ML) 
#> Estimate:   -0.1709
#> SE:         0.5165
#> Z:          -0.3306
#> p-value:    0.3705
#> Rate Ratio: 0.8429
#> CI (95%):  [0.3063, 2.3197]
#> Dispersion: Inf
#> 
#> Group Summary:
#>     treatment subjects events exposure
#>  Experimental       20      7 21.78820
#>       Control       20      8 20.98945
```
