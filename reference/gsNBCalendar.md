# Group sequential design for negative binomial outcomes

Creates a group sequential design for negative binomial outcomes based
on sample size calculations from
[`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md).

## Usage

``` r
gsNBCalendar(
  x,
  k = 3,
  test.type = 4,
  alpha = 0.025,
  beta = 0.1,
  astar = 0,
  delta = 0,
  sfu = gsDesign::sfHSD,
  sfupar = -4,
  sfl = gsDesign::sfHSD,
  sflpar = -2,
  sfharm = gsDesign::sfHSD,
  sfharmparam = -2,
  testUpper = TRUE,
  testLower = TRUE,
  testHarm = TRUE,
  tol = 1e-06,
  r = 18,
  usTime = NULL,
  lsTime = NULL,
  analysis_times = NULL
)
```

## Arguments

- x:

  An object of class `sample_size_nbinom_result` from
  [`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md).

- k:

  Number of analyses (interim + final). Default is 3.

- test.type:

  Test type as in
  [`gsDesign::gsDesign()`](https://keaven.github.io/gsDesign/reference/gsDesign.html):

  1

  :   One-sided

  2

  :   Two-sided symmetric

  3

  :   Two-sided, asymmetric, binding futility bound, beta-spending

  4

  :   Two-sided, asymmetric, non-binding futility bound, beta-spending

  5

  :   Two-sided, asymmetric, binding futility bound, lower spending

  6

  :   Two-sided, asymmetric, non-binding futility bound, lower spending

  7

  :   Two-sided, asymmetric, binding futility and binding harm bounds

  8

  :   Two-sided, asymmetric, non-binding futility and non-binding harm
      bounds

  Default is 4.

- alpha:

  Type I error (one-sided). Default is 0.025.

- beta:

  Type II error (1 - power). Default is 0.1.

- astar:

  For test.type 5 or 6, allocated Type I error for the lower bound. For
  test.type 7 or 8, total probability of crossing the harm bound under
  the null. Default is 0 (set to `1 - alpha` internally when needed).

- delta:

  Standardized effect size. Default is 0 (computed from design).

- sfu:

  Spending function for upper bound. Default is
  [`gsDesign::sfHSD`](https://keaven.github.io/gsDesign/reference/sfHSD.html).

- sfupar:

  Parameter for upper spending function. Default is -4.

- sfl:

  Spending function for lower bound. Default is
  [`gsDesign::sfHSD`](https://keaven.github.io/gsDesign/reference/sfHSD.html).

- sflpar:

  Parameter for lower spending function. Default is -2.

- sfharm:

  Spending function for the harm bound (test.type 7 or 8). Default is
  [`gsDesign::sfHSD`](https://keaven.github.io/gsDesign/reference/sfHSD.html).

- sfharmparam:

  Parameter for the harm spending function. Default is -2.

- testUpper:

  Logical scalar or vector of length `k` specifying which analyses
  include an upper (efficacy) bound. `TRUE` (default) means all
  analyses. Where `FALSE`, the upper bound is set to `+20` (effectively
  `Inf`) and displayed as `NA`. Must be `TRUE` at the final analysis.

- testLower:

  Logical scalar or vector of length `k` specifying which analyses
  include a lower (futility) bound. `TRUE` (default) means all analyses.
  Where `FALSE`, the lower bound is set to `-20` (effectively `-Inf`)
  and displayed as `NA`. Ignored for test.type 1.

- testHarm:

  Logical scalar or vector of length `k` specifying which analyses
  include a harm bound (test.type 7 or 8 only). `TRUE` (default) means
  all analyses. Where `FALSE`, the harm bound is set to `-20` and
  displayed as `NA`.

- tol:

  Tolerance for convergence. Default is 1e-06.

- r:

  Integer controlling grid size for numerical integration. Default is
  18.

- usTime:

  Spending time for upper bound (optional).

- lsTime:

  Spending time for lower bound (optional).

- analysis_times:

  Vector of calendar times for each analysis. Must have length k. These
  times are stored in the `T` element and displayed by
  [`gsDesign::gsBoundSummary()`](https://keaven.github.io/gsDesign/reference/gsBoundSummary.html).

## Value

An object of class `gsNB` which inherits from `gsDesign` and
`sample_size_nbinom_result`. While the final sample size would be
planned total enrollment, interim analysis sample sizes are the expected
number enrolled at the times specified in `analysis_times`. Output value
contains all elements from
[`gsDesign::gsDesign()`](https://keaven.github.io/gsDesign/reference/gsDesign.html)
plus:

- nb_design:

  The original `sample_size_nbinom_result` object

- n1:

  A vector with sample size per analysis for group 1

- n2:

  A vector with sample size per analysis for group 2

- n_total:

  A vector with total sample size per analysis

- events:

  A vector with expected total events per analysis

- events1:

  A vector with expected events per analysis for group 1

- events2:

  A vector with expected events per analysis for group 2

- exposure:

  A vector with expected average calendar exposure per analysis

- exposure_at_risk1:

  A vector with expected at-risk exposure per analysis for group 1

- exposure_at_risk2:

  A vector with expected at-risk exposure per analysis for group 2

- variance:

  A vector with variance of log rate ratio per analysis

- T:

  Calendar time at each analysis (if `analysis_times` provided)

- testUpper:

  Logical vector indicating which analyses have an efficacy bound

- testLower:

  Logical vector indicating which analyses have a futility bound

- testHarm:

  Logical vector indicating which analyses have a harm bound (test.type
  7 or 8 only)

Note that `n.I` in the returned object represents the statistical
information at each analysis.

## References

Jennison, C. and Turnbull, B.W. (2000), *Group Sequential Methods with
Applications to Clinical Trials*. Boca Raton: Chapman and Hall.

## Examples

``` r
# First create a sample size calculation
nb_ss <- sample_size_nbinom(
  lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
  accrual_rate = 10, accrual_duration = 20, trial_duration = 24
)

# Then create a group sequential design with analysis times
gs_design <- gsNBCalendar(nb_ss,
  k = 3, test.type = 4,
  analysis_times = c(10, 18, 24)
)

# Selective bound testing: futility only at IA1, efficacy deferred to IA2+
gs_selective <- gsNBCalendar(nb_ss,
  k = 3, test.type = 4,
  analysis_times = c(10, 18, 24),
  testUpper = c(FALSE, TRUE, TRUE),
  testLower = c(TRUE, TRUE, FALSE)
)
gs_selective
#> Asymmetric two-sided group sequential design with
#> 90 % power and 2.5 % Type I Error.
#> Upper bound spending computations assume
#> trial continues if lower bound is crossed.
#> 
#>                 ----Lower bounds----  ----Upper bounds-----
#>   Analysis N    Z   Nominal p Spend+  Z   Nominal p Spend++
#>          1 11 -0.24    0.4057 0.0148 3.01    0.0013  0.0013
#>          2 29  0.94    0.8267 0.0289 2.55    0.0054  0.0049
#>          3 44  2.00    0.9772 0.0563 2.00    0.0228  0.0188
#>      Total                    0.1000                 0.0250 
#> + lower bound beta spending (under H1):
#>  Hwang-Shih-DeCani spending function with gamma = -2.
#> ++ alpha spending:
#>  Hwang-Shih-DeCani spending function with gamma = -4.
#> 
#> Boundary crossing probabilities and expected sample size
#> assume any cross stops the trial
#> 
#> Upper boundary (power or Type I Error)
#>           Analysis
#>    Theta      1      2      3  Total E{N}
#>   0.0000 0.0013 0.0049 0.0171 0.0233 25.4
#>   0.5084 0.1412 0.4403 0.3185 0.9000 32.2
#> 
#> Lower boundary (futility or Type II Error)
#>           Analysis
#>    Theta      1      2      3  Total
#>   0.0000 0.4057 0.4290 0.1420 0.9767
#>   0.5084 0.0148 0.0289 0.0563 0.1000
```
