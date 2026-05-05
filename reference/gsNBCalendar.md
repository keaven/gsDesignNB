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
  timing = 1,
  sfu = gsDesign::sfHSD,
  sfupar = -4,
  sfl = gsDesign::sfHSD,
  sflpar = -2,
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

  Default is 4.

- alpha:

  Type I error (one-sided). Default is 0.025.

- beta:

  Type II error (1 - power). Default is 0.1.

- astar:

  Allocated Type I error for lower bound for test.type = 5 or 6. Default
  is 0.

- delta:

  Standardized effect size. Default is 0 (computed from design).

- timing:

  Timing of interim analyses. May be a vector of length k-1 with values
  between 0 and 1 representing information fractions. Default is 1
  (equally spaced).

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

  Optional vector of calendar times for each analysis. If provided, must
  have length k. These times are stored in the `T` element and displayed
  by
  [`gsDesign::gsBoundSummary()`](https://keaven.github.io/gsDesign/reference/gsBoundSummary.html).

## Value

An object of class `gsNB` which inherits from `gsDesign` and
`sample_size_nbinom_result`. Contains all elements from
[`gsDesign::gsDesign()`](https://keaven.github.io/gsDesign/reference/gsDesign.html)
plus:

- nb_design:

  The original `sample_size_nbinom_result` object

- n1:

  Sample size per analysis for group 1

- n2:

  Sample size per analysis for group 2

- T:

  Calendar time at each analysis (if `analysis_times` provided)

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
```
