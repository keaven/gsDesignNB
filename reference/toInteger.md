# Convert group sequential design to integer sample sizes

Generic function to round sample sizes in a group sequential design to
integers. This extends the
[`gsDesign::toInteger()`](https://keaven.github.io/gsDesign/reference/toInteger.html)
function from the gsDesign package to work with `gsNB` objects.

## Usage

``` r
toInteger(x, ...)

# S3 method for class 'gsDesign'
toInteger(x, ratio = x$ratio, roundUpFinal = TRUE, ...)

# S3 method for class 'gsNB'
toInteger(x, ratio = x$nb_design$inputs$ratio, roundUpFinal = TRUE, ...)
```

## Arguments

- x:

  An object of class `gsNB` or `gsDesign`.

- ...:

  Additional arguments passed to methods.

- ratio:

  Randomization ratio (n2/n1). If an integer is provided, rounding is
  done to a multiple of `ratio + 1`. If `ratio < 1` and `1/ratio` is an
  integer (e.g., 1:2 allocation, ratio = 0.5), rounding is done to a
  multiple of `1/ratio + 1`. Default uses the ratio from the original
  design.

- roundUpFinal:

  If `TRUE` (default), the final sample size is rounded up to ensure the
  target is met. If `FALSE`, rounding is to the nearest integer.

## Value

An object of the same class as input with integer sample sizes.

## Details

This function rounds sample sizes at each analysis to integers while
maintaining the randomization ratio and ensuring monotonically
increasing sample sizes across analyses. Only the final analysis sample
size is rounded to an integer; interim sample sizes remain as expected
(non-integer) values based on the information fraction.

When `analysis_times` were provided to
[`gsNBCalendar()`](https://keaven.github.io/gsDesignNB/reference/gsNBCalendar.md),
the statistical information (`n.I`) is recomputed at each analysis time
based on the new sample size and expected exposures.

## Methods (by class)

- `toInteger(gsDesign)`: Method for `gsDesign` objects (calls
  [`gsDesign::toInteger()`](https://keaven.github.io/gsDesign/reference/toInteger.html)).

- `toInteger(gsNB)`: Method for `gsNB` objects.

  Rounds sample sizes in a group sequential negative binomial design to
  integers, respecting the randomization ratio.

## Examples

``` r
nb_ss <- sample_size_nbinom(
  lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
  accrual_rate = 10, accrual_duration = 20, trial_duration = 24
)
gs_design <- gsNBCalendar(nb_ss, k = 3)
#> Error in gsNBCalendar(nb_ss, k = 3): analysis_times must be provided
gs_integer <- toInteger(gs_design)
#> Error: object 'gs_design' not found
```
