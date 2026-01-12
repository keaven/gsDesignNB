# Print method for gsNBsummary objects

Print method for gsNBsummary objects

## Usage

``` r
# S3 method for class 'gsNBsummary'
print(x, ...)
```

## Arguments

- x:

  An object of class `gsNBsummary`.

- ...:

  Additional arguments (currently ignored).

## Value

Invisibly returns the input object.

## Examples

``` r
nb_ss <- sample_size_nbinom(
  lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
  accrual_rate = 10, accrual_duration = 20, trial_duration = 24
)
gs_design <- gsNBCalendar(nb_ss, k = 3)
#> Error in gsNBCalendar(nb_ss, k = 3): analysis_times must be provided
s <- summary(gs_design)
#> Error: object 'gs_design' not found
print(s)
#> Error: object 's' not found
```
