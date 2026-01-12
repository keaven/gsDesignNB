# Print method for sample_size_nbinom_result objects

Prints a concise summary of the sample size calculation results.

## Usage

``` r
# S3 method for class 'sample_size_nbinom_result'
print(x, ...)
```

## Arguments

- x:

  An object of class `sample_size_nbinom_result`.

- ...:

  Additional arguments (currently ignored).

## Value

Invisibly returns the input object.

## Examples

``` r
x <- sample_size_nbinom(
  lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
  accrual_rate = 10, accrual_duration = 20, trial_duration = 24
)
print(x)
#> Sample size for negative binomial outcome
#> ==========================================
#> 
#> Sample size:     n1 = 19, n2 = 19, total = 38
#> Expected events: 212.8 (n1: 133.0, n2: 79.8)
#> Power: 80%, Alpha: 0.025 (1-sided)
#> Rates: control = 0.5000, treatment = 0.3000 (RR = 0.6000)
#> Dispersion: 0.1000, Avg exposure (calendar): 14.00
#> Accrual: 20.0, Trial duration: 24.0
```
