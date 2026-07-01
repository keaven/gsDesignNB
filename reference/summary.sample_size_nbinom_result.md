# Summary for sample_size_nbinom_result objects

Provides a textual summary of the sample size calculation for negative
binomial outcomes, similar to the summary for gsNB objects.

## Usage

``` r
# S3 method for class 'sample_size_nbinom_result'
summary(object, ...)
```

## Arguments

- object:

  An object of class `sample_size_nbinom_result`.

- ...:

  Additional arguments (currently ignored).

## Value

A character string summarizing the design (invisibly). The summary is
also printed to the console.

## Examples

``` r
x <- sample_size_nbinom(
  lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
  accrual_rate = 10, accrual_duration = 20, trial_duration = 24
)
class(x)
#> [1] "sample_size_nbinom_result" "list"                     
summary(x)
#> Fixed sample size design for negative binomial outcome, total sample size 38
#> (n1=19, n2=19), 80 percent power, 2.5 percent (1-sided) Type I error. Control
#> rate 0.5000, treatment rate 0.3000, risk ratio 0.6000, dispersion 0.1000.
#> Accrual duration 20.0, trial duration 24.0, average exposure 14.00. Expected
#> events 212.8. Randomization ratio 1:1.
#> 
```
