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
gs_design <- gsNBCalendar(nb_ss, k = 3, analysis_times = c(12, 18, 24))
s <- summary(gs_design)
print(s)
#> Asymmetric two-sided with non-binding futility bound group sequential design
#> for negative binomial outcomes, 3 analyses, total sample size 53.5, 90 percent
#> power, 2.5 percent (1-sided) Type I error. Control rate 0.5000, treatment rate
#> 0.3000, risk ratio 0.6000, dispersion 0.1000. Accrual duration 20.0, trial
#> duration 24.0, max follow-up Inf, average exposure 14.00. Randomization ratio
#> 1:1. Upper spending: Hwang-Shih-DeCani (gamma = -4) Lower spending:
#> Hwang-Shih-DeCani (gamma = -2)
#> Asymmetric two-sided with non-binding futility bound group sequential design
#> for negative binomial outcomes, 3 analyses, total sample size 53.5, 90 percent
#> power, 2.5 percent (1-sided) Type I error. Control rate 0.5000, treatment rate
#> 0.3000, risk ratio 0.6000, dispersion 0.1000. Accrual duration 20.0, trial
#> duration 24.0, max follow-up Inf, average exposure 14.00. Randomization ratio
#> 1:1. Upper spending: Hwang-Shih-DeCani (gamma = -4) Lower spending:
#> Hwang-Shih-DeCani (gamma = -2)
#> 
```
