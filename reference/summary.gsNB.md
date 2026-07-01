# Summary for gsNB objects

Provides a textual summary of a group sequential design for negative
binomial outcomes, similar to the summary provided by
[`gsDesign::gsDesign()`](https://keaven.github.io/gsDesign/reference/gsDesign.html).
For tabular output, use
[`gsDesign::gsBoundSummary()`](https://keaven.github.io/gsDesign/reference/gsBoundSummary.html)
directly on the gsNB object.

## Usage

``` r
# S3 method for class 'gsNB'
summary(object, ...)
```

## Arguments

- object:

  An object of class `gsNB`.

- ...:

  Additional arguments (currently ignored).

## Value

A character string summarizing the design (invisibly). The summary is
also printed to the console.

## Examples

``` r
nb_ss <- sample_size_nbinom(
  lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.9,
  accrual_rate = 10, accrual_duration = 20, trial_duration = 24
)
gs_design <- gsNBCalendar(nb_ss, k = 3, analysis_times = c(12, 18, 24))
summary(gs_design)
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

# For tabular bounds summary, use gsBoundSummary() directly:
gsBoundSummary(gs_design)
#>   Analysis                   Value Efficacy Futility
#>  IA 1: 32%                       Z   3.0107  -0.2387
#>      N: 14             p (1-sided)   0.0013   0.5943
#>  Month: 12         ~delta at bound  -0.8118   0.0644
#>                P(Cross) if delta=0   0.0013   0.4057
#>            P(Cross) if delta=-0.51   0.1412   0.0148
#>  IA 2: 64%                       Z   2.5465   0.9411
#>      N: 29             p (1-sided)   0.0054   0.1733
#>  Month: 18         ~delta at bound  -0.4834  -0.1786
#>                P(Cross) if delta=0   0.0062   0.8347
#>            P(Cross) if delta=-0.51   0.5815   0.0437
#>      Final                       Z   1.9992   1.9992
#>      N: 44             p (1-sided)   0.0228   0.0228
#>  Month: 24         ~delta at bound  -0.3046  -0.3046
#>                P(Cross) if delta=0   0.0233   0.9767
#>            P(Cross) if delta=-0.51   0.9000   0.1000
```
