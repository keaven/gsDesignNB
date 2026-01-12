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
gs_design <- gsNBCalendar(nb_ss, k = 3)
#> Error in gsNBCalendar(nb_ss, k = 3): analysis_times must be provided
summary(gs_design)
#> Error: object 'gs_design' not found

# For tabular bounds summary, use gsBoundSummary() directly:
gsBoundSummary(gs_design)
#> Error: object 'gs_design' not found
```
