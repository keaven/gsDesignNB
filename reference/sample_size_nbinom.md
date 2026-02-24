# Sample size calculation for negative binomial distribution

Computes the sample size for comparing two treatment groups assuming a
negative binomial distribution for the outcome.

## Usage

``` r
sample_size_nbinom(
  lambda1,
  lambda2,
  dispersion,
  power = NULL,
  alpha = 0.025,
  sided = 1,
  ratio = 1,
  rr0 = 1,
  accrual_rate,
  accrual_duration,
  trial_duration,
  dropout_rate = 0,
  max_followup = NULL,
  event_gap = NULL
)
```

## Arguments

- lambda1:

  Rate in group 1 (control).

- lambda2:

  Rate in group 2 (treatment).

- dispersion:

  Dispersion parameter `k` such that \\\mathrm{Var}(Y) = \mu + k
  \mu^2\\. Note that this is equivalent to `1/size` in R's
  [`stats::rnbinom()`](https://rdrr.io/r/stats/NegBinomial.html)
  parameterization.

- power:

  Power of the test (1 - beta). Default is 0.9.

- alpha:

  Significance level. Default is 0.025.

- sided:

  One-sided or two-sided test. 1 for one-sided, 2 for two-sided. Default
  is 1.

- ratio:

  Allocation ratio n2/n1. Default is 1.

- rr0:

  Rate ratio under the null hypothesis (lambda2/lambda1). Default is 1
  (superiority). For non-inferiority, use a value \> 1 (e.g., 1.1). For
  super-superiority, use a value \< 1 (e.g., 0.8).

- accrual_rate:

  Vector of accrual rates (patients per unit time).

- accrual_duration:

  Vector of durations for each accrual rate. Must be same length as
  `accrual_rate`.

- trial_duration:

  Total planned duration of the trial. If `trial_duration` is less than
  the sum of `accrual_duration`, accrual is truncated at
  `trial_duration`.

- dropout_rate:

  Dropout rate (hazard rate). Default is 0. Can be a vector of length 2.

- max_followup:

  Maximum follow-up time for any patient. Default is NULL (infinite).

- event_gap:

  Gap duration after each event during which no new events are counted.
  Default is NULL (no gap). If provided, the effective event rate is
  reduced.

## Value

An object of class `sample_size_nbinom_result`, which is a list
containing:

- inputs:

  Named list of the original function arguments.

- n1:

  Sample size for group 1

- n2:

  Sample size for group 2

- n_total:

  Total sample size

- alpha:

  Significance level

- sided:

  One-sided or two-sided test

- power:

  Power of the test

- exposure:

  Average exposure time used in calculation (calendar time). Vector of
  length 2.

- exposure_at_risk_n1:

  Average at-risk exposure time for group 1 (accounts for event gap)

- exposure_at_risk_n2:

  Average at-risk exposure time for group 2 (accounts for event gap)

- events_n1:

  Expected number of events in group 1

- events_n2:

  Expected number of events in group 2

- total_events:

  Total expected number of events

- variance:

  Variance of the log rate ratio

- accrual_rate:

  Accrual rate used in calculation

- accrual_duration:

  Accrual duration used in calculation

## References

Zhu, H., & Lakkis, H. (2014). Sample size calculation for comparing two
negative binomial rates. *Statistics in Medicine*, 33(3), 376–387.
[doi:10.1002/sim.5947](https://doi.org/10.1002/sim.5947)

Friede, T., & Schmidli, H. (2010). Blinded sample size reestimation with
negative binomial counts in superiority and non-inferiority trials.
*Methods of Information in Medicine*, 49(06), 618–624.
[doi:10.3414/ME09-02-0060](https://doi.org/10.3414/ME09-02-0060)

Mütze, T., Glimm, E., Schmidli, H., & Friede, T. (2019). Group
sequential designs for negative binomial outcomes. *Statistical Methods
in Medical Research*, 28(8), 2326–2347.
[doi:10.1177/0962280218773115](https://doi.org/10.1177/0962280218773115)

## See also

[`vignette("sample-size-nbinom", package = "gsDesignNB")`](https://keaven.github.io/gsDesignNB/articles/sample-size-nbinom.md)
for a detailed explanation of the methodology.

## Examples

``` r
# Calculate sample size for lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1
# with fixed recruitment of 10/month for 20 months, 24 month trial duration
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

# With piecewise accrual
# 5 patients/month for 3 months, then 10 patients/month for 3 months
# Trial ends at month 12.
x2 <- sample_size_nbinom(
  lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = 0.8,
  accrual_rate = c(5, 10), accrual_duration = c(3, 3),
  trial_duration = 12
)
summary(x2)
#> Fixed sample size design for negative binomial outcome, total sample size 52
#> (n1=26, n2=26), 80 percent power, 2.5 percent (1-sided) Type I error. Control
#> rate 0.5000, treatment rate 0.3000, risk ratio 0.6000, dispersion 0.1000.
#> Accrual duration 6.0, trial duration 12.0, average exposure 8.50. Expected
#> events 176.8. Randomization ratio 1:1.
#> 
```
