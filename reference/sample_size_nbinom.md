# Sample size calculation for negative binomial outcomes

Computes the sample size (or power) for comparing two treatment groups
assuming negative binomial distributed event counts. The test is based
on the Wald statistic for the log rate ratio, corresponding to Method 3
of Zhu & Lakkis (2014). The same underlying variance formula is used by
Friede & Schmidli (2010) and Mutze et al. (2019).

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

  Event rate for group 1 (control), in events per unit time.

- lambda2:

  Event rate for group 2 (treatment), in events per unit time.

- dispersion:

  Dispersion parameter \\k\\ such that \\\mathrm{Var}(Y) = \mu +
  k\mu^2\\. Equivalent to `1/size` in
  [`stats::rnbinom()`](https://rdrr.io/r/stats/NegBinomial.html). Can be
  a scalar (common dispersion) or a vector of length 2 (group-specific:
  control, treatment).

- power:

  Target power (\\1 - \beta\\). If `NULL`, power is computed for the
  given accrual rates (no sample size scaling). Default is 0.9.

- alpha:

  Significance level. Default is 0.025.

- sided:

  Number of sides for the test: 1 (one-sided) or 2 (two-sided). Default
  is 1.

- ratio:

  Allocation ratio \\r = n_2/n_1\\. Default is 1 (equal allocation).

- rr0:

  Rate ratio under the null hypothesis (\\\lambda_2 / \lambda_1\\).
  Default is 1 (superiority). For non-inferiority, use a value \> 1
  (e.g., 1.1). For super-superiority, use a value \< 1 (e.g., 0.8).

- accrual_rate:

  Vector of accrual rates (patients per unit time) for each recruitment
  segment.

- accrual_duration:

  Vector of durations for each accrual segment. Must be the same length
  as `accrual_rate`.

- trial_duration:

  Total planned duration of the trial. If `trial_duration` is less than
  the sum of `accrual_duration`, accrual is truncated at
  `trial_duration`.

- dropout_rate:

  Dropout hazard rate. Default is 0. Can be a vector of length 2 for
  group-specific dropout (control, treatment).

- max_followup:

  Maximum follow-up time for any patient. Default is `NULL` (infinite).
  Can be a vector of length 2 for group-specific caps.

- event_gap:

  Gap duration after each event during which no new events are counted
  (e.g., a recovery period). Default is `NULL` (no gap). When specified,
  the effective rate is reduced to \\\lambda\_{\mathrm{eff}} = \lambda /
  (1 + \lambda \cdot \mathrm{gap})\\.

## Value

An object of class `sample_size_nbinom_result`, which is a list
containing:

- inputs:

  Named list of the original function arguments.

- n1:

  Sample size for group 1 (control).

- n2:

  Sample size for group 2 (treatment).

- n_total:

  Total sample size (\\n_1 + n_2\\).

- alpha:

  Significance level used.

- sided:

  One-sided or two-sided test.

- power:

  Power of the test.

- exposure:

  Average calendar exposure \\\bar{t}\_g\\ (vector of length 2 for
  control and treatment).

- exposure_at_risk_n1:

  Average at-risk exposure for group 1 (adjusted for event gap).

- exposure_at_risk_n2:

  Average at-risk exposure for group 2 (adjusted for event gap).

- events_n1:

  Expected number of events in group 1.

- events_n2:

  Expected number of events in group 2.

- total_events:

  Total expected number of events.

- variance:

  Variance of the log rate ratio \\\mathrm{Var}(\hat\theta)\\.

- accrual_rate:

  Accrual rate(s) used (possibly scaled to achieve target power).

- accrual_duration:

  Accrual duration(s) used.

## Details

### Sample size formula

The sample size for group 1 is: \$\$n_1 = \frac{(z\_{\alpha/s} +
z\_\beta)^2 \tilde{V}}{(\theta - \theta_0)^2}\$\$ where \\\theta =
\log(\lambda_2/\lambda_1)\\, \\\theta_0 = \log(\mathrm{rr}\_0)\\, and:
\$\$\tilde{V} = \left(\frac{1}{\mu_1} + k_1\right) +
\frac{1}{r}\left(\frac{1}{\mu_2} + k_2\right)\$\$ with \\\mu_g =
\lambda_g \bar{t}\_g\\ the expected event count and \\\bar{t}\_g\\ the
average exposure for group \\g\\.

### Average exposure

The average exposure \\\bar{t}\_g\\ accounts for piecewise accrual,
exponential dropout, and maximum follow-up truncation. When dropout rate
\\\delta \> 0\\, the expected exposure for a patient with potential
follow-up \\u\\ is \\m(u) = (1 - e^{-\delta u})/\delta\\. The overall
average is a weighted mean across accrual segments.

### Variance inflation

When follow-up times are variable, the dispersion is inflated by a
factor \\Q_g = \mathrm{E}\[t_g^2\] / (\mathrm{E}\[t_g\])^2 \ge 1\\ (Zhu
& Lakkis, 2014) to account for the non-linear dependence of the NB
variance on exposure.

### Event gap correction (Jensen's inequality)

When `event_gap` \> 0, the naive effective rate \\\lambda / (1 + \lambda
g)\\ overestimates the true population-level effective rate because of
subject-level heterogeneity (frailty). In the Gamma-Poisson mixture,
each subject's rate \\\Lambda_i \sim \mathrm{Gamma}(1/k, k\lambda)\\ is
random. Since \\f(x) = x/(1+xg)\\ is concave, Jensen's inequality gives
\\\mathrm{E}\[f(\Lambda)\] \< f(\mathrm{E}\[\Lambda\])\\.

A second-order Taylor correction is applied: \$\$\lambda\_{\mathrm{eff}}
\approx \frac{\lambda}{1+\lambda g} \left(1 - \frac{k \lambda
g}{(1+\lambda g)^2}\right)\$\$ This uses \\f''(\lambda) = -2g/(1+\lambda
g)^3\\ and \\\mathrm{Var}(\Lambda) = k\lambda^2\\.

## References

Zhu, H., & Lakkis, H. (2014). Sample size calculation for comparing two
negative binomial rates. *Statistics in Medicine*, 33(3), 376–387.
[doi:10.1002/sim.5947](https://doi.org/10.1002/sim.5947)

Friede, T., & Schmidli, H. (2010). Blinded sample size reestimation with
negative binomial counts in superiority and non-inferiority trials.
*Methods of Information in Medicine*, 49(06), 618–624.
[doi:10.3414/ME09-02-0060](https://doi.org/10.3414/ME09-02-0060)

Mutze, T., Glimm, E., Schmidli, H., & Friede, T. (2019). Group
sequential designs for negative binomial outcomes. *Statistical Methods
in Medical Research*, 28(8), 2326–2347.
[doi:10.1177/0962280218773115](https://doi.org/10.1177/0962280218773115)

## See also

[`compute_info_at_time()`](https://keaven.github.io/gsDesignNB/reference/compute_info_at_time.md)
for computing statistical information at a given analysis time;
[`blinded_ssr()`](https://keaven.github.io/gsDesignNB/reference/blinded_ssr.md)
for blinded sample size reestimation;
[`gsNBCalendar()`](https://keaven.github.io/gsDesignNB/reference/gsNBCalendar.md)
for group sequential designs;
[`vignette("sample-size-nbinom", package = "gsDesignNB")`](https://keaven.github.io/gsDesignNB/articles/sample-size-nbinom.md)
for detailed methodology.

## Examples

``` r
# Basic sample size calculation
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

# Compute power for a fixed design (power = NULL)
sample_size_nbinom(
  lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.1, power = NULL,
  accrual_rate = 10, accrual_duration = 20, trial_duration = 24
)
#> Sample size for negative binomial outcome
#> ==========================================
#> 
#> Sample size:     n1 = 100, n2 = 100, total = 200
#> Expected events: 1120.0 (n1: 700.0, n2: 420.0)
#> Power: 100%, Alpha: 0.025 (1-sided)
#> Rates: control = 0.5000, treatment = 0.3000 (RR = 0.6000)
#> Dispersion: 0.1000, Avg exposure (calendar): 14.00
#> Accrual: 20.0, Trial duration: 24.0
```
