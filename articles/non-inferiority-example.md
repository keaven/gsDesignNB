# Non-inferiority and super-superiority designs

``` r
library(gsDesignNB)
```

This vignette demonstrates how to use
[`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md)
for non-inferiority and super-superiority designs by adjusting the `rr0`
parameter (rate ratio under the null hypothesis).

## Non-inferiority design

In a non-inferiority trial, the goal is to demonstrate that a new
treatment is not unacceptably worse than the control. We specify a
non-inferiority margin, which is the maximum acceptable increase in the
event rate for the new treatment.

Suppose we have a control group with an event rate \\\lambda_1 = 0.1\\.
We assume the new treatment has the same event rate \\\lambda_2 =
0.09\\. We want to rule out the possibility that the new treatment
increases the rate by more than 10% (i.e., a rate ratio of 1.1).

- \\\lambda_1 = 0.1\\
- \\\lambda_2 = 0.1\\
- `rr0` = 1.1 (Non-inferiority margin)
- Dispersion \\k = 0.5\\
- Power = 90%
- Alpha = 0.025 (one-sided)

``` r
ni_design <- sample_size_nbinom(
  lambda1 = 0.1,
  lambda2 = 0.09,
  rr0 = 1.1,
  dispersion = 0.5,
  power = 0.9,
  alpha = 0.025,
  trial_duration = 2,
  accrual_duration = 1,
  accrual_rate = 100 # Dummy value to get sample size
)

print(ni_design)
#> Sample size for negative binomial outcome
#> ==========================================
#> 
#> Sample size:     n1 = 3943, n2 = 3943, total = 7886
#> Expected events: 1123.8 (n1: 591.5, n2: 532.3)
#> Power: 90%, Alpha: 0.025 (1-sided)
#> Rates: control = 0.1000, treatment = 0.0900 (RR = 0.9000)
#> Null hypothesis RR: 1.1000
#> Dispersion: 0.5000, Avg exposure (calendar): 1.50
#> Accrual: 1.0, Trial duration: 2.0
```

The output shows the sample size required to reject the null hypothesis
\\H_0: \lambda_2 / \lambda_1 \ge 1.1\\ in favor of \\H_1: \lambda_2 /
\lambda_1 \< 1.1\\, assuming the true ratio is 0.9. This strategy is
useful when the new treatment may have other benefits (e.g., safety,
convenience) that justify a small potential increase in event rate. If
we think there is a slight advantage in efficacy (as in this example),
the required sample size will be much smaller than for a strict
superiority test.

## Super-superiority design

In a super-superiority trial, we want to demonstrate that the new
treatment is not just better, but better by a specific margin.

Suppose we have a control group with an event rate \\\lambda_1 = 0.1\\
and a highly effective new treatment with \\\lambda_2 = 0.02\\. For
instance, if the population is healthy we may want to show a large
reduction in event rates for a new preventive treatment that may have
side effects. This is in line with vaccine trials where large efficacy
margins are often required. We wish to test if the risk reduction is
greater than 50% (i.e., the rate ratio is less than 0.5).

- \\\lambda_1 = 0.1\\
- \\\lambda_2 = 0.02\\
- `rr0` = 0.5 (Super-superiority threshold)
- Dispersion \\k = 0.5\\
- Trial duration = 2
- Max follow-up = 1

``` r
ss_design <- sample_size_nbinom(
  lambda1 = 0.1,
  lambda2 = 0.02,
  rr0 = 0.5,
  dispersion = 0.5,
  power = 0.9,
  alpha = 0.025,
  trial_duration = 2,
  accrual_duration = 1, # Assumed accrual duration
  accrual_rate = 100,
  max_followup = 1
)

print(ss_design)
#> Sample size for negative binomial outcome
#> ==========================================
#> 
#> Sample size:     n1 = 764, n2 = 764, total = 1528
#> Expected events: 91.7 (n1: 76.4, n2: 15.3)
#> Power: 90%, Alpha: 0.025 (1-sided)
#> Rates: control = 0.1000, treatment = 0.0200 (RR = 0.2000)
#> Null hypothesis RR: 0.5000
#> Dispersion: 0.5000, Avg exposure (calendar): 1.00
#> Accrual: 1.0, Trial duration: 2.0
#> Max follow-up: 1.0
```

Here, we are testing \\H_0: \lambda_2 / \lambda_1 \ge 0.5\\ against
\\H_1: \lambda_2 / \lambda_1 \< 0.5\\. Since the assumed effect
(\\\lambda_2/\lambda_1 = 0.2\\) is much stronger than the threshold
(0.5), the power is high or the required sample size is relatively small
compared to a standard superiority test with a smaller effect size.
