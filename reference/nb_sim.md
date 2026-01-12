# Simulate recurrent events with fixed follow-up

Simulates recurrent events for a clinical trial with piecewise constant
enrollment, exponential failure rates (Poisson process), and piecewise
exponential dropout.

## Usage

``` r
nb_sim(
  enroll_rate,
  fail_rate,
  dropout_rate = NULL,
  max_followup = NULL,
  n = NULL,
  block = c(rep("Control", 2), rep("Experimental", 2)),
  event_gap = 0
)
```

## Arguments

- enroll_rate:

  A data frame with columns `rate` and `duration` defining the piecewise
  constant enrollment rates.

- fail_rate:

  A data frame with columns `treatment` and `rate` defining the
  exponential failure rate for each treatment group. Optionally, a
  `dispersion` column can be provided to generate data from a negative
  binomial distribution. The dispersion parameter `k` is such that
  \\\mathrm{Var}(Y) = \mu + k \mu^2\\.

- dropout_rate:

  A data frame with columns `treatment`, `rate`, and `duration` defining
  the piecewise constant dropout rates.

- max_followup:

  Numeric. Maximum duration of follow-up for each individual (relative
  to their randomization time).

- n:

  Total sample size. If NULL, it is estimated from `enroll_rate`. If
  provided, enrollment stops when `n` subjects are recruited.

- block:

  Block vector for treatment allocation. Default is
  `c(rep("Control", 2), rep("Experimental", 2))`. If NULL, simple
  randomization is used (treatments are assigned with equal
  probability). If provided, it specifies the block structure, for
  example, `c(rep("A", 2), rep("B", 2))` assigns 2 to group A and 2 to
  group B in each block.

- event_gap:

  Numeric. Gap duration after each event during which no new events are
  counted. Default is 0.

## Value

A data frame (tibble) with columns:

- id:

  Subject identifier

- treatment:

  Treatment group

- enroll_time:

  Time of enrollment relative to trial start

- tte:

  Time to event or censoring relative to randomization

- calendar_time:

  Calendar time of event or censoring (enroll_time + tte)

- event:

  Binary indicator: 1 for event, 0 for censoring

Multiple rows per subject are returned (one for each event, plus one for
the final censoring time).

## Details

The simulation generates data consistent with the negative binomial
models described by Friede and Schmidli (2010) and Mütze et al. (2019).
Specifically, it simulates a Gamma-distributed frailty variable for each
individual (if dispersion \> 0), which acts as a multiplier for that
individual's event rate. Events are then generated according to a
Poisson process with this subject-specific rate.

More explicitly, for a subject with baseline rate \\\lambda\\ and
exposure time \\t\\, the model used here is a Gamma–Poisson mixture:
\$\$\Lambda_i \sim \mathrm{Gamma}(\text{shape}=1/k,\\
\text{scale}=k\lambda), \quad Y_i \mid \Lambda_i \sim
\mathrm{Poisson}(\Lambda_i t).\$\$ Marginally, \\Y_i\\ follows a
negative binomial distribution with \\\mathrm{E}\[Y_i\]=\mu=\lambda t\\
and \\\mathrm{Var}(Y_i)=\mu + k\mu^2\\. This \\k\\ is the package
dispersion parameter (and corresponds to \\1/\theta\\ in
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html)
terminology).

## References

Friede, T., & Schmidli, H. (2010). Blinded sample size reestimation with
count data: methods and applications in multiple sclerosis. *Statistics
in Medicine*, 29(10), 1145–1156.
[doi:10.1002/sim.3861](https://doi.org/10.1002/sim.3861)

Mütze, T., Glimm, E., Schmidli, H., & Friede, T. (2019). Group
sequential designs for negative binomial outcomes. *Statistical Methods
in Medical Research*, 28(8), 2326–2347.
[doi:10.1177/0962280218773115](https://doi.org/10.1177/0962280218773115)

## Examples

``` r
enroll_rate <- data.frame(rate = 20 / (5 / 12), duration = 5 / 12)
fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3))
dropout_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0.1, 0.05), duration = c(100, 100)
)
sim <- nb_sim(enroll_rate, fail_rate, dropout_rate, max_followup = 2, n = 20)
head(sim)
#>   id id    treatment enroll_time       tte calendar_time event
#> 1  1  1 Experimental  0.09259987 1.6537641     1.7463640     1
#> 2  1  1 Experimental  0.09259987 2.0000000     2.0925999     0
#> 3  2  2 Experimental  0.12139553 2.0000000     2.1213955     0
#> 4  3  3      Control  0.13454967 1.8086915     1.9432412     1
#> 5  3  3      Control  0.13454967 2.0000000     2.1345497     0
#> 6  4  4      Control  0.16658556 0.1134942     0.2800798     1
```
