# Find calendar date for target completer count

Finds the calendar time (since start of randomization) at which a
specified number of subjects have completed their follow-up.

## Usage

``` r
cut_date_for_completers(data, target_completers)
```

## Arguments

- data:

  A data frame of simulated data, typically from
  [`nb_sim()`](https://keaven.github.io/gsDesignNB/reference/nb_sim.md)
  or
  [`nb_sim_seasonal()`](https://keaven.github.io/gsDesignNB/reference/nb_sim_seasonal.md).

- target_completers:

  Integer. The target number of completers.

## Value

Numeric. The calendar date when `target_completers` is achieved. If the
dataset contains fewer than `target_completers` completers, returns the
maximum calendar time in the dataset and prints a message.

## Examples

``` r
enroll_rate <- data.frame(rate = 20 / (5 / 12), duration = 5 / 12)
fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3))
dropout_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0.1, 0.05), duration = c(100, 100)
)
sim <- nb_sim(enroll_rate, fail_rate, dropout_rate, max_followup = 2, n = 20)
cut_date_for_completers(sim, target_completers = 5)
#> [1] 2.065678
```
