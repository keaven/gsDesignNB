# Find calendar date for target event count

Finds the calendar time (since start of randomization) at which a
specified total number of events is reached in the simulated dataset.

## Usage

``` r
get_analysis_date(data, planned_events, event_gap = 5/365.25)
```

## Arguments

- data:

  A data frame of simulated data, typically from
  [`nb_sim()`](https://keaven.github.io/gsDesignNB/reference/nb_sim.md).

- planned_events:

  Integer. The target number of events.

- event_gap:

  Gap duration after each event during which no new events are counted.
  Can be a numeric value (default `5 / 365.25`) or a function returning
  a numeric value.

## Value

Numeric. The calendar date when `planned_events` is achieved. If the
dataset contains fewer than `planned_events`, returns the maximum
calendar time in the dataset and prints a message.

## Examples

``` r
enroll_rate <- data.frame(rate = 20 / (5 / 12), duration = 5 / 12)
fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3))
dropout_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0.1, 0.05), duration = c(100, 100)
)
sim <- nb_sim(enroll_rate, fail_rate, dropout_rate, max_followup = 2, n = 40)
get_analysis_date(sim, planned_events = 15)
#> [1] 0.9853362
```
