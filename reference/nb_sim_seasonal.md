# Simulate recurrent events with seasonal rates

Simulates recurrent events where event rates depend on the season.

## Usage

``` r
nb_sim_seasonal(
  enroll_rate,
  fail_rate,
  dropout_rate = NULL,
  max_followup = NULL,
  randomization_start_date = NULL,
  n = NULL,
  block = c(rep("Control", 2), rep("Experimental", 2))
)
```

## Arguments

- enroll_rate:

  A data frame with columns `rate` and `duration`.

- fail_rate:

  A data frame with columns `treatment`, `season`, `rate`, and
  optionally `dispersion`. Seasons should be "Spring", "Summer", "Fall",
  "Winter".

- dropout_rate:

  A data frame with columns `treatment`, `rate`, `duration`.

- max_followup:

  Numeric. Max follow-up duration (years).

- randomization_start_date:

  Date. Start of randomization.

- n:

  Integer. Total sample size.

- block:

  Character vector for block randomization.

## Value

A data frame of class `nb_sim_seasonal` with columns: `id`, `treatment`,
`season`, `enroll_time`, `start`, `end`, `event`, `calendar_start`,
`calendar_end`. Rows represent intervals of risk or events. `event=1`
indicates an event at `end`. `event=0` indicates censoring or end of a
seasonal interval at `end`.
