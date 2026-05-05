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
  [`nb_sim()`](https://keaven.github.io/gsDesignNB/reference/nb_sim.md).

- target_completers:

  Integer. The target number of completers.

## Value

Numeric. The calendar date when `target_completers` is achieved. If the
dataset contains fewer than `target_completers` completers, returns the
maximum calendar time in the dataset and prints a message.
