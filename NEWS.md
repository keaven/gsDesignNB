# gsDesignNB 0.2.1

## Sample size and power

- `sample_size_nbinom()` computes sample size or power for fixed designs with
  two treatment groups. Supports piecewise accrual, dropout, maximum follow-up,
  and event gaps. Implements the Zhu and Lakkis (2014) and
  Friede and Schmidli (2010) methods.

## Group sequential designs

- `gsNBCalendar()` creates group sequential designs for negative binomial
  outcomes by adding calendar-time analysis schedules compatible with gsDesign.
  Inherits from both `gsDesign` and `sample_size_nbinom_result` classes.
- `compute_info_at_time()` computes statistical information for the log rate
  ratio at a given analysis time, accounting for staggered enrollment.
- `toInteger()` rounds sample sizes in a group sequential design to integers
  while respecting the randomization ratio.

## Simulation

- `nb_sim()` simulates recurrent events for trials with piecewise constant
  enrollment, exponential failure rates, and piecewise exponential dropout.
  Supports negative binomial overdispersion via gamma frailty and event gaps.
- `nb_sim_seasonal()` simulates recurrent events where event rates vary by
  season (Spring, Summer, Fall, Winter).

## Interim data handling

- `cut_data_by_date()` censors follow-up at a specified calendar time and
  aggregates events per subject, adjusting for event gaps.
- `get_analysis_date()` finds the calendar time at which a target event count
  is reached.
- `cut_completers()` subsets data to subjects randomized by a specified date.
- `cut_date_for_completers()` finds the calendar time at which a target number
  of subjects have completed their follow-up.

## Statistical inference

- `mutze_test()` fits a negative binomial (or Poisson) log-rate model and
  performs a Wald test for the treatment effect, following Mütze et al. (2019).

## Blinded sample size re-estimation

- `blinded_ssr()` estimates blinded dispersion and event rate from interim data
  and re-calculates sample size to maintain power, following Friede and
  Schmidli (2010).
- `calculate_blinded_info()` estimates blinded statistical information for
  the log rate ratio from aggregated interim data.

## Re-exports from gsDesign

- Re-exports `gsDesign()`, `gsBoundSummary()`, and common spending functions
  (`sfHSD()`, `sfLDOF()`, `sfLDPocock()`, and more) for convenience.
