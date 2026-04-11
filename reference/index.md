# Package index

## Sample size calculation

Plan fixed negative binomial designs for recurrent-event outcomes.

- [`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md)
  : Sample size calculation for negative binomial outcomes
- [`print(`*`<sample_size_nbinom_result>`*`)`](https://keaven.github.io/gsDesignNB/reference/print.sample_size_nbinom_result.md)
  : Print method for sample_size_nbinom_result objects
- [`summary(`*`<sample_size_nbinom_result>`*`)`](https://keaven.github.io/gsDesignNB/reference/summary.sample_size_nbinom_result.md)
  : Summary for sample_size_nbinom_result objects
- [`print(`*`<sample_size_nbinom_summary>`*`)`](https://keaven.github.io/gsDesignNB/reference/print.sample_size_nbinom_summary.md)
  : Print method for sample_size_nbinom_summary objects

## Simulation

Simulate recurrent-event data, seasonal intensity patterns, and NB
monitoring scenarios.

- [`nb_sim()`](https://keaven.github.io/gsDesignNB/reference/nb_sim.md)
  : Simulate recurrent events with fixed follow-up
- [`nb_sim_seasonal()`](https://keaven.github.io/gsDesignNB/reference/nb_sim_seasonal.md)
  : Simulate recurrent events with seasonal rates
- [`sim_gs_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_gs_nbinom.md)
  : Simulate group sequential clinical trial for negative binomial
  outcomes
- [`sim_ssr_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_ssr_nbinom.md)
  : Simulate adaptive group sequential trials with sample size
  re-estimation
- [`check_gs_bound()`](https://keaven.github.io/gsDesignNB/reference/check_gs_bound.md)
  : Check group sequential bounds
- [`summarize_gs_sim()`](https://keaven.github.io/gsDesignNB/reference/summarize_gs_sim.md)
  : Summarize group sequential simulation results
- [`summarize_ssr_sim()`](https://keaven.github.io/gsDesignNB/reference/summarize_ssr_sim.md)
  : Summarize adaptive SSR simulation results

## Analysis

Analyze interim negative binomial data and estimate information for
monitoring or SSR.

- [`blinded_ssr()`](https://keaven.github.io/gsDesignNB/reference/blinded_ssr.md)
  : Blinded sample size re-estimation for recurrent events
- [`unblinded_ssr()`](https://keaven.github.io/gsDesignNB/reference/unblinded_ssr.md)
  : Unblinded sample size re-estimation for recurrent events
- [`calculate_blinded_info()`](https://keaven.github.io/gsDesignNB/reference/calculate_blinded_info.md)
  : Calculate blinded statistical information
- [`estimate_nb_mom()`](https://keaven.github.io/gsDesignNB/reference/estimate_nb_mom.md)
  : Method of Moments Estimation for Negative Binomial Parameters
- [`mutze_test()`](https://keaven.github.io/gsDesignNB/reference/mutze_test.md)
  [`print(`*`<mutze_test>`*`)`](https://keaven.github.io/gsDesignNB/reference/mutze_test.md)
  : Wald test for treatment effect using negative binomial model (Mutze
  et al.)
- [`cut_data_by_date()`](https://keaven.github.io/gsDesignNB/reference/cut_data_by_date.md)
  : Cut simulated trial data at a calendar date
- [`cut_completers()`](https://keaven.github.io/gsDesignNB/reference/cut_completers.md)
  : Cut data for completers analysis
- [`cut_date_for_completers()`](https://keaven.github.io/gsDesignNB/reference/cut_date_for_completers.md)
  : Find calendar date for target completer count
- [`compute_info_at_time()`](https://keaven.github.io/gsDesignNB/reference/compute_info_at_time.md)
  : Compute statistical information at analysis time
- [`get_analysis_date()`](https://keaven.github.io/gsDesignNB/reference/get_analysis_date.md)
  : Find calendar date for target event count
- [`get_cut_date()`](https://keaven.github.io/gsDesignNB/reference/get_cut_date.md)
  : Determine analysis date based on criteria

## Group sequential design

Extend negative binomial designs to group sequential monitoring via
gsDesign boundaries.

- [`gsNBCalendar()`](https://keaven.github.io/gsDesignNB/reference/gsNBCalendar.md)
  : Group sequential design for negative binomial outcomes
- [`toInteger()`](https://keaven.github.io/gsDesignNB/reference/toInteger.md)
  : Convert group sequential design to integer sample sizes
- [`summary(`*`<gsNB>`*`)`](https://keaven.github.io/gsDesignNB/reference/summary.gsNB.md)
  : Summary for gsNB objects
- [`print(`*`<gsNBsummary>`*`)`](https://keaven.github.io/gsDesignNB/reference/print.gsNBsummary.md)
  : Print method for gsNBsummary objects

## Documentation site

Preview the pkgdown website locally over HTTP (avoids unstyled `file://`
pages).

- [`preview_pkgdown_site()`](https://keaven.github.io/gsDesignNB/reference/preview_pkgdown_site.md)
  : Preview built pkgdown site in the browser

## Interactive prototype

Launch the optional Shiny explorer for adaptive negative binomial SSR
scenarios.

- [`run_ssr_shiny()`](https://keaven.github.io/gsDesignNB/reference/run_ssr_shiny.md)
  : Launch the SSR Shiny prototype

## gsDesign re-exports

Spending functions and utilities re-exported from the gsDesign package
for convenience.

- [`reexports`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`gsDesign`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfBetaDist`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfCauchy`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfExponential`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfExtremeValue`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfExtremeValue2`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfGapped`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfHSD`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfLDOF`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfLDPocock`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfLinear`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfLogistic`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfNormal`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfPoints`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfPower`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfStep`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfTDist`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfTrimmed`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfTruncated`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfXG1`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfXG2`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`sfXG3`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  [`gsBoundSummary`](https://keaven.github.io/gsDesignNB/reference/reexports.md)
  : Objects exported from other packages
