# Changelog

## gsDesignNB 0.3.1 (development)

### Robust NB fallback: method-of-moments replaces Poisson under genuine overdispersion

- [`mutze_test()`](https://keaven.github.io/gsDesignNB/reference/mutze_test.md),
  [`calculate_blinded_info()`](https://keaven.github.io/gsDesignNB/reference/calculate_blinded_info.md),
  and
  [`unblinded_ssr()`](https://keaven.github.io/gsDesignNB/reference/unblinded_ssr.md)
  now fall back to method-of-moments (MoM) estimation via
  [`estimate_nb_mom()`](https://keaven.github.io/gsDesignNB/reference/estimate_nb_mom.md)
  when the maximum likelihood negative binomial fit does not converge or
  returns an extreme-overdispersion shape estimate. Previously, the
  Poisson fallback was used in both “near-Poisson” and
  “extreme-overdispersion” regimes; the latter is anti-conservative
  because the Poisson variance underestimates the true NB variance under
  genuine overdispersion. The MoM fallback computes the Wald standard
  error from the observed Fisher information formula \mathcal{I} =
  1/(1/W_1 + 1/W_2) with W_g = \sum_i \mu\_{g,i}/(1 +
  \hat{k}\mu\_{g,i}), preserving the NB variance structure without
  requiring ML convergence.
- [`mutze_test()`](https://keaven.github.io/gsDesignNB/reference/mutze_test.md)
  gains a `mom_threshold` argument (default `20`, corresponding to
  \hat{k} \> 20) that controls when the MoM branch is triggered. The
  existing `poisson_threshold` default is reduced from `1000` to `50`
  (\hat{k} \< 0.02) since NB and Poisson Wald standard errors are
  numerically indistinguishable at that point.
- All three functions now return an additional `fallback` element in
  their output (`"ml"`, `"mom"`, or `"poisson"`) so that downstream
  simulation engines can record which estimator was used at each
  interim.

## gsDesignNB 0.3.0

### Sample size methodology deep dive ([\#14](https://github.com/keaven/gsDesignNB/issues/14))

#### Jensen’s inequality correction for event gaps

- Applied a second-order Taylor correction to the effective rate formula
  when both dispersion (k \> 0) and event gap (g \> 0) are present. The
  naive formula \lambda/(1+\lambda g) overestimates the population-level
  effective rate due to Jensen’s inequality (subject-level frailty makes
  f(x)=x/(1+xg) concave). The corrected formula is \lambda\_{\text{eff}}
  \approx \frac{\lambda}{1+\lambda g}(1 - k\lambda g/(1+\lambda g)^2).
- Correction applied in both
  [`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md)
  and
  [`compute_info_at_time()`](https://keaven.github.io/gsDesignNB/reference/compute_info_at_time.md).
- Simulation study across multiple scenarios (10,000 replicates each)
  confirms the corrected design maintains nominal or conservative power,
  while the naive formula increasingly underpowers as k and g grow.

#### Documentation improvements

- Restructured the `sample-size-nbinom` vignette with a consistent
  notation table and comparison of Zhu-Lakkis, Friede-Schmidli, and
  Mutze et al. methods.
- Expanded average exposure derivation covering no-dropout, exponential
  dropout, max follow-up truncation, the Q variance inflation factor,
  and event gaps.
- Added a statistical information section covering per-subject Fisher
  information, total information, blinded vs unblinded estimation (ML
  and MoM), and the connection to sample size.
- Added simulation verification of average exposure in the vignette.
- Improved roxygen2 documentation for
  [`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md),
  [`calculate_blinded_info()`](https://keaven.github.io/gsDesignNB/reference/calculate_blinded_info.md),
  and
  [`compute_info_at_time()`](https://keaven.github.io/gsDesignNB/reference/compute_info_at_time.md)
  with `@details` sections, consistent notation, and cross-references.
- Updated `verification-simulation` vignette with a scenario sweep table
  and a discussion of why the correction is preferred despite partial
  cancellation with model-based SE bias.

#### Other

- Added DOIs to bibliography entries; added Schneider et al. (2013)
  reference.
- Switched pkgdown math rendering from MathJax (CDN) to KaTeX (bundled).

## gsDesignNB 0.2.6

CRAN release: 2026-02-16

- First CRAN release.

## gsDesignNB 0.2.5

- Added `rr0` parameter to
  [`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md)
  and
  [`blinded_ssr()`](https://keaven.github.io/gsDesignNB/reference/blinded_ssr.md)
  to support non-inferiority and super-superiority testing.
- Changed default `event_gap` to 0 in
  [`nb_sim()`](https://keaven.github.io/gsDesignNB/reference/nb_sim.md).

## gsDesignNB 0.2.4

- Fix
  [`cut_date_for_completers()`](https://keaven.github.io/gsDesignNB/reference/cut_date_for_completers.md)
  to support
  [`nb_sim_seasonal()`](https://keaven.github.io/gsDesignNB/reference/nb_sim_seasonal.md)
  output (no `tte` column).
- Correct
  [`calculate_blinded_info()`](https://keaven.github.io/gsDesignNB/reference/calculate_blinded_info.md)
  blinded information calculation to use subject-level exposure.

## gsDesignNB 0.2.3

- Fix
  [`toInteger.gsNB()`](https://keaven.github.io/gsDesignNB/reference/toInteger.md)
  to avoid unintended power changes by correctly recomputing information
  with `max_followup`, preserving `delta1`, and improving ratio-aware
  integer rounding.
- Vignette updates and documentation fixes.

## gsDesignNB 0.2.2

### Sample size and power

- [`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md)
  computes sample size or power for fixed designs with two treatment
  groups. Supports piecewise accrual, exponential dropout, maximum
  follow-up, and event gaps. Implements the Zhu and Lakkis (2014) and
  Friede and Schmidli (2010) methods.

### Group sequential designs

- [`gsNBCalendar()`](https://keaven.github.io/gsDesignNB/reference/gsNBCalendar.md)
  creates group sequential designs for negative binomial outcomes,
  optionally attaching calendar-time analysis schedules (via
  `analysis_times`) compatible with gsDesign. Inherits from both
  `gsDesign` and `sample_size_nbinom_result` classes.
- [`compute_info_at_time()`](https://keaven.github.io/gsDesignNB/reference/compute_info_at_time.md)
  computes statistical information for the log rate ratio at a given
  analysis time, accounting for staggered enrollment.
- [`toInteger()`](https://keaven.github.io/gsDesignNB/reference/toInteger.md)
  rounds sample sizes in a group sequential design to integers while
  respecting the randomization ratio.

### Simulation

- [`nb_sim()`](https://keaven.github.io/gsDesignNB/reference/nb_sim.md)
  simulates recurrent events for trials with piecewise constant
  enrollment, exponential failure rates, and piecewise exponential
  dropout. Supports negative binomial overdispersion via gamma frailty
  and event gaps.
- [`nb_sim_seasonal()`](https://keaven.github.io/gsDesignNB/reference/nb_sim_seasonal.md)
  simulates recurrent events where event rates vary by season (Spring,
  Summer, Fall, Winter).
- Group sequential simulation helpers:
  [`sim_gs_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_gs_nbinom.md)
  runs repeated simulations with flexible cut rules via
  [`get_cut_date()`](https://keaven.github.io/gsDesignNB/reference/get_cut_date.md),
  [`check_gs_bound()`](https://keaven.github.io/gsDesignNB/reference/check_gs_bound.md)
  updates spending bounds based on observed information, and
  [`summarize_gs_sim()`](https://keaven.github.io/gsDesignNB/reference/summarize_gs_sim.md)
  summarizes operating characteristics across analyses.

### Interim data handling

- [`cut_data_by_date()`](https://keaven.github.io/gsDesignNB/reference/cut_data_by_date.md)
  censors follow-up at a specified calendar time and aggregates events
  per subject, adjusting for event gaps.
- [`get_analysis_date()`](https://keaven.github.io/gsDesignNB/reference/get_analysis_date.md)
  finds the calendar time at which a target event count is reached.
- [`cut_completers()`](https://keaven.github.io/gsDesignNB/reference/cut_completers.md)
  subsets data to subjects randomized by a specified date.
- [`cut_date_for_completers()`](https://keaven.github.io/gsDesignNB/reference/cut_date_for_completers.md)
  finds the calendar time at which a target number of subjects have
  completed their follow-up.

### Statistical inference

- [`mutze_test()`](https://keaven.github.io/gsDesignNB/reference/mutze_test.md)
  fits a negative binomial (or Poisson) log-rate model and performs a
  Wald test for the treatment effect, following Mütze et al. (2019).

### Blinded sample size re-estimation

- [`blinded_ssr()`](https://keaven.github.io/gsDesignNB/reference/blinded_ssr.md)
  estimates blinded dispersion and event rate from interim data and
  re-calculates sample size to maintain power, following Friede and
  Schmidli (2010).
- [`calculate_blinded_info()`](https://keaven.github.io/gsDesignNB/reference/calculate_blinded_info.md)
  estimates blinded statistical information for the log rate ratio from
  aggregated interim data.

### Re-exports from gsDesign

- Re-exports
  [`gsDesign()`](https://keaven.github.io/gsDesign/reference/gsDesign.html),
  [`gsBoundSummary()`](https://keaven.github.io/gsDesign/reference/gsBoundSummary.html),
  and common spending functions
  ([`sfHSD()`](https://keaven.github.io/gsDesign/reference/sfHSD.html),
  [`sfLDOF()`](https://keaven.github.io/gsDesign/reference/sfLDOF.html),
  [`sfLDPocock()`](https://keaven.github.io/gsDesign/reference/sfLDOF.html),
  and more) for convenience.
