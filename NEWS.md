# gsDesignNB 0.3.2 (development)

## Score-test sizing and inference guidance

- Expanded the paper, sample-size vignette, and score-vs-Wald simulation
  vignette with recommendations for when the Zhu-Lakkis / Friede-Schmidli /
  Mutze Wald sample-size formula is appropriate and when to use score-test
  sizing with separate null and alternative variance factors.
- Added an `ai-skills` vignette demonstrating how the package `SKILL.md` and
  generated `llms.txt` can guide package-native score-test and simulation
  workflows without replacing statistical review.
- Refined score-test recommendations to distinguish the final analysis test
  from the sample-size formula: Wald/Zhu-Lakkis sizing remains a useful
  practical baseline and may provide a power margin when paired with the score
  test in fixed-design superiority settings, while SSR designs should verify
  any starting-size margin inside the planned adaptation rule.
- Added a targeted SSR starting-size sensitivity cache comparing Wald- and
  score-sized starting designs under the score final test in a low-event stress
  setting. Both starts preserved near-nominal score-test Type I error, and the
  larger Wald-sized start did not produce a clear SSR power advantage.
- Documented the cached score-vs-Wald simulation results: the Wald test was
  mildly anti-conservative in several finite-sample scenarios, while the score
  test preserved Type I error more conservatively and should be paired with
  simulation-based power checks.
- Corrected the reported `variance_null` field from `sample_size_nbinom()` so
  it is on the same final-analysis scale as `variance`; the score sample-size
  calculation itself was already using the null variance factor correctly.

## Group sequential simulation vignette

- Corrected the group sequential simulation vignette and cached results so that
  simulations use the rounded final group sequential sample size, monthly
  dropout hazard, and helper functions for boundary checking and summary.
- `toInteger.gsNB()` now preserves calendar-time enrollment quantities at
  interim analyses and recomputes expected events, exposures, and information
  after rounding the final sample size.
- `summarize_gs_sim()` now reports optional sample-size and exposure summaries
  when available and uses finite trimmed means for information estimates.
- Updated SSR simulation reporting in the manuscript and SSR simulation article
  to use the production score-test cache: 3,600 replicates per power scenario,
  20,000 per RR = 1 main-grid scenario, 1,000 per RR > 1 scenario, and 20,000
  per dispersion/test-statistic cell for the non-binding Type I tables.
- Updated the SSR simulation study to compare Wald and score tests at the same
  nominal one-sided alpha of 0.025, and aligned the SSR power simulations with
  the score final-analysis recommendation.
- Added a checkpointed production generator for the SSR score-test cache so the
  long-running power and Type I simulations can be resumed from saved chunks.
- Updated the blinded-information diagnostics article to distinguish historical
  raw ML pathologies from the current MoM fallback behavior.

## Missing data and imputation documentation

- Added documentation clarifying that the primary recurrent-event analysis is an
  observed-exposure negative binomial likelihood analysis. Under ignorable
  censoring / MAR assumptions, partially followed subjects contribute their
  observed events and exposure, and multiple imputation is not required simply
  because follow-up is censored.
- Expanded MNAR sensitivity-analysis guidance for recurrent-event endpoints,
  including the need to preserve censoring reason and planned remaining follow-up
  before post-dropout outcomes can be imputed.
- Added Keene, Roger, Hartley, and Kenward (2014) as the recurrent-event
  controlled-imputation reference for de facto / reference-based sensitivity
  analyses, and aligned the manuscript, slides, and simulation vignette around
  that framing.
- Updated the paper and slide materials to cite package articles as executable
  supplementary material for the fuller simulation grids and reporting.

# gsDesignNB 0.3.1

## Robust NB fallback: method-of-moments replaces Poisson under genuine overdispersion

- `mutze_test()`, `calculate_blinded_info()`, and `unblinded_ssr()` now fall
  back to method-of-moments (MoM) estimation via `estimate_nb_mom()` when the
  maximum likelihood negative binomial fit does not converge or returns an
  extreme-overdispersion shape estimate. Previously, the Poisson fallback was
  used in both "near-Poisson" and "extreme-overdispersion" regimes; the
  latter is anti-conservative because the Poisson variance underestimates the
  true NB variance under genuine overdispersion. The MoM fallback computes
  the Wald standard error from the observed Fisher information formula
  $\mathcal{I} = 1/(1/W_1 + 1/W_2)$ with
  $W_g = \sum_i \mu_{g,i}/(1 + \hat{k}\mu_{g,i})$, preserving the NB variance
  structure without requiring ML convergence.
- `mutze_test()` gains a `mom_threshold` argument (default `20`, corresponding
  to $\hat{k} > 20$) that controls when the MoM branch is triggered. The
  existing `poisson_threshold` default is reduced from `1000` to `50`
  ($\hat{k} < 0.02$) since NB and Poisson Wald standard errors are
  numerically indistinguishable at that point.
- All three functions now return an additional `fallback` element in their
  output (`"ml"`, `"mom"`, or `"poisson"`) so that downstream simulation
  engines can record which estimator was used at each interim.

# gsDesignNB 0.3.0

## Sample size methodology deep dive (#14)

### Jensen's inequality correction for event gaps

- Applied a second-order Taylor correction to the effective rate formula when
  both dispersion ($k > 0$) and event gap ($g > 0$) are present. The naive
  formula $\lambda/(1+\lambda g)$ overestimates the population-level effective
  rate due to Jensen's inequality (subject-level frailty makes $f(x)=x/(1+xg)$
  concave). The corrected formula is
  $\lambda_{\text{eff}} \approx \frac{\lambda}{1+\lambda g}(1 - k\lambda g/(1+\lambda g)^2)$.
- Correction applied in both `sample_size_nbinom()` and `compute_info_at_time()`.
- Simulation study across multiple scenarios (10,000 replicates each) confirms
  the corrected design maintains nominal or conservative power, while the naive
  formula increasingly underpowers as $k$ and $g$ grow.

### Documentation improvements

- Restructured the `sample-size-nbinom` vignette with a consistent notation
  table and comparison of Zhu-Lakkis, Friede-Schmidli, and Mutze et al. methods.
- Expanded average exposure derivation covering no-dropout, exponential dropout,
  max follow-up truncation, the $Q$ variance inflation factor, and event gaps.
- Added a statistical information section covering per-subject Fisher
  information, total information, blinded vs unblinded estimation (ML and MoM),
  and the connection to sample size.
- Added simulation verification of average exposure in the vignette.
- Improved roxygen2 documentation for `sample_size_nbinom()`,
  `calculate_blinded_info()`, and `compute_info_at_time()` with `@details`
  sections, consistent notation, and cross-references.
- Updated `verification-simulation` vignette with a scenario sweep table and
  a discussion of why the correction is preferred despite partial cancellation
  with model-based SE bias.

### Other

- Added DOIs to bibliography entries; added Schneider et al. (2013) reference.
- Switched pkgdown math rendering from MathJax (CDN) to KaTeX (bundled).

# gsDesignNB 0.2.6

- First CRAN release.

# gsDesignNB 0.2.5

- Added `rr0` parameter to `sample_size_nbinom()` and `blinded_ssr()` to support non-inferiority and super-superiority testing.
- Changed default `event_gap` to 0 in `nb_sim()`.

# gsDesignNB 0.2.4

- Fix `cut_date_for_completers()` to support `nb_sim_seasonal()` output (no `tte` column).
- Correct `calculate_blinded_info()` blinded information calculation to use subject-level exposure.

# gsDesignNB 0.2.3

- Fix `toInteger.gsNB()` to avoid unintended power changes by correctly recomputing
  information with `max_followup`, preserving `delta1`, and improving ratio-aware
  integer rounding.
- Vignette updates and documentation fixes.

# gsDesignNB 0.2.2

## Sample size and power

- `sample_size_nbinom()` computes sample size or power for fixed designs with
  two treatment groups. Supports piecewise accrual, exponential dropout,
  maximum follow-up, and event gaps. Implements the Zhu and Lakkis (2014) and
  Friede and Schmidli (2010) methods.

## Group sequential designs

- `gsNBCalendar()` creates group sequential designs for negative binomial
  outcomes, optionally attaching calendar-time analysis schedules (via
  `analysis_times`) compatible with gsDesign.
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
- Group sequential simulation helpers: `sim_gs_nbinom()` runs repeated
  simulations with flexible cut rules via `get_cut_date()`,
  `check_gs_bound()` updates spending bounds based on observed information,
  and `summarize_gs_sim()` summarizes operating characteristics across analyses.

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
