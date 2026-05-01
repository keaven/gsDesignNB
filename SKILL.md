# gsDesignNB — Comprehensive AI Skill

> **Last updated:** 2026-04-30 — matches package version 0.3.2

## Package purpose

**gsDesignNB** is an R package for planning, monitoring, and adapting
clinical trials whose primary endpoint is a recurrent event count
analyzed with a negative binomial (NB) model. It extends the
**gsDesign** ecosystem to handle the special features of NB designs:
information that depends on exposure and dispersion (not just event
counts), protocol event gaps, calendar-time interim looks, and
blinded/unblinded sample size re-estimation.

## Statistical model

Conditional on a subject-level gamma frailty, events follow a Poisson
process with constant rate. Marginally, the event count Y_i for subject
i with exposure t_i is:

- Mean: \mu_i = \lambda_g t_i
- Variance: \mu_i + k \mu_i^2

where \lambda_g is the group rate and k is the dispersion parameter. The
treatment effect is the log rate ratio \theta = \log(\lambda_2 /
\lambda_1), estimated via NB GLM with log link and exposure offset.

Fisher information per subject: \mathcal{I}\_i = \mu_i / (1 + k \mu_i).

### Wald vs score test

[`mutze_test()`](https://keaven.github.io/gsDesignNB/reference/mutze_test.md)
supports `test_type = "wald"` (default) and `test_type = "score"`:

- **Wald test:** Fits arm-specific NB GLMs; z = \hat\theta /
  \text{SE}(\hat\theta).
- **Score test:** Fits a pooled null NB GLM via
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html); score
  statistic U = \sum (y_i - \hat\mu_i) / (1 + \hat k_0 \hat\mu_i);
  Fisher information I_0 = W_1 W_2 / (W_1 + W_2); z = U / \sqrt{I_0}.

The score test is slightly more conservative under H_0 and better
controls Type I error in small-sample / low-rate scenarios.

### Two-variance sample size formula

[`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md)
supports `test_type = "wald"` (default) and `test_type = "score"`:

- **Wald (single variance):** n_1 = (z\_\alpha + z\_\beta)^2 V_1 /
  \theta^2
- **Score (two variance):** n_1 = (z\_\alpha \sqrt{V_0} + z\_\beta
  \sqrt{V_1})^2 / \theta^2

where V_1 = (1/\mu_1 + k_1)(1 + 1/r) + (1/\mu_2 + k_2)/r is the
alternative-hypothesis variance, and V_0 = (1/\mu_0 + k_0)(1 + 1/r) uses
a pooled rate \lambda_0 = (\lambda_1 + r \lambda_2) / (1 + r) under H_0.

The score formula is matched to the score test’s null-variance reference
distribution, but it is not automatically more powerful. In the package
superiority grid, Wald/Zhu–Lakkis sizing paired with the score test
preserved Type I error and gave slightly higher score-test power than
score sizing because it retained a small sample-size margin. Compare
both sizing rules and verify Type I error and power by simulation for
the design setting.

## Core workflow

    sample_size_nbinom() → gsNBCalendar() → sim_gs_nbinom() or sim_ssr_nbinom() → summarize

1.  **Fixed design:**
    [`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md)
    computes sample size/power under planned rates, dispersion,
    piecewise accrual, exponential dropout, max follow-up, and event
    gaps.
2.  **GS design:**
    [`gsNBCalendar()`](https://keaven.github.io/gsDesignNB/reference/gsNBCalendar.md)
    converts the fixed design to a `gsNB` object (inherits `gsDesign`)
    with calendar-time analysis looks and NB-specific information.
3.  **Simulation:**
    [`sim_gs_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_gs_nbinom.md)
    (fixed GS timing) or
    [`sim_ssr_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_ssr_nbinom.md)
    (adaptive SSR) generates trial-level operating characteristics.
4.  **Summary:**
    [`check_gs_bound()`](https://keaven.github.io/gsDesignNB/reference/check_gs_bound.md) +
    [`summarize_gs_sim()`](https://keaven.github.io/gsDesignNB/reference/summarize_gs_sim.md),
    or
    [`summarize_ssr_sim()`](https://keaven.github.io/gsDesignNB/reference/summarize_ssr_sim.md).

## Function reference

### Sample size and design

| Function | Purpose |
|:---|:---|
| [`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md) | Fixed-design sample size/power; `test_type = "wald"` (V₁ only) or `"score"` (V₀+V₁) |
| [`compute_info_at_time()`](https://keaven.github.io/gsDesignNB/reference/compute_info_at_time.md) | Information at a given calendar time (internal to sample_size_nbinom) |
| [`gsNBCalendar()`](https://keaven.github.io/gsDesignNB/reference/gsNBCalendar.md) | Calendar-time group sequential design from fixed design result |
| [`toInteger()`](https://keaven.github.io/gsDesignNB/reference/toInteger.md) | Round design to integer sample sizes (S3 method for `gsNB`) |

### Data simulation

| Function | Purpose |
|:---|:---|
| [`nb_sim()`](https://keaven.github.io/gsDesignNB/reference/nb_sim.md) | Simulate recurrent events — piecewise accrual, exponential dropout, gamma-Poisson mixture |
| [`nb_sim_seasonal()`](https://keaven.github.io/gsDesignNB/reference/nb_sim_seasonal.md) | Same as [`nb_sim()`](https://keaven.github.io/gsDesignNB/reference/nb_sim.md) but with season-varying event rates |

### Interim data handling

| Function | Purpose |
|:---|:---|
| [`cut_data_by_date()`](https://keaven.github.io/gsDesignNB/reference/cut_data_by_date.md) | Censor simulated data at a calendar date; aggregate events/exposure per subject |
| [`get_cut_date()`](https://keaven.github.io/gsDesignNB/reference/get_cut_date.md) | Find analysis date from calendar, event, completer, or information targets |
| [`get_analysis_date()`](https://keaven.github.io/gsDesignNB/reference/get_analysis_date.md) | Find calendar time when target event count is reached |
| [`cut_date_for_completers()`](https://keaven.github.io/gsDesignNB/reference/cut_date_for_completers.md) | Find calendar time when target completer count is reached |
| [`cut_completers()`](https://keaven.github.io/gsDesignNB/reference/cut_completers.md) | Convenience wrapper: find completer date + cut data |

### Statistical analysis and information

| Function | Purpose |
|:---|:---|
| [`mutze_test()`](https://keaven.github.io/gsDesignNB/reference/mutze_test.md) | Wald or score test for treatment effect (`test_type = "wald"` / `"score"`) |
| [`calculate_blinded_info()`](https://keaven.github.io/gsDesignNB/reference/calculate_blinded_info.md) | Blinded dispersion/rate estimation and statistical information |
| [`estimate_nb_mom()`](https://keaven.github.io/gsDesignNB/reference/estimate_nb_mom.md) | Method-of-moments NB parameter estimation (blinded or unblinded) |

### Group sequential simulation

| Function | Purpose |
|:---|:---|
| [`sim_gs_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_gs_nbinom.md) | Monte Carlo GS trial simulation with calendar-time cuts |
| [`check_gs_bound()`](https://keaven.github.io/gsDesignNB/reference/check_gs_bound.md) | Apply GS boundaries to simulation results on chosen info scale |
| [`summarize_gs_sim()`](https://keaven.github.io/gsDesignNB/reference/summarize_gs_sim.md) | Operating characteristic summaries (power, futility, per-analysis) |

### Sample size re-estimation

| Function | Purpose |
|:---|:---|
| [`blinded_ssr()`](https://keaven.github.io/gsDesignNB/reference/blinded_ssr.md) | Blinded SSR using pooled NB fit (Friede & Schmidli 2010) |
| [`unblinded_ssr()`](https://keaven.github.io/gsDesignNB/reference/unblinded_ssr.md) | Unblinded SSR using arm-specific rate/dispersion estimates |
| [`sim_ssr_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_ssr_nbinom.md) | Full adaptive simulation engine (No adaptation / Blinded / Unblinded) |
| [`summarize_ssr_sim()`](https://keaven.github.io/gsDesignNB/reference/summarize_ssr_sim.md) | Summary of SSR simulation results |

### Utilities

| Function | Purpose |
|:---|:---|
| [`run_ssr_shiny()`](https://keaven.github.io/gsDesignNB/reference/run_ssr_shiny.md) | Interactive Shiny interface for SSR exploration |
| [`preview_pkgdown_site()`](https://keaven.github.io/gsDesignNB/reference/preview_pkgdown_site.md) | Serve docs locally with correct CSS (not [`pkgdown::preview_site()`](https://pkgdown.r-lib.org/reference/preview_site.html)) |

### Re-exported from gsDesign

[`gsDesign()`](https://keaven.github.io/gsDesign/reference/gsDesign.html),
[`gsBoundSummary()`](https://keaven.github.io/gsDesign/reference/gsBoundSummary.html),
and all `sf*()` spending functions are re-exported so users can work in
one namespace.

## Jensen’s inequality correction for event gaps

When a protocol defines a minimum gap \gamma between counted events, the
naive effective rate \lambda\_{\text{eff}} = \lambda / (1 +
\lambda\gamma) is biased upward due to Jensen’s inequality applied over
the gamma-distributed subject-level rates. The second-order correction
is:

\lambda\_{\text{eff}} \approx \frac{\lambda}{1 + \lambda\gamma}
\left(1 - \frac{k\lambda\gamma}{(1 + \lambda\gamma)^2}\right)

**When it matters:** k \geq 0.5 and gap \geq 20 days. A 64-scenario
simulation sweep showed the naive formula underpowers by 1–2.5 pp in
these regimes; the correction reduces effective-rate approximation error
from ~8% to ~1%.

**When it is negligible:** k \leq 0.1 or gap \leq 10 days.

The same gap rule flows through
[`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md),
[`nb_sim()`](https://keaven.github.io/gsDesignNB/reference/nb_sim.md),
and
[`cut_data_by_date()`](https://keaven.github.io/gsDesignNB/reference/cut_data_by_date.md)
so that planning, simulation, and interim cutting are consistent.

## Three kinds of “gaps”

| Concept | Scale | Package interface |
|:---|:---|:---|
| Analysis gap \Delta_i | Trial calendar | `analysis_times` in [`gsNBCalendar()`](https://keaven.github.io/gsDesignNB/reference/gsNBCalendar.md), search criteria in [`get_cut_date()`](https://keaven.github.io/gsDesignNB/reference/get_cut_date.md) |
| Timing heterogeneity | Accrual, dropout, follow-up | `accrual_rate/duration`, `dropout_rate`, `max_followup`, `trial_duration` |
| Event gap g | Patient time | `event_gap` argument, carried through planning → simulation → interim cut |

## SSR guidance

- [`sim_ssr_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_ssr_nbinom.md)
  supports strategies: `"No adaptation"`, `"Blinded SSR"`,
  `"Unblinded SSR"`.
- `bound_info` controls which information scale drives boundaries:
  `"unblinded_ml"` (default), `"blinded_ml"`, `"unblinded_mom"`,
  `"blinded_mom"`.
- For SSR studies, treat the final test statistic as central. The
  current production SSR study uses the score test for power simulations
  and compares Wald vs score Type I error at the same nominal one-sided
  \alpha = 0.025.
- In that non-binding Type I comparison, Wald is mildly
  anti-conservative (about 2.8%–3.1%) and score is conservative (about
  2.1%–2.3%). Prefer the score final test when Type I calibration is the
  priority, then check whether Wald sizing, a higher target power, or a
  modest information margin is needed for power.
- SSR adaptation is driven mainly by interim nuisance estimates and the
  adaptation cap; the score test reduces anti-conservative rejection
  after adaptation rather than eliminating the need to simulate adapted
  sample-size behavior.
- Blinded SSR preserves masking. Unblinded SSR may be more
  sample-efficient, but the relative behavior is design-specific and
  should be compared under the trial’s nuisance range, adaptation cap,
  and final test statistic.
- Spending is not accelerated when observed information exceeds planned:
  uses `min(planned IF, actual IF)`.
- Futility assessment is deferred until ≥ 30% of planned information.

## Time-scale checklist

Before writing or reviewing code, confirm all of these use the **same
time unit**:

- `enroll_rate$duration` and `enroll_rate$rate`
- `fail_rate$rate`
- `dropout_rate$rate` and `dropout_rate$duration`
- `trial_duration`
- `max_followup`
- Calendar cut times / `analysis_times`
- `event_gap`

If units are mixed, fix that first.

## Key arguments shared across functions

| Argument | Used in | Meaning |
|:---|:---|:---|
| `lambda1`, `lambda2` | [`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md) | Control and experimental event rates |
| `dispersion` / `k` | Multiple | NB dispersion (scalar or arm-specific vector) |
| `accrual_rate`, `accrual_duration` | [`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md) | Piecewise enrollment |
| `enroll_rate` | [`nb_sim()`](https://keaven.github.io/gsDesignNB/reference/nb_sim.md), [`sim_gs_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_gs_nbinom.md), [`sim_ssr_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_ssr_nbinom.md) | `data.frame(rate, duration)` |
| `fail_rate` | Simulation functions | `data.frame(treatment, rate, dispersion)` |
| `dropout_rate` | Multiple | `data.frame(treatment, rate, duration)` for exponential dropout |
| `max_followup` | Multiple | Maximum per-subject follow-up |
| `trial_duration` | [`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md) | Total calendar time of the trial |
| `event_gap` | [`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md), [`nb_sim()`](https://keaven.github.io/gsDesignNB/reference/nb_sim.md), [`cut_data_by_date()`](https://keaven.github.io/gsDesignNB/reference/cut_data_by_date.md) | Minimum gap between counted events |
| `rr0` | [`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md) | Null rate ratio (default 1; set \< 1 for non-inferiority) |
| `test_type` | [`mutze_test()`](https://keaven.github.io/gsDesignNB/reference/mutze_test.md), [`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md), [`sim_gs_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_gs_nbinom.md), [`sim_ssr_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_ssr_nbinom.md) | `"wald"` (default) or `"score"` |

## S3 class hierarchy

- `sample_size_nbinom_result`: output of
  [`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md);
  has [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html).
- `gsNB`: output of
  [`gsNBCalendar()`](https://keaven.github.io/gsDesignNB/reference/gsNBCalendar.md);
  inherits from `gsDesign`; has
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`toInteger()`](https://keaven.github.io/gsDesignNB/reference/toInteger.md).
- `nb_sim_data` / `nb_sim_seasonal`: output of
  [`nb_sim()`](https://keaven.github.io/gsDesignNB/reference/nb_sim.md)
  /
  [`nb_sim_seasonal()`](https://keaven.github.io/gsDesignNB/reference/nb_sim_seasonal.md);
  dispatches
  [`cut_data_by_date()`](https://keaven.github.io/gsDesignNB/reference/cut_data_by_date.md).
- `mutze_test`: output of
  [`mutze_test()`](https://keaven.github.io/gsDesignNB/reference/mutze_test.md);
  has [`print()`](https://rdrr.io/r/base/print.html).

## Vignettes

| Vignette | Topic |
|:---|:---|
| `sample-size-nbinom` | Fixed-design NB sample size planning |
| `simulation-example` | Basic recurrent-event data generation |
| `seasonal-simulation` | Seasonal rate simulation |
| `group-sequential-simulation` | Calendar-time GS design workflow |
| `ssr-example` | Short SSR workflow |
| `ssr-simulation-study` | Large SSR simulation case study |
| `completers-interim-example` | Completer-based interim looks |
| `non-inferiority-example` | Non-inferiority design |
| `verification-simulation` | Simulation verification of sample size and test |
| `score-vs-wald-simulation` | 2×2 comparison: Wald/score sizing × Wald/score test (interactive DT/plotly) |
| `ai-skills` | Demonstrates how to use `SKILL.md` and `llms.txt` to guide package-native workflows |
| `blinded-info-diagnostics` | Blinded information estimation diagnostics |
| `blinded-info-diagnostics` | Edge cases for blinded information |
| `verification-simulation` | Design-formula verification against simulation |

## Common pitfalls

- Do not re-implement SSR logic if
  [`sim_ssr_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_ssr_nbinom.md)
  covers the workflow.
- [`mutze_test()`](https://keaven.github.io/gsDesignNB/reference/mutze_test.md)
  may fall back to Poisson based on `poisson_threshold`; this is
  defensive, not a model choice.
- `event_gap` changes both counted events and at-risk exposure; it is
  not cosmetic.
- [`check_gs_bound()`](https://keaven.github.io/gsDesignNB/reference/check_gs_bound.md)
  accepts `info_col` to select the information column; use that rather
  than renaming.
- For pkgdown preview, use
  [`gsDesignNB::preview_pkgdown_site()`](https://keaven.github.io/gsDesignNB/reference/preview_pkgdown_site.md),
  not
  [`pkgdown::preview_site()`](https://pkgdown.r-lib.org/reference/preview_site.html).

## Disease-area calibration

Published NB trial parameters for reference:

| Therapeutic area | ~λ (events/month) | ~k | Gap |
|:---|:--:|:--:|:---|
| Heart failure hospitalizations | 0.3 | 0.5 | Extended hospitalization-free window |
| COPD exacerbations | 0.5–0.8 | 0.5–1.0 | ≥ 28 symptom-free days |
| MS relapses | 1–2/year | moderate–high | Variable |
| Asthma exacerbations | Variable | ~0.8 | Protocol-specific |

## Dependencies

- **Runtime:** gsDesign (≥ 3.9.0), data.table, simtrial, MASS, stats
- **Suggested:** ggplot2, dplyr, gt, DT, shiny, foreach, doFuture,
  future, testthat, knitr
