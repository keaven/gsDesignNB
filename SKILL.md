---
name: gsDesignNB
description: >
  Design, simulate, and adapt clinical trials with negative binomial recurrent event endpoints
  using the gsDesignNB R package. Use this skill when the task involves: NB sample size or power,
  event gaps and Jensen correction, calendar-time group sequential design, blinded or unblinded
  information estimation, sample size re-estimation (SSR), recurrent-event simulation (constant or
  seasonal rates), interim data cuts, completer-based analyses, or non-inferiority/super-superiority
  designs for recurrent events.
---

# gsDesignNB — Comprehensive AI Skill

## Package purpose

**gsDesignNB** is an R package for planning, monitoring, and adapting clinical trials whose primary
endpoint is a recurrent event count analyzed with a negative binomial (NB) model. It extends the
**gsDesign** ecosystem to handle the special features of NB designs: information that depends on
exposure and dispersion (not just event counts), protocol event gaps, calendar-time interim looks,
and blinded/unblinded sample size re-estimation.

## Statistical model

Conditional on a subject-level gamma frailty, events follow a Poisson process with constant rate.
Marginally, the event count $Y_i$ for subject $i$ with exposure $t_i$ is:

- Mean: $\mu_i = \lambda_g t_i$
- Variance: $\mu_i + k \mu_i^2$

where $\lambda_g$ is the group rate and $k$ is the dispersion parameter. The treatment effect is
the log rate ratio $\theta = \log(\lambda_2 / \lambda_1)$, estimated via NB GLM with log link and
exposure offset.

Fisher information per subject: $\mathcal{I}_i = \mu_i / (1 + k \mu_i)$.

## Core workflow

```
sample_size_nbinom() → gsNBCalendar() → sim_gs_nbinom() or sim_ssr_nbinom() → summarize
```

1. **Fixed design:** `sample_size_nbinom()` computes sample size/power under planned rates,
   dispersion, piecewise accrual, exponential dropout, max follow-up, and event gaps.
2. **GS design:** `gsNBCalendar()` converts the fixed design to a `gsNB` object (inherits
   `gsDesign`) with calendar-time analysis looks and NB-specific information.
3. **Simulation:** `sim_gs_nbinom()` (fixed GS timing) or `sim_ssr_nbinom()` (adaptive SSR)
   generates trial-level operating characteristics.
4. **Summary:** `check_gs_bound()` + `summarize_gs_sim()`, or `summarize_ssr_sim()`.

## Function reference

### Sample size and design

| Function | Purpose |
|:---------|:--------|
| `sample_size_nbinom()` | Fixed-design sample size/power for NB recurrent events |
| `compute_info_at_time()` | Information at a given calendar time (internal to sample_size_nbinom) |
| `gsNBCalendar()` | Calendar-time group sequential design from fixed design result |
| `toInteger()` | Round design to integer sample sizes (S3 method for `gsNB`) |

### Data simulation

| Function | Purpose |
|:---------|:--------|
| `nb_sim()` | Simulate recurrent events — piecewise accrual, exponential dropout, gamma-Poisson mixture |
| `nb_sim_seasonal()` | Same as `nb_sim()` but with season-varying event rates |

### Interim data handling

| Function | Purpose |
|:---------|:--------|
| `cut_data_by_date()` | Censor simulated data at a calendar date; aggregate events/exposure per subject |
| `get_cut_date()` | Find analysis date from calendar, event, completer, or information targets |
| `get_analysis_date()` | Find calendar time when target event count is reached |
| `cut_date_for_completers()` | Find calendar time when target completer count is reached |
| `cut_completers()` | Convenience wrapper: find completer date + cut data |

### Statistical analysis and information

| Function | Purpose |
|:---------|:--------|
| `mutze_test()` | Wald test for treatment effect (NB GLM, Poisson fallback) |
| `calculate_blinded_info()` | Blinded dispersion/rate estimation and statistical information |
| `estimate_nb_mom()` | Method-of-moments NB parameter estimation (blinded or unblinded) |

### Group sequential simulation

| Function | Purpose |
|:---------|:--------|
| `sim_gs_nbinom()` | Monte Carlo GS trial simulation with calendar-time cuts |
| `check_gs_bound()` | Apply GS boundaries to simulation results on chosen info scale |
| `summarize_gs_sim()` | Operating characteristic summaries (power, futility, per-analysis) |

### Sample size re-estimation

| Function | Purpose |
|:---------|:--------|
| `blinded_ssr()` | Blinded SSR using pooled NB fit (Friede & Schmidli 2010) |
| `unblinded_ssr()` | Unblinded SSR using arm-specific rate/dispersion estimates |
| `sim_ssr_nbinom()` | Full adaptive simulation engine (No adaptation / Blinded / Unblinded) |
| `summarize_ssr_sim()` | Summary of SSR simulation results |

### Utilities

| Function | Purpose |
|:---------|:--------|
| `run_ssr_shiny()` | Interactive Shiny interface for SSR exploration |
| `preview_pkgdown_site()` | Serve docs locally with correct CSS (not `pkgdown::preview_site()`) |

### Re-exported from gsDesign

`gsDesign()`, `gsBoundSummary()`, and all `sf*()` spending functions are re-exported so users
can work in one namespace.

## Jensen's inequality correction for event gaps

When a protocol defines a minimum gap $\gamma$ between counted events, the naive effective rate
$\lambda_{\text{eff}} = \lambda / (1 + \lambda\gamma)$ is biased upward due to Jensen's inequality
applied over the gamma-distributed subject-level rates. The second-order correction is:

$$
\lambda_{\text{eff}} \approx \frac{\lambda}{1 + \lambda\gamma}
\left(1 - \frac{k\lambda\gamma}{(1 + \lambda\gamma)^2}\right)
$$

**When it matters:** $k \geq 0.5$ and gap $\geq 20$ days. A 64-scenario simulation sweep showed
the naive formula underpowers by 1–2.5 pp in these regimes; the correction reduces effective-rate
approximation error from ~8% to ~1%.

**When it is negligible:** $k \leq 0.1$ or gap $\leq 10$ days.

The same gap rule flows through `sample_size_nbinom()`, `nb_sim()`, and `cut_data_by_date()` so
that planning, simulation, and interim cutting are consistent.

## Three kinds of "gaps"

| Concept | Scale | Package interface |
|:--------|:------|:------------------|
| Analysis gap $\Delta_i$ | Trial calendar | `analysis_times` in `gsNBCalendar()`, search criteria in `get_cut_date()` |
| Timing heterogeneity | Accrual, dropout, follow-up | `accrual_rate/duration`, `dropout_rate`, `max_followup`, `trial_duration` |
| Event gap $g$ | Patient time | `event_gap` argument, carried through planning → simulation → interim cut |

## SSR guidance

- `sim_ssr_nbinom()` supports strategies: `"No adaptation"`, `"Blinded SSR"`, `"Unblinded SSR"`.
- `bound_info` controls which information scale drives boundaries: `"unblinded_ml"` (default),
  `"blinded_ml"`, `"unblinded_mom"`, `"blinded_mom"`.
- Blinded SSR preserves masking but requests larger sample sizes; unblinded SSR is more efficient
  when operationally acceptable.
- Spending is not accelerated when observed information exceeds planned: uses `min(planned IF, actual IF)`.
- Futility assessment is deferred until ≥ 30% of planned information.

## Time-scale checklist

Before writing or reviewing code, confirm all of these use the **same time unit**:

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
|:---------|:--------|:--------|
| `lambda1`, `lambda2` | `sample_size_nbinom()` | Control and experimental event rates |
| `dispersion` / `k` | Multiple | NB dispersion (scalar or arm-specific vector) |
| `accrual_rate`, `accrual_duration` | `sample_size_nbinom()` | Piecewise enrollment |
| `enroll_rate` | `nb_sim()`, `sim_gs_nbinom()`, `sim_ssr_nbinom()` | `data.frame(rate, duration)` |
| `fail_rate` | Simulation functions | `data.frame(treatment, rate, dispersion)` |
| `dropout_rate` | Multiple | `data.frame(treatment, rate, duration)` for exponential dropout |
| `max_followup` | Multiple | Maximum per-subject follow-up |
| `trial_duration` | `sample_size_nbinom()` | Total calendar time of the trial |
| `event_gap` | `sample_size_nbinom()`, `nb_sim()`, `cut_data_by_date()` | Minimum gap between counted events |
| `rr0` | `sample_size_nbinom()` | Null rate ratio (default 1; set < 1 for non-inferiority) |

## S3 class hierarchy

- `sample_size_nbinom_result`: output of `sample_size_nbinom()`; has `print()` and `summary()`.
- `gsNB`: output of `gsNBCalendar()`; inherits from `gsDesign`; has `summary()`, `toInteger()`.
- `nb_sim_data` / `nb_sim_seasonal`: output of `nb_sim()` / `nb_sim_seasonal()`; dispatches `cut_data_by_date()`.
- `mutze_test`: output of `mutze_test()`; has `print()`.

## Vignettes

| Vignette | Topic |
|:---------|:------|
| `sample-size-nbinom` | Fixed-design NB sample size planning |
| `simulation-example` | Basic recurrent-event data generation |
| `seasonal-simulation` | Seasonal rate simulation |
| `group-sequential-simulation` | Calendar-time GS design workflow |
| `ssr-example` | Short SSR workflow |
| `ssr-simulation-study` | Large SSR simulation case study |
| `completers-interim-example` | Completer-based interim looks |
| `non-inferiority-example` | Non-inferiority design |
| `blinded-info-diagnostics` | Edge cases for blinded information |
| `verification-simulation` | Design-formula verification against simulation |

## Common pitfalls

- Do not re-implement SSR logic if `sim_ssr_nbinom()` covers the workflow.
- `mutze_test()` may fall back to Poisson based on `poisson_threshold`; this is defensive, not a model choice.
- `event_gap` changes both counted events and at-risk exposure; it is not cosmetic.
- `check_gs_bound()` accepts `info_col` to select the information column; use that rather than renaming.
- For pkgdown preview, use `gsDesignNB::preview_pkgdown_site()`, not `pkgdown::preview_site()`.

## Disease-area calibration

Published NB trial parameters for reference:

| Therapeutic area | ~λ (events/month) | ~k | Gap |
|:-----------------|:------------------:|:--:|:----|
| Heart failure hospitalizations | 0.3 | 0.5 | Extended hospitalization-free window |
| COPD exacerbations | 0.5–0.8 | 0.5–1.0 | ≥ 28 symptom-free days |
| MS relapses | 1–2/year | moderate–high | Variable |
| Asthma exacerbations | Variable | ~0.8 | Protocol-specific |

## Dependencies

- **Runtime:** gsDesign (≥ 3.9.0), data.table, simtrial, MASS, stats
- **Suggested:** ggplot2, dplyr, gt, DT, shiny, foreach, doFuture, future, testthat, knitr
