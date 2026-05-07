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

> **Last updated:** 2026-04-30 — matches package version 0.3.2

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

### Wald vs score test

`mutze_test()` supports `test_type = "wald"` (default) and `test_type = "score"`:

- **Wald test:** Fits arm-specific NB GLMs; $z = \hat\theta / \text{SE}(\hat\theta)$.
- **Score test:** Fits a pooled null NB GLM via `MASS::glm.nb()`; score statistic
  $U = \sum (y_i - \hat\mu_i) / (1 + \hat k_0 \hat\mu_i)$; Fisher information
  $I_0 = W_1 W_2 / (W_1 + W_2)$; $z = U / \sqrt{I_0}$.

The score test is slightly more conservative under $H_0$ and better controls Type I error
in small-sample / low-rate scenarios.

### Two-variance sample size formula

`sample_size_nbinom()` supports `test_type = "wald"` (default) and `test_type = "score"`:

- **Wald (single variance):** $n_1 = (z_\alpha + z_\beta)^2 V_1 / \theta^2$
- **Score (two variance):** $n_1 = (z_\alpha \sqrt{V_0} + z_\beta \sqrt{V_1})^2 / \theta^2$

where $V_1 = (1/\mu_1 + k_1)(1 + 1/r) + (1/\mu_2 + k_2)/r$ is the alternative-hypothesis
variance, and $V_0 = (1/\mu_0 + k_0)(1 + 1/r)$ uses a pooled rate
$\lambda_0 = (\lambda_1 + r \lambda_2) / (1 + r)$ under $H_0$.

The score formula is matched to the score test's null-variance reference distribution, but it is
not automatically more powerful. In the package superiority grid, Wald/Zhu--Lakkis sizing paired
with the score test preserved Type I error and gave slightly higher score-test power than score
sizing because it retained a small sample-size margin. Compare both sizing rules and verify Type I
error and power by simulation for the design setting.

## Natural-language sample size: extraction and prompting

When the user provides a protocol excerpt, publication quote, or plain-language description of a
trial design, follow this checklist to extract parameters, identify gaps, and compute the sample
size without requiring the user to write code.

### Step 1: Extract what is stated

Parse the text for these quantities. Record each as "found" or "missing":

| Parameter | What to look for | `sample_size_nbinom()` argument |
|:----------|:-----------------|:-------------------------------|
| Control event rate | "annual rate", "event rate", "exacerbation rate", etc. | `lambda1` |
| Treatment effect | Rate ratio, hazard ratio, % reduction, or statement like "no difference" | `lambda2` (= `lambda1 × RR`) |
| Dispersion | "overdispersion", "NB shape", "index parameter k", "theta" | `dispersion` |
| Power | "80% power", "90% power", "> 95% power" | `power` |
| Alpha | "significance level", "one-sided 0.025", "two-sided 0.05" | `alpha` |
| Null hypothesis | "non-inferiority margin", "rule out X-fold increase", "superiority" | `rr0` (default 1.0 for superiority) |
| Follow-up | "52-week treatment period", "12 months follow-up" | `max_followup` |
| Dropout | "30% dropout", "15% discontinuation rate" | See dropout guidance below |
| Event gap | "28 symptom-free days", "new event window" | `event_gap` |
| Enrollment | Accrual rate, enrollment period, ramp-up | `accrual_rate`, `accrual_duration` |
| Trial duration | Total study length | `trial_duration` |

### Step 2: Prompt for essentials

If any of the following are missing, **ask the user** before computing. These cannot be
defaulted safely:

1. **Control event rate and time unit** — "What is the assumed control-arm event rate, and in what
   time unit (events per month, per year)?"
2. **Treatment effect** — "What rate ratio or % reduction is assumed? Or is this a non-inferiority
   design with no assumed difference?"
3. **Dispersion** — "What NB dispersion parameter is assumed? (If unknown, common starting values
   are 0.3–1.0 depending on disease area.)"
4. **Power and alpha** — Usually stated; confirm if ambiguous (e.g., one-sided vs two-sided alpha).
5. **Follow-up duration** — "How long is each patient followed (maximum)?"

### Step 3: Apply sensible defaults for non-essentials

These can be defaulted with a note to the user:

- **Event gap:** If not mentioned, assume `event_gap = 0`. State: "No event gap specified;
  assuming 0. If the protocol requires a minimum inter-event interval, provide the gap in days."
- **Accrual rate and trial duration:** For a fixed-design sample size where every patient gets
  the same maximum follow-up (subject to dropout), the accrual pattern does not affect the
  per-subject information or total N. Supply placeholder values (e.g., `accrual_rate = 100`,
  `accrual_duration = 1`, `trial_duration = max_followup + 1`). State: "Accrual pattern does
  not affect the fixed-design N when all patients have the same max follow-up. Placeholder
  values used. These become important for calendar-time group sequential designs."
- **Test type:** Default to `test_type = "wald"` for replication of published designs (most
  literature uses the Wald/Zhu--Lakkis formula). Offer to compare with score sizing.

### Step 4: Handle dropout correctly

Protocol language about dropout is often ambiguous. Apply these rules:

1. **"X% dropout" or "X% discontinuation"** usually means cumulative incidence over the study
   duration, **not** a rate per time unit. The most common approach in published protocols is to
   compute N assuming full follow-up, then inflate by $1/(1 - X/100)$. Replicate this by setting
   `dropout_rate = 0` and multiplying the result: `N_total <- ceiling(n_evaluable / (1 - X/100))`.

2. **"X% annual dropout rate"** is still typically cumulative incidence (what fraction drop out in
   one year), not a hazard. Convert to instantaneous rate: `dropout_rate = -log(1 - X/100) / T`
   where $T$ is the time period in your chosen unit.

3. **"Dropout + protocol violators"** (as in FLAME) inflates the sample but does not reduce
   exposure in the same way as true dropout. Use the inflation-factor approach.

4. When in doubt, **ask:** "Does the X% dropout figure represent cumulative incidence over the
   study, or a per-unit-time hazard rate? Should it be applied as a sample-size inflation factor
   or modeled through the exposure integral?"

### Step 5: Watch for k parameterization

The NB dispersion parameter has multiple conventions:

- **gsDesignNB convention:** $k$ is the dispersion such that $\text{Var}(Y) = \mu + k\mu^2$.
  Larger $k$ means more overdispersion. `MASS::glm.nb()` returns `theta = 1/k`.
- **Some protocols** use $k$ or $\theta$ to mean the Gamma shape parameter $1/k_{\text{ours}}$.
  If the stated value seems very large (e.g., $k = 50$) or very small for the disease area,
  ask: "Is the dispersion parameter $k$ the overdispersion ($\text{Var} = \mu + k\mu^2$) or
  the Gamma shape parameter ($1/k$)?"

### Step 6: Compute and report

After extraction and defaults, call `sample_size_nbinom()` and present:
- The evaluable sample size (from the formula)
- The inflated sample size (if dropout inflation is used)
- A comparison table if the user wants Wald vs score
- A note about any defaulted values

If the user provided a published target N, show the comparison and note that exact replication
may not be expected due to differences in the specific formula, parameterization, or variance
correction used by the original authors.

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
| `sample_size_nbinom()` | Fixed-design sample size/power; `test_type = "wald"` (V₁ only) or `"score"` (V₀+V₁) |
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
| `mutze_test()` | Wald or score test for treatment effect (`test_type = "wald"` / `"score"`) |
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
- For SSR studies, treat the final test statistic as central. The current production SSR study uses
  the score test for power simulations and compares Wald vs score Type I error at the same nominal
  one-sided $\alpha = 0.025$.
- In that non-binding Type I comparison, Wald is mildly anti-conservative (about 2.8%--3.1%) and
  score is conservative (about 2.1%--2.3%). Prefer the score final test when Type I calibration is
  the priority, then check whether Wald sizing, a higher target power, or a modest information
  margin is needed for power.
- SSR adaptation is driven mainly by interim nuisance estimates and the adaptation cap; the score
  test reduces anti-conservative rejection after adaptation rather than eliminating the need to
  simulate adapted sample-size behavior.
- Blinded SSR preserves masking. Unblinded SSR may be more sample-efficient, but the relative
  behavior is design-specific and should be compared under the trial's nuisance range, adaptation
  cap, and final test statistic.
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
| `test_type` | `mutze_test()`, `sample_size_nbinom()`, `sim_gs_nbinom()`, `sim_ssr_nbinom()` | `"wald"` (default) or `"score"` |

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
| `verification-simulation` | Simulation verification of sample size and test |
| `score-vs-wald-simulation` | 2×2 comparison: Wald/score sizing × Wald/score test (interactive DT/plotly) |
| `ai-skills` | Demonstrates how to use `SKILL.md` and `llms.txt` to guide package-native workflows |
| `blinded-info-diagnostics` | Blinded information estimation diagnostics |
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
