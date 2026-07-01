# Response to Co-author Comments on gsDesignNB Paper

This document tracks how each comment from Honghong Zhou (HZ) and Andrea Maes (AM) was addressed in the revised manuscript (`gsDesignNB-paper.qmd`).

---

## Honghong Zhou Comments

### Introduction / Background

| # | Comment | Response | Location |
|:--|:--------|:---------|:---------|
| HZ-1 | "Dispersion is not always reported. How to obtain it from the medical paper?" | Added paragraph in Introduction explaining that $k$ can be back-calculated from published mean/variance or from the NB shape parameter $\theta_{\text{NB}} = 1/k$. Referenced `estimate_nb_mom()` for method-of-moments estimation from summary statistics. | Introduction, para 2 (after dispersion range) |
| HZ-2 | "Do we need to fact check these?" (literature support for dispersion values) | Expanded and nuanced the disease-area calibration paragraph. Added note that lower dispersion ($k < 0.4$) tends to co-occur with larger gaps in COPD, and higher values ($k \geq 2$) are observed with severe COPD/asthma. Added note about IMPACT using start-to-start counting with always-at-risk offset. | Introduction para 2; Broader parameter sweep "Rationale" paragraph |
| HZ-3 | "Need ref for this claim?" (reference needed) | Citations for @friede2010blinded and @mutze2019group already present; added them explicitly to the Introduction bullet point 3 where monitoring and SSR are first mentioned. | Introduction, bullet 3 |
| HZ-4 | "Need example" | The worked design example (Section 2.2) already provides a concrete illustration. No additional example added since the existing one covers all features. | Section 2.2 |
| HZ-5 | "Should write 'with interim analysis(es)'" | Changed "beyond a fixed final analysis" to "beyond a fixed final analysis to include interim analyses." | Introduction, para 3 |
| HZ-6 | "Need ref for monitoring and SSR" | Added @friede2010blinded and @mutze2019group citations to the Introduction bullet point about SSR. | Introduction, bullet 3 |

### Package overview

| # | Comment | Response | Location |
|:--|:--------|:---------|:---------|
| HZ-7 | "This function for fixed design? Also need horizontal dividers in this table" | Restructured the package overview table with category headers: **Planning**, **Data generation and interim construction**, **Inference and monitoring**, **Adaptive designs**. Renamed "Fixed design planning" to "Fixed design" for clarity. | Section 2.1 table |

### Gap terminology / Missing data

| # | Comment | Response | Location |
|:--|:--------|:---------|:---------|
| HZ-8 | "No 'gap' is used here" | The event gap section (2.3) now explicitly defines what the gap represents, with concrete COPD examples. The terminology is now introduced with a clear definition before being used elsewhere. | Section 2.3 |
| HZ-9 | "Do we have to have this? ZL paper used expected/average FU" | Retained the event gap feature but clarified that it is optional (`event_gap = 0` for protocols without gaps). Added AAR guidance for protocols that don't need gap machinery. | Section 2.3 |
| HZ-10 | "What's this?" (clarification) | Added explicit definitions throughout. See subscript 1/2 definitions and $k_g$ explanation after the Wald formula. | Section 3.2 |
| HZ-11 | "Dropout/MI section is poorly written, contents outdated or untrue. Package does not offer MI." | Rewrote entire section. Renamed to "Dropout and missing data". Removed all MI discussion. Three tight paragraphs: (1) dropout = exposure truncation, all subjects contribute including those with 0 events; (2) MAR validity including treatment-policy/retrieved-dropout estimand; (3) MNAR as future extension with concrete example (experimental dropouts revert to control rates). | Section 2.6 |

### Sample size methodology

| # | Comment | Response | Location |
|:--|:--------|:---------|:---------|
| HZ-12 | "Define footnote 1 and 2" (undefined subscript notation) | Added: "Throughout, subscript $g = 1$ denotes the control arm and $g = 2$ the experimental arm; $k_g$ is the arm-specific dispersion, which reduces to a common $k$ when specified as a scalar." | Section 3.2, after Wald formula |

### Calendar-time GS section

| # | Comment | Response | Location |
|:--|:--------|:---------|:---------|
| HZ-13 | "Inconsistent notation: $i$ was used for subject, but here for looks" | The formula already uses $j$ for looks ($C_j$, $\delta_j$, $\mathcal{I}_j$, $\tau_j$) and $i$ for subjects. Verified notation is consistent --- the concern may have been about a draft that was subsequently corrected. No change needed. | Section 4.1 |

### Reproducible workflows

| # | Comment | Response | Location |
|:--|:--------|:---------|:---------|
| HZ-14 | "Do we need to show the results in the paper?" | Retained current approach: code shown in paper, results referenced to package vignettes. This follows the JSS software-paper convention where supplementary articles serve as executable material. | Sections 5.1--5.2 |

---

## Andrea Maes Comments

### Abstract

| # | Comment | Response | Location |
|:--|:--------|:---------|:---------|
| AM-1 | "Can we say constant rate if we are concerned about seasonality?" | Added clause to abstract: "a seasonal-rate extension relaxes the constant-rate assumption for respiratory and other calendar-time-dependent applications." | Abstract |
| AM-2 | "Is event gap interval-censoring or protocol-defined counting rule?" | Clarified in abstract: replaced "protocol-defined gaps between counted events" with "protocol-defined minimum intervals between countable events (combining event duration and required symptom-free days)." | Abstract |

### Introduction

| # | Comment | Response | Location |
|:--|:--------|:---------|:---------|
| AM-3 | "1.5 is moderate dispersion. For severe COPD and asthma, much higher values." | Expanded range from "0.3 to 1.5" to "0.3 to 3.0 or higher". Added calibration scale: $k < 0.3$ near-Poisson, $0.3$–$1.0$ moderate, $> 1.5$ high, $> 5$ extreme. | Introduction, para 2 |
| AM-4 | "Can we provide guidelines for what might be extreme overdispersion?" | Added to Introduction: $k > 5$ is extreme (where MoM fallback becomes relevant). The MoM fallback threshold in the code is $\hat{k} > 20$. | Introduction, para 2 |
| AM-5 | "Reconcile constant rate with seasonality" | Addressed via abstract clause about `nb_sim_seasonal()` and in the dropout section noting that `nb_sim_seasonal()` captures calendar-time variation while retaining the observed-exposure likelihood structure. | Abstract; Section 2.6 |
| AM-6 | "Clarify what event gap means" | Rewrote Section 2.3 with explicit definition: event gap = total exclusion window = event duration + protocol-defined separation period. Concrete COPD example (7-day steroid + 7-day separation = 14-day gap). | Section 2.3 |
| AM-7 | "Do we just assume the user will correctly adjust follow-up time?" | Added explicit guidance in Section 2.3: user is responsible for summing event-duration and separation components. Added AAR vs ERT distinction with guidance on when to set `event_gap = 0`. | Section 2.3 |

### Worked example / Event gap specification

| # | Comment | Response | Location |
|:--|:--------|:---------|:---------|
| AM-8 | "Can we discuss minimum gap without getting into event duration? Does tool simulate event duration?" | Section 2.3 now explicitly addresses this. The tool does not simulate event duration separately; `event_gap` is the total exclusion window and the user combines event duration and separation period. Simulation enforces the gap during event generation (not post-hoc). | Section 2.3 |
| AM-9 | "Should we be clear that $k = 1/\text{shape}$?" | Added to Section 3.1: "$\theta_{\text{NB}} = 1/k$; thus $k = 1/\theta_{\text{NB}}$." Also added notation convention paragraph. | Section 3.1 |
| AM-10 | "Assumes 'Always at Risk'?" | Clarified that the package implements ERT (excluding recovery time) by default. AAR users should set `event_gap = 0`. | Section 2.3 |

### Calendar timing / Event gaps

| # | Comment | Response | Location |
|:--|:--------|:---------|:---------|
| AM-11 | "Event gap removes event duration as well (e.g., steroid course)" | Addressed in Section 2.3: gap = event duration + separation period. User sums components. | Section 2.3 |
| AM-12 | "Gap definition under different estimands (ERT vs AAR)" | Added explicit ERT vs AAR distinction. ERT: gap affects both counting and exposure offset (package default). AAR: set `event_gap = 0`, combine events externally. | Section 2.3 |

### Missing data section

| # | Comment | Response | Location |
|:--|:--------|:---------|:---------|
| AM-13 | "Including 0 events" | Added: "including subjects who discontinue before experiencing any event (contributing observed exposure $> 0$ and event count $= 0$)." | Section 2.6, para 1 |
| AM-14 | "MAR assumption challenged by seasonality" | Added note that `nb_sim_seasonal()` captures within-subject rate variation while the likelihood retains observed-exposure structure. | Section 2.6, para 2 |
| AM-15 | "Treatment policy scenario with retrieved dropout" | Added: "This includes the common treatment-policy estimand scenario in which subjects who discontinue treatment are followed for the planned study duration ('retrieved dropout')." | Section 2.6, para 2 |
| AM-16 | "Active would be expected to wash out and behave like placebo after discontinuation" | Addressed in MNAR paragraph: "imputing post-dropout events under reference-based assumptions (e.g., experimental-arm dropouts revert to control-arm rates)." | Section 2.6, para 3 |
| AM-17 | "Observed exposure time?" | Changed wording from "observed time at risk" to "observed exposure time" consistently. | Section 2.6 |
| AM-18 | "Versus planned" | Clarified that the analysis uses observed (not planned) exposure throughout. | Section 2.6 |

### Score test discussion

| # | Comment | Response | Location |
|:--|:--------|:---------|:---------|
| AM-19 | "Sentence hard to read. 'Modest' is vague." | Rewrote recommendation as a bullet list of specific conditions. Replaced "sample size is modest" with "total sample size is small (e.g., fewer than $\sim$200 per arm)". | Section 3.2 |

### Jensen correction / Event gaps

| # | Comment | Response | Location |
|:--|:--------|:---------|:---------|
| AM-20 | "Hard to think about without event duration. Important for bronchiectasis." | Revised opening of Section 3.3 to define $\gamma$ as "total exclusion window that may include event duration plus a protocol-defined separation period." | Section 3.3, opening sentence |
| AM-21 | "Rate interpretation: ERT = rate per at-risk time vs calendar time" | Added: "The rate $\lambda$ is therefore interpreted as the rate per unit of *at-risk* time (excluding recovery), not per unit of calendar time." | Section 3.3, implementation paragraph |
| AM-22 | "What if you just reduced the time at risk offset? How does that compare?" | Added: "An alternative approach --- simply reducing the exposure offset by the gap periods without modelling the gap in the event-generation process --- would understate the effective rate and produce a different (generally less accurate) power calculation." | Section 3.3, implementation paragraph |
| AM-23 | "Events should be simulated, not just discarded inside the gap" | Clarified that `nb_sim()` generates events sequentially with the gap enforced during generation: "after each counted event, the next event cannot occur until at least $\gamma$ time units later (this is enforced during event generation, not applied as a post-hoc filter)." | Section 3.3, implementation paragraph |

### Parameter sweep

| # | Comment | Response | Location |
|:--|:--------|:---------|:---------|
| AM-24 | "Dispersion values too low (can get into 2s and 3s)" | Added sweep table caption note: sweep extends to $k = 1.0$; for $k \geq 2$ the correction is even larger but the second-order approximation loses accuracy and simulation is the primary validation tool. | Table caption for sweep grid |
| AM-25 | "In exacerbation trials, lower dispersion with larger gaps. 28 days is largest gap." | Added to disease-area calibration: "lower dispersion ($k < 0.4$) tends to co-occur with larger gap definitions." | Broader sweep "Rationale" paragraph |
| AM-26 | "IMPACT assumed start-to-start and always-at-risk (duration = 0)" | Added: "Some published analyses (e.g., IMPACT) effectively assume zero event duration by using start-date-to-start-date counting and an always-at-risk exposure offset; the `event_gap` parameter can accommodate either convention, as described in Section 2.3." | Broader sweep "Rationale" paragraph |

### Example code

| # | Comment | Response | Location |
|:--|:--------|:---------|:---------|
| AM-27 | "Inputs hard to follow" | Added inline comments to all data frame definitions in Example 1 code. | Section 5.1 |
| AM-28 | "Is this the accrual rate? Or event rate?" | Added comment: `# enrollment: 12 subjects/month for 6 months`. | Section 5.1 |
| AM-29 | "Not clear what fail_rate means" | Added comment: `# event rates (fail_rate follows simtrial convention)`. Added per-column comments for `rate` and `dispersion`. | Section 5.1 |

### Discussion

| # | Comment | Response | Location |
|:--|:--------|:---------|:---------|
| AM-30 | "More detail requested" | The Discussion references the vignettes for fuller treatment. Updated the dropout/MAR paragraph in the Discussion to include treatment-policy and MNAR reference-based language. Detailed methodology is in the vignettes. | Discussion |

---

## Summary of Changes by Section

| Section | Changes |
|:--------|:--------|
| **Abstract** | Added seasonal-rate extension note; clarified event gap definition |
| **Introduction** | Expanded dispersion range to 0.3–3.0+; added calibration scale; added guidance for obtaining $k$ from literature; added `estimate_nb_mom()` reference; fixed "interim analyses" wording; added @friede2010blinded/@mutze2019group citations to SSR bullet |
| **Section 2.1 (Scope)** | Restructured overview table with category headers |
| **Section 2.3 (Calendar timing)** | Major expansion: event duration + gap framing; ERT vs AAR distinction; start-to-start vs end-to-start convention; simulation mechanism clarification; concrete COPD example |
| **Section 2.6 (Dropout)** | Renamed from "Dropout, MAR, and imputation scope"; complete rewrite removing MI discussion; added 0-event contributor note; treatment-policy/retrieved-dropout framing; seasonality note; reference-based MNAR future direction |
| **Section 3.1 (NB model)** | Added $k = 1/\theta_{\text{NB}}$ explicit inversion; added notation convention for $k$ vs $K$ |
| **Section 3.2 (Wald/score)** | Added subscript 1/2 and $k_g$ definitions; rewrote score-test recommendation as bullet list with specific thresholds |
| **Section 3.3 (Jensen correction)** | Redefined $\gamma$ as total exclusion window; clarified rate as per-at-risk-time; addressed offset-reduction alternative; clarified simulation is generative not post-hoc; expanded sweep table caption for high-$k$ limitation |
| **Section 3.3 (Broader sweep)** | Expanded disease-area calibration; added COPD gap/dispersion co-occurrence note; added IMPACT start-to-start note |
| **Section 5.1 (Example 1)** | Added inline comments to all data frame definitions |
| **Discussion** | Updated dropout paragraph with treatment-policy and MNAR reference-based language |

---

## Items Not Changed (with Rationale)

1. **Notation $i$ vs $j$ for looks (HZ-13)**: The formula already uses $j$ for looks and $i$ for subjects. No inconsistency found in the current text.

2. **Show results in paper (HZ-14)**: Retained JSS convention of code-in-paper + results-in-vignettes. The manuscript would become unwieldy with full output tables.

3. **Rename `fail_rate` variable (AM-29)**: Added comments instead. The name follows the `simtrial` package convention used throughout the gsDesign ecosystem; renaming would break consistency.

4. **Expand simulation sweep beyond $k = 1.0$ (AM-24)**: Would require new computation. Added caption note acknowledging the limitation and pointing to simulation for $k \geq 2$.

5. **Simulate event duration separately (AM-8)**: The package models the total exclusion window, not separate event duration and gap components. This is a deliberate simplification: the user combines components into `event_gap`. A future extension could model event duration explicitly.
