# Sample Size for Group Sequential Survival Trials with Calendar Gaps Between Interim Analyses and Piecewise Event Dynamics

**Draft manuscript** (target length: 10–15 journal pages; ~6,000 words)

**Authors:** Keaven M. Anderson, Hongtao Zhang, Andrea Maes  

**Affiliation:** Merck & Co., Inc., Rahway, NJ, USA

**Keywords:** group sequential design; survival; sample size; interim analysis; alpha spending; sample size re-estimation; software

---

## Abstract

Interim analysis schedules for time-to-event trials are often specified using a mixture of calendar-time rules, minimum follow-up, and event-driven criteria. Naive sample size calculations that ignore **mandatory gaps between successive interim looks**—for example, minimum time from one analysis to the next while events continue to accrue—can materially misstate expected information, power, and spending paths. We also treat **patient-level minimum gaps between successive events** (e.g. at least 20 days between one qualifying event and the next on the same patient): in simulation, any event that would fall inside that window is **discarded** and does **not** count toward protocol-defined event totals. We describe a **unified expected-events framework** under piecewise exponential failure, dropout, and enrollment (Lachin–Foulkes style), show how **calendar floors** and **minimum spacing between analyses** enter the timing solution, and link these to **error-spending** group sequential designs, including **gapped spending functions** that formalize intervals of information time in which no bound is tested. We summarize **simulation-based verification** used for both fixed sample size and **sample size re-estimation** (conditional-power adaptation), including scenarios with **inter-event censoring** on the patient timeline. Implementation references include the **gsDesign** R package (survival / group sequential core) and complementary **gsDesignNB** tools for recurrent and overdispersed count outcomes with **event gaps**. Finally, we discuss how **AI-assisted development** accelerated implementation, with verification emphasizing large-scale simulation and regression testing—consistent with a shift toward **delegation and verification** of agent-generated code [Dohmke & Kalliamvakou, 2025](https://ashtom.github.io/developers-reinvented). Compared with many commercial or narrow-purpose tools, open-source stacks can provide extensive documentation, vignettes, and reproducible examples.

---

## 1. Introduction

Group sequential methods for time-to-event endpoints are standard in confirmatory oncology and other therapeutic areas [Jennison & Turnbull, 2000]. Software typically combines (i) a **fixed-design** sample size or event count for a proportional hazards alternative, (ii) a **spending function** to allocate Type I/II error across information fractions, and (iii) boundaries obtained by numerical inversion. Regulatory practice often adds constraints that are **not** captured by a pure “fraction of final events” rule: analyses may be tied to **calendar dates**, require **minimum patient exposure**, or impose **minimum elapsed time** between successive database locks while events continue to accumulate.

We use the phrase **“gaps between events”** in three complementary senses that arise in applications:

1. **Operational / calendar gaps between interim analyses:** successive looks must occur no earlier than \(T_{i-1}+\Delta_i\) after the previous analysis, possibly concurrent with additional accrual and ongoing event generation.
2. **Structural gaps in the intensity model:** piecewise-constant hazards with change-points (e.g., distinct post-randomization periods) imply **piecewise** accumulation of expected events; “gaps” in spending can also be imposed on the **information scale** via **gapped spending functions** that allocate no error spending over specified information-time intervals.
3. **Patient-level minimum gap between successive events:** a protocol may require that **only** events separated by at least \(\delta\) **patient time** (e.g. \(\delta = 20\) days) count toward the analysis event total. **Simulated** events that fall inside \((t_{\text{prev}},\, t_{\text{prev}}+\delta]\) after a previous qualifying event at \(t_{\text{prev}}\) for the same patient are **excluded** from the count—they do not contribute to information. This rule mimics adjudication windows, avoidance of double-counting clinically indistinguishable events, or operational definitions of “new” progression.

The methodological contribution of this paper is not a single new asymptotic result but a **clear definition and derivation** of the **expected event** and **timing** calculations that underpin sample size when such gaps are specified, and their integration with **alpha/beta spending** and **SSR** (sample size re-estimation) tools. Standard **closed-form** expected-event formulas (Section 2) typically characterize **first** or **total** events under a PH model; **patient-level inter-event gaps** are easiest to incorporate rigorously via **simulation** (Section 5), where the counting rule is applied explicitly. Reference implementations for the piecewise calendar and group sequential machinery appear in the **gsDesign** package [Anderson and contributors], including `nSurv()`, `gsSurv()`, `gsSurvCalendar()`, `gsSurvPower()`, `sfGapped()`, and `ssrCP()`. The **gsDesignNB** package builds on **gsDesign** for negative binomial / recurrent-event settings with **event gaps** in sample size and simulation [Zhu & Lakkis, 2014; package documentation].

**Outline.** Section 2 sets notation and derives expected events under a piecewise model, and defines **per-patient inter-event spacing**. Section 3 connects calendar and event-driven rules to interim timing. Section 4 reviews group sequential integration and gapped spending. Section 5 outlines simulation studies: **§5.1** validates analytical expected events against Monte Carlo under the documented model; further subsections cover operating characteristics, SSR, and regression testing. Section 6 describes AI’s role and the **verification strategy**. Section 7 concludes.

---

## 2. Piecewise exponential model and expected events

### 2.1 Notation

Consider a two-arm trial with **control** failure hazard \(\lambda_C(t)\) modeled as **piecewise constant** on intervals defined by cumulative cutpoints \(\mathbf{S}=(S_1,\ldots)\) (possibly empty for a single exponential piece). **Dropout** hazards \(\eta_C(t)\) (and experimental-arm analogs) may share the same or analogous structure. **Enrollment** follows a piecewise constant rate \(\gamma_j\) over deterministic periods of length \(R_j\) with \(\sum_j R_j\) possibly exceeding planned total trial duration; enrollment may be truncated by calendar time \(T\).

Let **HR** denote the hazard ratio under the alternative (experimental vs control), and **HR\(_0\)** under the null hypothesis for testing. The experimental arm hazards are obtained by multiplying control hazards by HR (or HR\(_0\)) as appropriate for the design.

### 2.2 Expected events in one arm (derivation sketch)

For a single stratum and treatment arm, conditional on the piecewise structure, the **expected number of events** by calendar time \(T\) follows from integrating the product of **risk sets** and **event intensity** over sub-intervals where \((\lambda(\cdot),\eta(\cdot),\gamma(\cdot))\) are constant. Partition \([0,T]\) at all change-points from failure/dropout **and** enrollment. On each sub-interval of length \(s_\ell\), the contribution to the expected event count takes the schematic form

\[
\mathbb{E}[\text{events in interval } \ell] \;=\; \int_{\text{interval }\ell} \text{(risk} \times \lambda)(u)\,du
\]

which, under piecewise-constant assumptions, reduces to **closed-form contributions** proportional to terms such as

\[
\frac{\lambda}{\lambda+\eta}\,\Big(1 - e^{-(\lambda+\eta)s}\Big)
\]

weighted by **survival of enrollment and failure** carried forward across intervals (see `eEvents1()` in gsDesign: the implementation follows the algebraic decomposition in the package technical documentation, aggregating period-by-period contributions).

**Algorithmic construction.** Let \(0=t_0<t_1<\cdots<t_M=T\) denote the ordered distinct cutpoints from \(\mathbf{S}\), enrollment boundaries derived from \(\mathbf{R}\), and \(T\) itself. For each half-open interval \((t_{m-1},t_m]\), constants \(\lambda_m,\eta_m\) apply to failure and dropout; the **enrollment rate** is piecewise on the **time-since-study-start** axis, which induces a **backward** construction of the mass of patients still at risk for follow-up at \(t_{m-1}\). The expected proportion failing in the interval is the product of (i) the **conditional** probability of failure before dropout given entry before \(t_{m-1}\), and (ii) the **expected** number still enrolled and not yet failed or dropped out at \(t_{m-1}\). Summing over \(m\) yields \(\mathbb{E}[N_{\text{events}}(T)]\) for that arm. This is the **same structural decomposition** used by Lachin and Foulkes [1986], extended to **vector** \(\lambda\) and \(\eta\) across periods.

The **total expected events** for the trial is obtained by summing over arms and strata (independent strata), and scaling by randomization ratio.

### 2.3 Fixed \(L\)-sample size (events) under PH alternative

Given Type I error \(\alpha\), power \(1-\beta\), and one-sided or two-sided testing conventions, the **fixed-design** event requirement \(d\) solves the standard Schoenfeld-style relationship using the **log-rank** asymptotics, with non-centrality driven by \(\log(\text{HR})\) and the expected event count. The `nSurv()` routine solves for enrollment scale and event targets consistent with the **piecewise** enrollment and failure model—not only the scalar exponential case.

### 2.4 Patient-level minimum gap between successive events (inter-event window)

Some protocols define the **analysis event count** using a **minimum separation** between successive **qualifying** events on the **same patient**. Let \(\delta > 0\) denote that minimum (e.g. **20 days** on the study or calendar timeline from randomization, or from the previous event—definitions should match the SAP). Suppose patient \(j\) would generate a **sequence** of potential event times \(T_{j,1} < T_{j,2} < \cdots\) under a simulation model (e.g. from a semi-Markov or recurrent-event process, or from multiple draws in a composite endpoint).

**Counting rule (simulation).** Sort potential event times for patient \(j\) as \(T_{j,1} < T_{j,2} < \cdots\). Let \(C_{j,1}=T_{j,1}\) be the first **counted** time. Recursively: the next **counted** time is the **first** \(T_{j,r}\) such that \(T_{j,r} \ge C_{j,s} + \delta\), where \(C_{j,s}\) is the previous counted time; all **earlier** \(T_{j,r}\) in \((C_{j,s},\, C_{j,s}+\delta)\) are **discarded**—they **do not count** toward the trial’s cumulative event total. Equivalently: after a counted event at calendar time \(t\), **no** simulated event strictly before \(t+\delta\) can become the next counted event.

This rule **thins** recurrent or densely simulated events relative to a model with no minimum gap. The **operational** information accrual is therefore **slower** than the naive expected event total from Section 2.2 would suggest whenever multiple events per patient are possible under the generative model. **Analytical** adjustment requires a recurrent-events model with **gap time** constraints (e.g. gap-time Cox or frailty formulations); a practical approach for design is to **calibrate** \(d\) or enrollment inflation by **Monte Carlo** under the \(\delta\)-rule until targeted power is met.

**Example.** If \(\delta = 20\) days and the simulator proposes a second event 12 days after the first, that second event **does not count**. Only events occurring **at least 20 days** after the previously **counted** event on that patient can enter the information denominator.

---

## 3. Calendar constraints and “gaps” between successive interim analyses

### 3.1 Information versus calendar time

In event-driven designs, the **information fraction** at analysis \(i\) is often proxied by cumulative events divided by the planned final events. However, when interim analyses are **also** constrained by calendar rules, the realized **sequence \((T_1,\ldots,T_K)\)** of analysis times is the solution of **feasibility inequalities**:

- **Planned calendar time** (if any): \(T_i \ge C_i\).
- **Minimum gap** after the previous analysis: \(T_i \ge T_{i-1}+\Delta_i\) (months or days on study clock).
- **Minimum exposure:** enough patients with \(\ge\) minimal follow-up.
- **Event threshold:** if the \(i\)-th look is triggered by **target events** \(d_i\), then \(T_i\) is the **smallest** time \(t\) such that **expected** (or observed) events \(\ge d_i\), subject to the floors above.

This matches the **priority/floor** logic documented for `gsSurvPower()`: take the **maximum** of applicable **floor times** (calendar, gap-from-previous, enrollment+follow-up), then adjust upward if event targets require a later time, optionally capped by **maxExtension** rules.

### 3.2 Why ignoring gaps biases sample size

If one assumes interim analyses occur at fixed fractions of final events **without** enforcing \(T_i \ge T_{i-1}+\Delta_i\), the **expected information path** is systematically misspecified whenever \(\Delta_i>0\) and enrollment continues between looks: **events accrue during the gap**, shifting the spending time axis relative to a naive “events-only” schedule. The proposed **expected-events + floor** construction aligns the **timing** vector fed to `gsDesign()` with trial operations.

### 3.3 Interaction with patient-level inter-event gaps (Section 2.4)

**Trial-level** gaps \(\Delta_i\) between database locks (Section 3.1) and **patient-level** gaps \(\delta\) (Section 2.4) operate on **different** time scales: the former constrain **when** looks occur; the latter constrain **which** patient-level events increment the event tally. A design that sets target events \(d_i\) assuming every simulated failure counts will be **overpowered** or **under-accrue** information if the SAP applies a **δ-rule** that reduces the realized event rate. Reporting simulations **with and without** \(\delta\)-thinning is therefore essential when recurrent or composite events are modeled.

---

## 4. Group sequential integration and gapped spending

### 4.1 Error spending on information time

Let \(t \in [0,1]\) denote **information fraction** (often event fraction). Standard spending functions \(f_\alpha(t)\) allocate cumulative Type I error; analogs allocate \(\beta\) for futility under the alternative. Given **timing** \(t_i\) at look \(i\), group sequential boundaries are computed by inverting spending constraints jointly with the correlation structure of the sequential log-rank (or mean-difference) statistics.

### 4.2 Gapped spending (`sfGapped`)

Regulatory or operational considerations sometimes require **no formal efficacy test** over an interior range of information times—e.g., **no interim efficacy** between 20% and 60% of information, while still allowing a later interim. The **gapped spending function** composes a baseline spending shape on a sub-interval \((\tau_1,\tau_2)\) and forces **full spending** to occur by the first analysis **after** the gap, as formalized in `sfGapped()` (see package documentation: spending is **flat** across the “hole” on the information scale, then resumes). This is distinct from calendar gaps in Section 3 but **interacts** with them when information fractions are driven by **expected events** computed under the calendar-consistent schedule.

### 4.3 Putting pieces together

A complete **survival group sequential** specification combines:

1. **Piecewise** enrollment, failure, dropout (Section 2).
2. **Timing** \(T_i\) from Section 3 (calendar gaps + events).
3. **Spending** \(f_\alpha, f_\beta\) (and harm-bound spending for three-boundary designs), possibly **gapped**.
4. Inflation of enrollment or trial duration to meet power at the final analysis.

`gsSurv()` and `gsSurvCalendar()` implement complementary paths: event-driven timing versus **calendar-first** construction with implied information fractions.

---

## 5. Simulation studies

We recommend reporting **three layers** of simulation-related evidence in a submission: (0) **consistency** of analytical expected events with Monte Carlo under the **documented** model; (1) broader operating characteristics; (2) SSR.

### 5.1 Accuracy of the sample size method (expected total events)

The **gsDesign** documentation for `nSurv()` / `eEvents()` specifies a **Lachin–Foulkes**-style construction: piecewise constant enrollment, failure, and dropout; proportional hazards between arms; **total expected events** under the alternative given by closed-form period contributions (`eEvents1()` and aggregation in `eEvents()`).

**Evaluation via simulation:** under the **same** assumptions used analytically—uniform enrollment on \([0,\sum R]\), independent exponential failure and dropout (**competing risks**), failure counted only if it occurs by calendar time \(T\)—one can simulate many replicate trials at fixed sample size `round(nSurv()$n)` and compare the **Monte Carlo mean** of the total failure count to **`nSurv()$d`**. Agreement within a small multiple of the **Monte Carlo standard error of the mean** supports correct implementation of the integral formulas for the **scalar-hazard** case.

**Reproducible script:** `papers/simulation_sample_size_accuracy.R` in the **gsDesignNB** repository (run with `Rscript`; optional `SIM_NSIM` for the number of replicates). A short guide is in `papers/README-simulation.md`. The script also exports `thin_interevent_gap()` for the patient-level \(\delta\) rule (Section 2.4) when used with recurrent-event simulators.

**Regression testing:** `tests/testthat/test-independent-test-simulation-accuracy.R` checks that the Monte Carlo mean stays within tolerance of `d` for a representative design.

**Piecewise failure hazards** (`S` not `NULL`): the same analytical `d` remains valid; a piecewise **simulator** from enrollment time must respect interval boundaries (extend the scalar simulator or compare arm/stratum totals to `eEvents()`). The repository script documents this extension point.

### 5.2 Monte Carlo operating characteristics (fixed design / GSD)

- **Data generating process:** piecewise exponential or Weibull deviations to stress PH assumptions; optional non-proportional hazards; **recurrent or multi-event** generators when patient-level \(\delta\) is relevant.
- **Seasonality:** when event intensity is expected to vary over calendar seasons (for example, respiratory or infectious-disease endpoints), **gsDesignNB** supports direct Monte Carlo generation with season-specific rates via `nb_sim_seasonal()`, after which the same cut-and-analysis workflow can be applied.
- **Patient-level \(\delta\)-rule:** after drawing event times per patient, apply the **Section 2.4** filter so events within \(\delta\) of the previous **counted** event do not count; accumulate trial-level events from the thinned streams.
- **Outputs:** empirical Type I error (where feasible), power, distributions of interim analysis times, distribution of cumulative spending paths, early stopping rates, **distribution of time to reach \(d\) counted events** with versus without thinning.
- **Comparison:** designs **with** versus **without** enforced minimum gaps between looks; **with** versus **without** gapped spending when relevant; **\(\delta = 0\)** versus **\(\delta = 20\)** days (or other SAP values).

### 5.3 Sample size re-estimation (SSR)

The `ssrCP()` framework implements **conditional-power**–based SSR for two-stage adaptations following a group sequential plan [Lehmacher & Wassmer, 1999; related references in package help]. Simulations should report:

- unconditional power (`Power.ssrCP()`),
- distributions of adapted stage-2 sample size,
- control of Type I error under the combination-test conventions implemented.

These simulations complement **analytical** expected-event calculations: they validate **behavior under model misspecification** and **discrete** event accrual.

For **negative binomial recurrent-event** designs, **gsDesignNB** extends the same idea: interim analyses are naturally keyed to **statistical information** (variance of the log rate-ratio under a working NB model), not only to calendar milestones. **Sample size re-estimation** is most valuable when planning assumptions for **nuisance parameters**—especially **overdispersion** and accrual or follow-up intensity—are uncertain; updating the projected total enrollment or follow-up after an interim look limits unconditional power loss without rewriting the group-sequential boundary logic, provided the adaptation rule is **pre-specified** and operating characteristics are simulated.

**Blinded versus unblinded nuisance updates** trade scientific and regulatory goals. **Blinded** SSR pools arms to estimate a common dispersion and pooled rate, which preserves masking of the treatment effect at the adaptation; in simulation studies reported in the package, blinded re-estimation can yield **larger** adapted sample sizes than **unblinded** updates that use arm-specific fits, because the blinded working model is intentionally conservative about heterogeneity. Where a data monitoring charter and regulators accept **unblinded** nuisance-parameter updates with appropriate error-spending or combination-test discipline, that path is often **more sample-efficient**; the choice should be justified with **simulated Type I error and power** under the intended SSR rule, not by analytical shortcut alone.

### 5.4 Illustrative simulation factors (reporting template)

| Factor | Levels (example) |
|--------|-------------------|
| Hazard ratio under alternative | 0.65, 0.75, 0.85 |
| Minimum gap \(\Delta\) between looks | 0, 3, 6 (months) |
| Patient-level min gap \(\delta\) between counted events | 0, 20 (days) |
| Enrollment slowdown | nominal \(\gamma\) vs 80% of nominal |
| Gapped spending interval \((\tau_1,\tau_2)\) | none; (0.2, 0.5) on information scale |
| SSR rule | `ssrCP()` with `maxinc` capped |

For each configuration: **10,000+** trial replicates (more for Type I near \(\alpha\)). Summarize **median** and **IQR** of calendar times at looks; **empirical power**; **mean** stage-2 \(n\) under SSR.

### 5.5 Software regression as a complement to Monte Carlo

Beyond standalone simulation scripts, the **gsDesign** and **gsDesignNB** repositories use **automated** tests (e.g., `testthat`, snapshot tests for plots where applicable) to catch **drift** in computations and plotting when refactoring. These are **not** substitutes for methodological simulation but **guardrails** ensuring released code matches stored reference outputs—especially important when AI-assisted edits increase commit velocity.

---

## 6. Role of AI in methods development, coding, and verification

Recent accounts of **AI-experienced developers** emphasize a shift from **line-by-line coding** to **delegation** (framing tasks, curating context, reviewing plans) and **verification** (testing, review, security) [Dohmke & Kalliamvakou, 2025](https://ashtom.github.io/developers-reinvented). That perspective aligns with statistical software for regulated contexts: **correctness** is established by **independent checks**, not by authority of the generator.

In the **gsDesign** / **gsDesignNB** ecosystem, AI tools can accelerate:

- exploratory algebra for **piecewise** interval bookkeeping;
- refactoring for readability;
- documentation and vignette expansion;
- test authoring for edge cases.

**Verification** prioritized:

- extensive **testthat** coverage and visual regression where appropriate;
- **simulation** studies as the primary external validity check for complex integrals and adaptive paths;
- **R CMD check**, vignette rebuilds, and **pkgdown** consistency.

This also suggests a **layered software delivery model**. Core estimands,
design calculations, simulation engines, and summaries belong in the
**package**, where they are versioned, scriptable, and covered by tests. A
**Shiny** interface can sit above the same functions for interactive scenario
exploration (e.g., accrual, dispersion, event gaps, SSR caps), but should not
re-implement the statistical core. **AI skills** or similar IDE helpers can
route users to the right entry points, check unit consistency, and scaffold
workflows; they are useful accelerators for contributors and power users, but
they are **not** substitutes for validated APIs or manuscript-level method
descriptions.

This mirrors the “realistic optimist” stance in [Dohmke & Kalliamvakou, 2025](https://ashtom.github.io/developers-reinvented): adopt AI for throughput while **raising the bar** on validation discipline.

---

## 7. Discussion and conclusions

1. **Modeling gaps:** Enforcing **minimum times between interim analyses** and using **piecewise** hazard/enrollment specifications materially affects **expected events** and **information growth**. The **expected-events + floor** construction gives a coherent account suitable for **gsDesign**-class implementations.

2. **Spending gaps:** **Gapped spending** (`sfGapped`) formalizes **intervals of information time** without testing, complementary to calendar operational constraints.

3. **Simulations:** Operating characteristics should be reported for both **fixed GSD** and **SSR** pipelines; these serve as the principal **external** check on complex code paths.

4. **AI & software:** Open-source stacks (**gsDesign**, **gsDesignNB**) can pair **validated package APIs**, optional **interactive** layers, and **AI-assisted** workflow support with simulation-backed verification.

5. **Limitations:** The piecewise exponential framework is an approximation; deviations (delayed effects, cure fractions) require extended models or more aggressive simulation. Calendar rules in real trials involve **random** accrual jitter and operational delays not fully captured by deterministic \(\Delta_i\). **Patient-level** \(\delta\)-rules are handled straightforwardly in **simulation**; closed-form **expected counted events** under general recurrent processes with hard gaps is **not** included in standard PH fixed-design formulas and may require dedicated recurrent-event theory or empirical calibration.

6. **Mutze-type interim tests and dispersion:** For recurrent **count** outcomes, the **Mutze et al.** Wald test implemented in **gsDesignNB** as `mutze_test()` fits a negative binomial log-linear model with an offset for exposure. When maximum-likelihood estimation implies **essentially Poisson** variation (very large estimated NB shape \(\theta\)) or **numerically unstable** extreme overdispersion (very small \(\theta\)), the implementation **falls back to Poisson** regression for the Wald statistic. The switch is governed by `poisson_threshold` \(T\) (default \(T=1000\)): Poisson is used if \(\theta \notin [1/T, T]\). **Lowering \(T\)** triggers the Poisson fallback **sooner**—i.e., whenever you prefer to treat the data as **without extra-Poisson dispersion** once \(\theta\) is only moderately large—and is a practical way to align the interim model with a **Poisson** SAP without forcing `method = "poisson"` for every analysis. **Raising \(T\)** keeps the NB fit over a wider \(\theta\) range. Sensitivity analyses across \(T\) or explicit `method = "poisson"` versus `"nb"` should accompany any change from package defaults.

7. **Method of moments as a complement:** Group-sequential simulation in **gsDesignNB** (`sim_gs_nbinom()`) records **maximum-likelihood** and **method-of-moments** dispersion estimates and the associated blinded and unblinded **information** paths. If ML-based nuisance estimates behave poorly at small interim sample sizes, comparing MoM and ML paths in simulation supports whether MoM (or blinded helpers that fall back when a blinded NB fit fails) is preferable for **SSR** or information monitoring. This does not replace protocol pre-specification but gives a **data-driven** way to stress-test default switching rules before locking them.

---

## Acknowledgments

The authors thank colleagues at Merck & Co., Inc., Rahway, NJ, USA for support and feedback.

---

## References

Anderson, K. et al. (ongoing). **gsDesign**: Group Sequential Design. R package.  
URL: `https://github.com/keaven/gsDesign` and `https://keaven.github.io/gsDesign/`.

Anderson, K. et al. (ongoing). **gsDesignNB**: Sample Size and Simulation for Negative Binomial Outcomes. R package.  
URL: `https://github.com/keaven/gsDesignNB` and `https://keaven.github.io/gsDesignNB/`.

Bauer, P. & Köhne, K. (1994). Evaluation of experiments with adaptive interim analyses. *Biometrics*.

Chen, Y. H., DeMets, D. L., & Lan, K. K. G. (2004). Increasing the sample size when the unblinded interim result is promising. *Statistics in Medicine*.

Cui, L., Hung, H. M. J., & Wang, S. J. (1999). Modification of sample size in group sequential trials. *Statistics in Medicine*.

Dohmke, T. & Kalliamvakou, E. (2025). **Developers, Reinvented.** GitHub blog.  
`https://ashtom.github.io/developers-reinvented`

Jennison, C. & Turnbull, B. W. (2000). *Group Sequential Methods with Applications to Clinical Trials*. Chapman & Hall/CRC.

Lachin, J. M. & Foulkes, M. A. (1986). Evaluation of sample size and power for analyses of survival with allowance for nonuniform patient entry, losses to follow-up, noncompliance, and stratification. *Biometrics*.

Lan, K. K. G. & DeMets, D. L. (1989). Discrete sequential boundaries for clinical trials. *Biometrika*.

Lehmacher, W. & Wassmer, G. (1999). Adaptive sample size calculations in group sequential trials. *Biometrics*.

Mehta, C. R. & Pocock, S. J. (2011). Adaptive increase in sample size when interim results are promising: practical guidance and theoretical issues. *Clinical Trials*.

Proschan, M. A. & Hunsberger, S. A. (1995). Designed extension of studies based on conditional power. *Biometrics*.

Zhu, H. & Lakkis, H. (2014). Sample size calculation for comparing two negative binomial rates in noninferiority and equivalence trials with variable follow-up. *Statistics in Medicine* 33, 376–387.

---

## Appendix A. Mapping to software functions (non-exhaustive)

| Concept | Implementation |
|--------|----------------|
| MC validation of `nSurv()$d` vs simulated trials | **gsDesignNB** repo: `papers/simulation_sample_size_accuracy.R`; test `test-independent-test-simulation-accuracy.R` |
| Patient-level \(\delta\) between counted events | *Simulation post-processing*; apply Section 2.4 filter to simulated streams before summing trial events |
| Piecewise expected events | `eEvents1()`, `nEventsIA()` (**gsDesign**) |
| Fixed survival sample size | `nSurv()` (**gsDesign**) |
| Group sequential survival | `gsSurv()`, `gsSurvCalendar()` (**gsDesign**) |
| Power / timing with floors & gaps | `gsSurvPower()` (**gsDesign**) |
| Gapped spending | `sfGapped()` (**gsDesign**) |
| SSR / conditional power | `ssrCP()`, `Power.ssrCP()`, `condPower()` (**gsDesign**) |
| NB sample size with event gaps, recurrent / seasonal simulation | `sample_size_nbinom()`, `nb_sim()`, `nb_sim_seasonal()` (**gsDesignNB**) |
| Mutze Wald test (NB with Poisson fallback by \(\theta\) range) | `mutze_test()` (**gsDesignNB**) |
| NB GSD simulation (ML vs MoM information paths) | `sim_gs_nbinom()` (**gsDesignNB**) |
| Blinded SSR / blinded information | `blinded_ssr()`, `calculate_blinded_info()` (**gsDesignNB**) |

---

## Appendix B. Suggested journal fit

Possible outlets: *Statistics in Medicine*, *Biometrics*, *Pharmaceutical Statistics*, *Clinical Trials*, or a **software-focused** track (e.g., *Journal of Statistical Software*) if emphasis shifts to implementation. The present draft mixes methods and software; editors may request splitting.

---

*End of draft.*
