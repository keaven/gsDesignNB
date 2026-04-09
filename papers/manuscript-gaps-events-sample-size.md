# Negative Binomial Recurrent-Event Group Sequential Designs with Calendar Gaps, Event Spacing, and Sample Size Re-estimation

**Draft manuscript** (target length: 10–15 journal pages; ~6,000 words)

**Authors:** Keaven M. Anderson, Hongtao Zhang, Andrea Maes  

**Affiliation:** Merck & Co., Inc., Rahway, NJ, USA

**Keywords:** negative binomial; recurrent events; group sequential design; sample size re-estimation; information time; software

---

## Abstract

Negative binomial recurrent-event trials are often monitored using a mixture of calendar-time rules, minimum follow-up, and information-based interim criteria. In **gsDesignNB**, the working recurrent-event model is simple: conditional on treatment arm and a subject-level gamma random effect, recurrent events follow a Poisson process with constant within-subject rate, yielding marginal negative binomial counts. Naive planning that ignores **mandatory gaps between successive interim looks**—for example, minimum time from one analysis to the next while events continue to accrue—can materially misstate expected information, power, and spending paths. We also treat **patient-level minimum gaps between successive events** (e.g. at least 20 days between one qualifying event and the next on the same patient): in simulation, any event that would fall inside that window is **discarded** and does **not** count toward protocol-defined totals. We describe a **timing and information framework** that combines this recurrent-event model with piecewise enrollment, dropout, and supporting expected-event calculations; show how **calendar floors** and **minimum spacing between analyses** enter the timing solution; and connect these calculations to **gsDesignNB** workflows for negative binomial sample size, simulation, information monitoring, and **sample size re-estimation** (SSR). Because statistical information depends jointly on exposure, event rates, and dispersion, it is **not** directly proportional to raw event counts. We summarize **simulation-based verification** for fixed designs and SSR, including the **event-gap correction** needed to align planned power with Monte Carlo under recurrent-event counting rules, scenarios with **seasonal** event rates, and comparison of **maximum-likelihood** versus **method-of-moments** information paths. In the SSR vignette configuration studied here, nominal one-sided \(\alpha = 0.025\) showed slight inflation, a smaller nominal alpha improved calibration, and blinded updates tended to request larger adapted sample sizes than unblinded updates. The **gsDesign** package remains the sequential-design backend for spending and boundary calculations. Finally, we note that **AI-assisted development** helped implement and validate a broader software package around these methods, with verification emphasizing large-scale simulation and regression testing [Dohmke & Kalliamvakou, 2025](https://ashtom.github.io/developers-reinvented).

---

## 1. Introduction

Negative binomial models are natural for recurrent-event endpoints when subjects may experience multiple overdispersed events and follow-up varies across patients. In **gsDesignNB**, the default generative story is simple: conditional on treatment arm and a subject-specific gamma frailty, each patient follows a Poisson process with a constant event rate over follow-up, so repeated exponential waiting times produce recurrent events while the marginal event count distribution is negative binomial. The package also supports **season-specific** event rates as an explicit time-varying extension. In confirmatory development, these trials still inherit many operational rules familiar from event-driven monitoring: interim analyses tied to **calendar dates**, **minimum patient exposure**, target information thresholds, and **minimum elapsed time** between successive database locks while events continue to accrue. If those constraints are not reflected in planned information growth, the resulting sample size, spending path, and final power can be materially misstated.

We use the phrase **“gaps between events”** in three complementary senses that arise in applications:

1. **Operational / calendar gaps between interim analyses:** successive looks must occur no earlier than \(T_{i-1}+\Delta_i\) after the previous analysis, possibly concurrent with additional accrual and ongoing event generation.
2. **Structural gaps in the timing infrastructure:** piecewise accrual, dropout, and related supporting expected-event calculations imply nonuniform accumulation of exposure and events.
3. **Patient-level minimum gap between successive events:** a protocol may require that **only** events separated by at least \(\delta\) **patient time** (e.g. \(\delta = 20\) days) count toward the analysis event total. **Simulated** events that fall inside \((t_{\text{prev}},\, t_{\text{prev}}+\delta]\) after a previous qualifying event at \(t_{\text{prev}}\) for the same patient are **excluded** from the count—they do not contribute to information. This rule mimics adjudication windows, avoidance of double-counting clinically indistinguishable events, or operational definitions of “new” progression.

The methodological contribution of this paper is not a single new asymptotic result but a **clear definition and derivation** of the timing and information calculations that underpin sample size when such gaps are specified, and their integration with **alpha/beta spending** and **SSR** tools for recurrent-event endpoints. The primary package story is the negative binomial recurrent-event model above: simple exponential event generation within subject, gamma frailty across subjects, and information for the log rate ratio determined jointly by exposure, rates, and dispersion. Supporting **closed-form** expected-event formulas (Section 2) remain useful for timing infrastructure and repository validation under simpler event-time assumptions; **patient-level inter-event gaps** are easiest to incorporate rigorously via **simulation** (Section 5), where the counting rule is applied explicitly. In software, **gsDesignNB** is the user-facing layer for negative binomial sample size, recurrent-event simulation, information monitoring, **Mutze-type** interim testing, and blinded or unblinded SSR. The **gsDesign** package remains the supporting backend for the sequential-design pieces: spending functions, correlation calculations, and boundary inversion.

**Outline.** Section 2 sets notation for the recurrent-event model, then reviews supporting expected-event infrastructure and **per-patient inter-event spacing**. Section 3 connects calendar and information-based rules to interim timing. Section 4 shows how those calculations feed negative binomial group sequential monitoring while keeping **gsDesign** in a backend role. Section 5 outlines simulation studies: **§5.1** verifies the NB fixed-design calculation and the event-gap correction; later subsections cover **negative binomial** operating characteristics, SSR, and regression testing. Section 6 describes AI’s role in implementing and validating the software. Section 7 concludes.

---

## 2. Recurrent-event model and timing infrastructure

### 2.1 Recurrent-event model and notation

Within **gsDesignNB**, the default recurrent-event model assumes that, conditional on treatment arm, each subject follows a Poisson process with a constant event rate over follow-up. Subject-level heterogeneity is represented by a gamma frailty random effect, so the marginal event count over exposure time \(t\) follows a negative binomial distribution. Let \(\lambda_C\) and \(\lambda_E\) denote the control and experimental mean event rates, and let \(\mathrm{RR}=\lambda_E/\lambda_C\) denote the treatment **rate ratio**.

For timing calculations, enrollment may be **piecewise constant** over deterministic periods of length \(R_j\), dropout may be **piecewise exponential**, and maximum follow-up may truncate exposure at calendar time \(T\). The package's main built-in time-varying extension for the event rate itself is **seasonality**, implemented through season-specific base rates in `nb_sim_seasonal()`.

When we refer to the supporting survival-style expected-event infrastructure used in Section 2.2, we use standard piecewise exponential notation with control hazard \(\lambda_C(t)\), experimental hazard \(\lambda_E(t)\), and hazard ratio **HR**.

### 2.2 Supporting expected-event calculations (derivation sketch)

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

The **total expected events** for the trial is obtained by summing over arms and strata (independent strata), and scaling by randomization ratio. In **gsDesignNB**, analogous exposure bookkeeping is paired with rate-model information calculations for negative binomial recurrent-event designs; total events remain useful summaries, but they are not by themselves the design target.

### 2.3 Supporting fixed-\(L\) event calculation and NB target information

As a supporting comparator, the **fixed-design** event requirement \(d\) under a PH alternative solves the standard Schoenfeld-style relationship using the **log-rank** asymptotics, with non-centrality driven by \(\log(\text{HR})\) and the expected event count. The negative binomial sample-size calculation implemented in **gsDesignNB**, however, follows **Zhu and Lakkis (2014, Method 3)** and shares the same large-sample Wald variance structure used by **Friede and Schmidli (2010)** for blinded SSR and by **Mutze et al. (2019)** for group sequential monitoring.

For the negative binomial log rate-ratio estimator \(\hat\theta\), statistical information is \(\mathcal{I}=1/\mathrm{Var}(\hat\theta)\). For subject \(i\) in group \(g\) with exposure \(t_i\) and expected count \(\mu_i = \lambda_g t_i\), the per-subject Fisher information contribution is
\[
\mathcal{I}_i = \frac{\mu_i}{1 + k_g \mu_i}.
\]
Aggregating over subjects gives \(W_g = \sum_{i \in g} \mu_i/(1 + k_g \mu_i)\) and
\[
\mathcal{I} = \frac{1}{1/W_1 + 1/W_2}.
\]
When exposures are homogeneous within arm, this reduces to
\[
\mathrm{Var}(\hat\theta) = \frac{1/\mu_1 + k_1}{n_1} + \frac{1/\mu_2 + k_2}{n_2}.
\]
For balanced allocation, the corresponding fixed-design approximation can be written as
\[
n_{\text{total}} \approx
\frac{(z_{1-\alpha^\star}+z_{1-\beta})^2
\left[\left(\frac{1}{\mu_1}+k_1\right)+\left(\frac{1}{\mu_2}+k_2\right)\right]}
{(\theta-\theta_0)^2},
\]
where \(\theta=\log(\mathrm{RR})\), \(\theta_0\) is the null value, and \(\alpha^\star\) denotes the one-sided or two-sided tail area convention.
Accordingly, target information depends jointly on exposure, event rates, and overdispersion; it is **not** directly proportional to the raw number of observed events. Two analyses with similar event totals can therefore carry different information if follow-up or dispersion differs.
This is a practical distinction from usual first-event formulas: when a protocol imposes an event gap \(g\), the expected count \(\mu_g\) must reflect the **gap-corrected effective rate** rather than a naive event count proxy.

### 2.4 Patient-level minimum gap between successive events (inter-event window)

Some protocols define the **analysis event count** using a **minimum separation** between successive **qualifying** events on the **same patient**. Let \(\delta > 0\) denote that minimum (e.g. **20 days** on the study or calendar timeline from randomization, or from the previous event—definitions should match the SAP). Suppose patient \(j\) would generate a **sequence** of potential event times \(T_{j,1} < T_{j,2} < \cdots\) under a simulation model (e.g. from a semi-Markov or recurrent-event process, or from multiple draws in a composite endpoint).

**Counting rule (simulation).** Sort potential event times for patient \(j\) as \(T_{j,1} < T_{j,2} < \cdots\). Let \(C_{j,1}=T_{j,1}\) be the first **counted** time. Recursively: the next **counted** time is the **first** \(T_{j,r}\) such that \(T_{j,r} \ge C_{j,s} + \delta\), where \(C_{j,s}\) is the previous counted time; all **earlier** \(T_{j,r}\) in \((C_{j,s},\, C_{j,s}+\delta)\) are **discarded**—they **do not count** toward the trial’s cumulative event total. Equivalently: after a counted event at calendar time \(t\), **no** simulated event strictly before \(t+\delta\) can become the next counted event.

This rule **thins** recurrent or densely simulated events relative to a model with no minimum gap. The **operational** information accrual is therefore **slower** than the naive expected event total from Section 2.2 would suggest whenever multiple events per patient are possible under the generative model. **Analytical** adjustment requires a recurrent-events model with **gap time** constraints (e.g. gap-time Cox or frailty formulations); a practical approach for design is to **calibrate** target information, \(d\), or enrollment inflation by **Monte Carlo** under the \(\delta\)-rule until targeted power is met.

**Example.** If \(\delta = 20\) days and the simulator proposes a second event 12 days after the first, that second event **does not count**. Only events occurring **at least 20 days** after the previously **counted** event on that patient can enter the information denominator.

---

## 3. Calendar constraints and “gaps” between successive interim analyses

### 3.1 Information versus calendar time

In supporting event-driven designs, the **information fraction** at analysis \(i\) is often proxied by cumulative events divided by the planned final events. For **negative binomial recurrent-event** designs, that proxy is incomplete: statistical information for the log rate ratio depends jointly on exposure, event rates, and dispersion. Events therefore contribute to information growth, but they are not its sole determinant. When interim analyses are **also** constrained by calendar rules, the realized **sequence \((T_1,\ldots,T_K)\)** of analysis times is the solution of **feasibility inequalities**:

- **Planned calendar time** (if any): \(T_i \ge C_i\).
- **Minimum gap** after the previous analysis: \(T_i \ge T_{i-1}+\Delta_i\) (months or days on study clock).
- **Minimum exposure:** enough patients with \(\ge\) minimal follow-up.
- **Information threshold:** if the \(i\)-th look is triggered by **target information** \(\mathcal{I}_i\), then \(T_i\) is the **smallest** time \(t\) such that **expected** (or observed) information \(\ge \mathcal{I}_i\), subject to the floors above. In supporting event-driven approximations, event thresholds \(d_i\) may be used as proxies.

This matches the backend **priority/floor** timing logic: take the **maximum** of applicable **floor times** (calendar, gap-from-previous, enrollment+follow-up), then adjust upward if target information—or, in supporting survival-style approximations, event targets—requires a later time, optionally capped by **maxExtension** rules.

### 3.2 Why ignoring gaps biases sample size

If one assumes interim analyses occur at fixed fractions of final events—or other crude event-count proxies—**without** enforcing \(T_i \ge T_{i-1}+\Delta_i\), the **expected information path** is systematically misspecified whenever \(\Delta_i>0\) and enrollment continues between looks: **events accrue during the gap**, but the information gained also depends on follow-up and dispersion. The proposed **timing + floor** construction aligns the information-time vector fed to `gsDesign()` with trial operations.

### 3.3 Interaction with patient-level inter-event gaps (Section 2.4)

**Trial-level** gaps \(\Delta_i\) between database locks (Section 3.1) and **patient-level** gaps \(\delta\) (Section 2.4) operate on **different** time scales: the former constrain **when** looks occur; the latter constrain **which** patient-level events increment the event tally and therefore the accumulated information. A design that sets information targets assuming every simulated failure counts will be **overpowered** or **under-accrue** information if the SAP applies a **δ-rule** that reduces the realized event rate. Reporting simulations **with and without** \(\delta\)-thinning is therefore essential when recurrent or composite events are modeled.

---

## 4. Group sequential monitoring and backend spending machinery

### 4.1 Error spending on information time

Let \(t \in [0,1]\) denote **information fraction**. For negative binomial recurrent-event designs, this is naturally based on the inverse variance of the log rate-ratio estimator and is only approximately related to event fraction because exposure and dispersion matter. Standard spending functions \(f_\alpha(t)\) allocate cumulative Type I error; analogs allocate \(\beta\) for futility under the alternative. Given **timing** \(t_i\) at look \(i\), group sequential boundaries are computed by inverting spending constraints jointly with the correlation structure of the sequential test statistics (e.g. log-rank or log rate-ratio statistics).

### 4.2 Putting the NB workflow together

A complete **negative binomial recurrent-event** group sequential workflow combines:

1. A recurrent-event model with simple exponential within-subject event generation, gamma frailty, and when needed **patient-level event-spacing** rules (Section 2).
2. **Timing** \(T_i\) from Section 3 (calendar gaps + target information).
3. **Negative binomial** design and information calculations (`sample_size_nbinom()`, `compute_info_at_time()`, `calculate_blinded_info()`).
4. **Spending** \(f_\alpha, f_\beta\) and related boundary calculations supplied by the **gsDesign** backend.
5. Inflation of enrollment or trial duration to meet power at the final analysis.

Within the broader ecosystem, survival-oriented timing tools remain useful as supporting comparators. For recurrent-event count endpoints, however, **gsDesignNB** wraps the same backend boundary logic around rate-ratio information and recurrent-event simulation.

---

## 5. Simulation studies

We recommend reporting **three layers** of simulation-related evidence in a submission: (0) verification that the **NB sample-size calculation** matches Monte Carlo under the documented model; (1) broader **negative binomial recurrent-event** operating characteristics; (2) SSR.

### 5.1 Simulation verification of NB sample size and event-gap correction

The fixed-design calculation in `sample_size_nbinom()` should be checked directly against Monte Carlo under the same recurrent-event data-generating model used in planning. The key validation target is not just the average exposure calculation but also the **event-gap correction** when a protocol requires a dead time \(g\) after each counted event. Without that correction, the naive renewal approximation \(\lambda/(1+\lambda g)\) overstates the effective event rate under gamma frailty and therefore makes the planned design look slightly more powerful than it really is.

Using the random frailty \(\Lambda\) from the gamma--Poisson mixture, the relevant population quantity is \(\mathrm{E}[\Lambda/(1+\Lambda g)]\), not \( \lambda/(1+\lambda g)\). A second-order correction gives
\[
\mathrm{E}\!\left[\frac{\Lambda}{1+\Lambda g}\right]
\approx
\frac{\lambda}{1+\lambda g}
\left(1 - \frac{k\lambda g}{(1+\lambda g)^2}\right),
\]
which is the adjustment implemented in the package. During package development, this modification was needed for **simulation-based power** to agree with the analytical planning calculation when event gaps were non-negligible.
In that sense, the recurrent-event sample-size formula differs from more familiar first-event formulas not only by using NB rate-ratio information, but also by requiring the event-gap correction when the protocol counting rule induces dead time after events.

The sample-size vignette documents this verification in two complementary ways: by comparing theoretical and simulated average exposure under piecewise accrual/dropout, and by checking that the corrected NB planning formula reproduces the intended power in Monte Carlo. It also includes graphical illustrations of the event-count distribution as dispersion increases, emphasizing why information per expected event decreases as \(k\) grows.

### 5.2 Monte Carlo operating characteristics (fixed design / GSD)

- **Data generating process:** simple exponential recurrent-event generation with gamma frailty as the base case; optional **seasonality** for calendar-time variation; and **recurrent or multi-event** generators when patient-level \(\delta\) is relevant.
- **Seasonality:** when event intensity is expected to vary over calendar seasons (for example, respiratory or infectious-disease endpoints), **gsDesignNB** supports direct Monte Carlo generation with season-specific rates via `nb_sim_seasonal()`, after which the same cut-and-analysis workflow can be applied.
- **Patient-level \(\delta\)-rule:** after drawing event times per patient, apply the **Section 2.4** filter so events within \(\delta\) of the previous **counted** event do not count; accumulate trial-level events from the thinned streams.
- **Outputs:** empirical Type I error (where feasible), power, distributions of interim analysis times, blinded and unblinded **information fractions**, distribution of cumulative spending paths, early stopping rates, and **distribution of time to reach target information** with versus without thinning.
- **Comparison:** designs **with** versus **without** enforced minimum gaps between looks; **with** versus **without** event-gap correction when \(g>0\); **\(\delta = 0\)** versus **\(\delta = 20\)** days (or other SAP values).

### 5.3 Sample size re-estimation (SSR)

Classical conditional-power formulations such as `ssrCP()` are useful reference points for two-stage adaptation, but the practical package workflow here is the **negative binomial recurrent-event** layer in **gsDesignNB**. For blinded negative binomial SSR, **Friede and Schmidli (2010)** and **Schneider, Schmidli, and Friede (2013)** are the closest methodological reference points. Simulations from `sim_ssr_nbinom()` should report:

- unconditional power,
- distributions of adapted total enrollment or extended follow-up,
- blinded versus unblinded nuisance-parameter estimates at the adaptation,
- control of Type I error under the pre-specified boundary rule.

These simulations complement **analytical** expected-event calculations: they validate **behavior under model misspecification** and **discrete** recurrent-event accrual.

For **negative binomial recurrent-event** designs, **gsDesignNB** extends the same idea: interim analyses are naturally keyed to **statistical information** (variance of the log rate-ratio under a working NB model), not only to calendar milestones. **Sample size re-estimation** is most valuable when planning assumptions for **nuisance parameters**—especially **overdispersion** and accrual or follow-up intensity—are uncertain; updating the projected total enrollment or follow-up after an interim look limits unconditional power loss without rewriting the group-sequential boundary logic, provided the adaptation rule is **pre-specified** and operating characteristics are simulated. In this setting, target information depends jointly on exposure, observed event rates, and dispersion, so raw event totals alone do not determine how much evidence has accrued at the interim.

In the SSR vignette configuration developed for the package, nominal one-sided \(\alpha = 0.025\) showed **slight Type I inflation** under the planned SSR rule, consistent with the cautionary findings reported by **Mutze et al. (2019)** for adaptive group sequential settings. Rebuilding the design at a slightly smaller nominal alpha provided a practical improvement in calibration for that example. This should be viewed as an **empirical design-specific sensitivity result**, not as a universal correction.

**Blinded versus unblinded nuisance updates** trade scientific and regulatory goals. **Blinded** SSR pools arms to estimate a common dispersion and pooled rate, which preserves masking of the treatment effect at the adaptation; in simulation studies reported in the package, blinded re-estimation can yield **larger** adapted sample sizes than **unblinded** updates that use arm-specific fits, because the blinded working model is intentionally conservative about heterogeneity. Where a data monitoring charter and regulators accept **unblinded** nuisance-parameter updates with appropriate error-spending or combination-test discipline, that path is often **more sample-efficient** and may therefore be preferable in practice; the choice should be justified with **simulated Type I error and power** under the intended SSR rule, not by analytical shortcut alone.

### 5.4 Illustrative simulation factors (reporting template)

| Factor | Levels (example) |
|--------|-------------------|
| Rate ratio under alternative | 0.65, 0.75, 0.85 |
| Dispersion \(k\) | planned vs inflated |
| Minimum gap \(\Delta\) between looks | 0, 3, 6 (months) |
| Patient-level min gap \(\delta\) between counted events | 0, 20 (days) |
| Enrollment slowdown | nominal \(\gamma\) vs 80% of nominal |
| Seasonality | none vs winter peak |
| SSR rule | `sim_ssr_nbinom()` with blinded or unblinded updates and capped increase |

For each configuration: **10,000+** trial replicates (more for Type I near \(\alpha\)). Summarize **median** and **IQR** of calendar times at looks; **empirical power**; **mean** stage-2 \(n\) under SSR.

### 5.5 Software regression as a complement to Monte Carlo

Beyond standalone simulation scripts, the **gsDesignNB** package and its supporting materials use **automated** tests (e.g., `testthat`, snapshot tests for plots where applicable) to catch **drift** in computations and plotting when refactoring. These are **not** substitutes for methodological simulation but **guardrails** ensuring released code matches stored reference outputs—especially important when AI-assisted edits increase commit velocity.

---

## 6. Role of AI in methods development, coding, and verification

Recent accounts of **AI-experienced developers** emphasize a shift from **line-by-line coding** to **delegation** (framing tasks, curating context, reviewing plans) and **verification** (testing, review, security) [Dohmke & Kalliamvakou, 2025](https://ashtom.github.io/developers-reinvented). That perspective aligns with statistical software for regulated contexts: **correctness** is established by **independent checks**, not by authority of the generator.

In a **gsDesignNB**-centered workflow, AI tools can accelerate implementation of the methods and the surrounding software package by helping with:

- algebra and bookkeeping for recurrent-event information, accrual, and timing calculations;
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

1. **NB-first workflow:** The main public story is **negative binomial recurrent-event** design and monitoring in **gsDesignNB**; **gsDesign** is the supporting backend for spending functions and boundary calculations.

2. **Model and information target:** The default recurrent-event model is simple exponential within subject with gamma frailty across subjects, and the design target is **statistical information**, which depends jointly on exposure, event rates, and dispersion rather than on raw event totals alone.

3. **Sample size details:** The fixed-design NB formula is not just a relabeling of survival methods; it uses rate-ratio information, and when protocol event gaps are present the **event-gap correction** is needed for Monte Carlo power to agree with planning.

4. **Simulation and SSR:** Operating characteristics should be reported for both **fixed GSD** and **SSR** pipelines, including recurrent-event features such as **seasonality** or event thinning; these are the principal **external** checks on complex code paths. In the vignette example, blinded SSR tended to request larger adapted sample sizes than unblinded SSR, and a slightly smaller nominal alpha improved Type I calibration.

5. **Mutze-type interim tests and dispersion:** For recurrent **count** outcomes, the **Mutze et al.** Wald test implemented in **gsDesignNB** as `mutze_test()` fits a negative binomial log-linear model with an offset for exposure. When maximum-likelihood estimation implies **essentially Poisson** variation (very large estimated NB shape \(\theta\)) or **numerically unstable** extreme overdispersion (very small \(\theta\)), the implementation **falls back to Poisson** regression for the Wald statistic. The switch is governed by `poisson_threshold` \(T\) (default \(T=1000\)): Poisson is used if \(\theta \notin [1/T, T]\). **Lowering \(T\)** triggers the Poisson fallback **sooner**—i.e., whenever you prefer to treat the data as **without extra-Poisson dispersion** once \(\theta\) is only moderately large—and is a practical way to align the interim model with a **Poisson** SAP without forcing `method = "poisson"` for every analysis. **Raising \(T\)** keeps the NB fit over a wider \(\theta\) range. Sensitivity analyses across \(T\) or explicit `method = "poisson"` versus `"nb"` should accompany any change from package defaults.

6. **Method of moments and delivery layers:** Group-sequential simulation in **gsDesignNB** (`sim_gs_nbinom()`) records **maximum-likelihood** and **method-of-moments** dispersion estimates and the associated blinded and unblinded **information** paths. If ML-based nuisance estimates behave poorly at small interim sample sizes, comparing MoM and ML paths in simulation supports whether MoM (or blinded helpers that fall back when a blinded NB fit fails) is preferable for **SSR** or information monitoring. Validated package APIs should remain primary; optional **Shiny** and **AI** layers are useful only when they sit above the same tested statistical core.

7. **Limitations:** The supporting piecewise exponential timing framework is an approximation; deviations (delayed effects, cure fractions) require extended models or more aggressive simulation. Calendar rules in real trials involve **random** accrual jitter and operational delays not fully captured by deterministic \(\Delta_i\). **Patient-level** \(\delta\)-rules are handled straightforwardly in **simulation**; closed-form **expected counted events** under general recurrent processes with hard gaps is **not** included in standard PH fixed-design formulas and may require dedicated recurrent-event theory or empirical calibration.

---

## Acknowledgments

The authors thank colleagues at Merck & Co., Inc., Rahway, NJ, USA for support and feedback.

---

## References

Anderson, K. et al. (ongoing). **gsDesignNB**: Sample Size and Simulation for Negative Binomial Outcomes. R package.  
URL: `https://github.com/keaven/gsDesignNB` and `https://keaven.github.io/gsDesignNB/`.

Anderson, K. et al. (ongoing). **gsDesign**: Group Sequential Design. R package.  
URL: `https://github.com/keaven/gsDesign` and `https://keaven.github.io/gsDesign/`.

Bauer, P. & Köhne, K. (1994). Evaluation of experiments with adaptive interim analyses. *Biometrics*.

Chen, Y. H., DeMets, D. L., & Lan, K. K. G. (2004). Increasing the sample size when the unblinded interim result is promising. *Statistics in Medicine*.

Cui, L., Hung, H. M. J., & Wang, S. J. (1999). Modification of sample size in group sequential trials. *Statistics in Medicine*.

Dohmke, T. & Kalliamvakou, E. (2025). **Developers, Reinvented.** GitHub blog.  
`https://ashtom.github.io/developers-reinvented`

Friede, T. & Schmidli, H. (2010). Blinded sample size reestimation with negative binomial counts in superiority and non-inferiority trials. *Methods of Information in Medicine* 49, 618–624.

Jennison, C. & Turnbull, B. W. (2000). *Group Sequential Methods with Applications to Clinical Trials*. Chapman & Hall/CRC.

Lachin, J. M. & Foulkes, M. A. (1986). Evaluation of sample size and power for analyses of survival with allowance for nonuniform patient entry, losses to follow-up, noncompliance, and stratification. *Biometrics*.

Lan, K. K. G. & DeMets, D. L. (1989). Discrete sequential boundaries for clinical trials. *Biometrika*.

Lehmacher, W. & Wassmer, G. (1999). Adaptive sample size calculations in group sequential trials. *Biometrics*.

Mehta, C. R. & Pocock, S. J. (2011). Adaptive increase in sample size when interim results are promising: practical guidance and theoretical issues. *Clinical Trials*.

Mütze, T., Glimm, E., Schmidli, H., & Friede, T. (2019). Group sequential designs for negative binomial outcomes. *Statistical Methods in Medical Research* 28, 2326–2347.

Proschan, M. A. & Hunsberger, S. A. (1995). Designed extension of studies based on conditional power. *Biometrics*.

Schneider, S., Schmidli, H., & Friede, T. (2013). Blinded sample size re-estimation for recurrent event data with time trends. *Statistics in Medicine* 32, 5448–5457.

Zhu, H. & Lakkis, H. (2014). Sample size calculation for comparing two negative binomial rates in noninferiority and equivalence trials with variable follow-up. *Statistics in Medicine* 33, 376–387.

---

## Appendix A. Software map (NB-first, non-exhaustive)

| Concept | Implementation |
|--------|----------------|
| NB fixed sample size with event gaps | `sample_size_nbinom()` (**gsDesignNB**) |
| Recurrent-event simulation, including seasonal rates | `nb_sim()`, `nb_sim_seasonal()` (**gsDesignNB**) |
| NB calendar-aware group sequential design | `gsNBCalendar()` (**gsDesignNB**) |
| NB GSD simulation (ML vs MoM information paths) | `sim_gs_nbinom()`, `summarize_gs_sim()` (**gsDesignNB**) |
| NB SSR simulation and summaries | `sim_ssr_nbinom()`, `summarize_ssr_sim()` (**gsDesignNB**) |
| Blinded / unblinded SSR and information helpers | `blinded_ssr()`, `unblinded_ssr()`, `calculate_blinded_info()` (**gsDesignNB**) |
| Mutze Wald test (NB with Poisson fallback by \(\theta\) range) | `mutze_test()` (**gsDesignNB**) |
| Patient-level \(\delta\) between counted events | *Simulation post-processing*; apply Section 2.4 filter to simulated streams before summing trial events |
| Supporting event-time comparators and timing workflows | `eEvents1()`, `nEventsIA()`, `nSurv()`, `gsSurv()`, `gsSurvCalendar()`, `gsSurvPower()` (**gsDesign**) |
| Supporting conditional-power utilities | `ssrCP()`, `Power.ssrCP()`, `condPower()` (**gsDesign**) |

---

## Appendix B. Suggested journal fit

Possible outlets: *Statistics in Medicine*, *Biometrics*, *Pharmaceutical Statistics*, *Clinical Trials*, or a **software-focused** track (e.g., *Journal of Statistical Software*) if emphasis shifts to implementation. The present draft mixes methods and software; editors may request splitting.

---

*End of draft.*
