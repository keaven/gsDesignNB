# Simulation studies (`papers/`)

Manuscript materials for **gsDesignNB** (draft: `manuscript-gaps-events-sample-size.md`).

## Sample size accuracy (`simulation_sample_size_accuracy.R`)

This script evaluates **Monte Carlo agreement** between:

- **Analytical** total expected events under the alternative, `gsDesign::nSurv()$d` (Lachin–Foulkes / `eEvents()` machinery in **gsDesign**), and  
- **Simulated** total failure counts from independent trial replicates under the **same** model: uniform enrollment on `[0, \sum R]`, proportional hazards by arm, exponential dropout, failures counted only if they occur by calendar time `T`.

Run from the **gsDesignNB** repository root (with **gsDesign** installed):

```sh
Rscript papers/simulation_sample_size_accuracy.R
```

Optional environment variables:

- `SIM_NSIM` — number of Monte Carlo trials (default `2000`).

Example:

```sh
SIM_NSIM=5000 Rscript papers/simulation_sample_size_accuracy.R
```

The script also defines `thin_interevent_gap()` for the **patient-level minimum gap** rule described in the manuscript (Section 2.4), for use with recurrent-event simulators.

**Piecewise hazards** (`S` not `NULL`): the default simulator targets the **scalar-hazard** case. Extending the simulator to piecewise failure rates from enrollment is straightforward but not duplicated here; analytical checks still use `nSurv()$d` / `eEvents()`.

Automated regression: `tests/testthat/test-independent-test-simulation-accuracy.R`.

**Authors (manuscript):** Keaven M. Anderson, Hongtao Zhang, Andrea Maes (Merck & Co., Inc., Rahway, NJ, USA).
