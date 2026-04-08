---
name: gsdesignnb-simulation-entrypoint
description: Route gsDesignNB simulation requests to the right package functions and vignettes. Use when working on recurrent-event simulation, group sequential simulation, sample size re-estimation, blinded information, or when the user needs help choosing between nb_sim, sim_gs_nbinom, sim_ssr_nbinom, and related workflows.
---
# gsDesignNB Simulation Entry Point

## Quick routing

Choose the smallest package entry point that matches the user request:

- **Single recurrent-event data set**:
  Use `nb_sim()`.
- **Seasonal enrollment or event rates**:
  Use `nb_sim_seasonal()`.
- **Calendar-based group sequential operating characteristics**:
  Use `sim_gs_nbinom()`, then `check_gs_bound()` and `summarize_gs_sim()`.
- **Adaptive SSR with information-based interims**:
  Use `sim_ssr_nbinom()`, then `summarize_ssr_sim()`.
- **One interim cut from an existing data set**:
  Use `cut_data_by_date()`, then `blinded_ssr()`, `unblinded_ssr()`, `calculate_blinded_info()`, or `mutze_test()`.

## Time-scale checklist

Before writing or reviewing code, confirm that all of these use the same time unit:

- `enroll_rate$duration`
- `fail_rate$rate`
- `dropout_rate$rate` and `dropout_rate$duration`
- `trial_duration`
- `max_followup`
- calendar cut times / analysis times
- `event_gap`

If units are mixed, fix that first.

## Recommended user path

For most user-facing work:

1. Build a fixed design with `sample_size_nbinom()`.
2. Convert to a GS design with `gsNBCalendar()` and usually `toInteger()`.
3. Pick the simulator:
   - `sim_gs_nbinom()` for fixed GS timing
   - `sim_ssr_nbinom()` for adaptive SSR
4. Summarize with:
   - `summarize_gs_sim()`
   - `summarize_ssr_sim()`

## SSR-specific guidance

When the request mentions **blinded vs unblinded SSR**, **adaptive sample size**, or
**information-based interim timing**:

- default to `sim_ssr_nbinom()`
- mention that it returns:
  - `trial_results`
  - `analysis_results`
  - expected participants with events
  - expected events observed
- use `bound_info` deliberately:
  - default is `unblinded_ml`
  - alternatives include blinded / MoM information

## Common pitfalls

- Do not re-implement the SSR engine inside a vignette or notebook if `sim_ssr_nbinom()` covers the workflow.
- `mutze_test()` uses the Mutze-style Wald model and may fall back to Poisson depending on `poisson_threshold`.
- `check_gs_bound()` can now use an explicit `info_col`; use that instead of renaming columns.
- `event_gap` changes both event counts and exposure at risk; do not treat it as only a plotting artifact.

## Best vignette targets

- `vignettes/simulation-example.Rmd`: basic recurrent-event simulation
- `vignettes/group-sequential-simulation.Rmd`: calendar-based GS workflow
- `vignettes/ssr-example.Rmd`: short SSR workflow
- `vignettes/ssr-simulation-study.Rmd`: advanced SSR case study
- `vignettes/blinded-info-diagnostics.Rmd`: edge cases for blinded information
