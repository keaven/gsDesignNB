# Paper replication

Run the lightweight manuscript reproduction from the package root:

```sh
Rscript data-raw/reproduce-paper.R
```

The script writes the worked-design summary and boundaries, the Jensen-correction
table, and the two 50-replicate workflow examples to `paper/replication-output/`.
The SSR example is run after the group sequential example because it uses the
planning objects defined by that workflow.

The planned-calendar-time information validation study is separate from the
lightweight reproduction and routine tests. Run a small workflow check with:

```sh
Rscript data-raw/gs_validation_study.R --debug
```

Run the production study explicitly with:

```sh
Rscript data-raw/gs_validation_study.R --production
```

The production run uses 20,000 null and 5,000 alternative replicates per
scenario and timing method. At every planned analysis date it compares mean
simulated information with expected information under the generating rates,
dispersion, accrual, dropout, and follow-up. Under the design alternative,
this expectation equals the design's planned information; under null rates,
the expectation is recomputed rather than assuming the alternative's target.
All replicates contribute at every date, without boundary-based stopping.
The summary reports each ML and MoM estimator separately, Monte Carlo standard
errors and intervals for its mean, and counts of finite and failed estimates.
It does not assess boundary crossings, Type I error, or power.

The script writes `information-summary.csv`, `gs-information-validation.rds`,
and `session-info.txt` to `paper/gs-validation-results/` by default. Set
`--output-dir=/path/to/output` to choose another location. Each completed
scenario cell is saved under `cells/` in the output directory; rerunning the
same command reuses those cells and runs only the missing ones, so an
interrupted run can be resumed. Information-only cells use the prefix
`information-v1__` and do not reuse earlier crossing-calibration caches.

To use another output location:

```sh
Rscript data-raw/reproduce-paper.R --output-dir=/path/to/output
```

The manuscript's long Monte Carlo studies are supplied as versioned RDS caches
in `inst/extdata/`; the default script verifies their presence and records the
commands that regenerate them. Those commands are deliberately not run by
default because the production studies use up to 20,000 replicates per scenario.
