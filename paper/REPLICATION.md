# Paper replication

Run the lightweight manuscript reproduction from the package root:

```sh
Rscript data-raw/reproduce-paper.R
```

The script writes the worked-design summary and boundaries, the Jensen-correction
table, and the two 50-replicate workflow examples to `paper/replication-output/`.
The SSR example is run after the group sequential example because it uses the
planning objects defined by that workflow.

To use another output location:

```sh
Rscript data-raw/reproduce-paper.R --output-dir=/path/to/output
```

The manuscript's long Monte Carlo studies are supplied as versioned RDS caches
in `inst/extdata/`; the default script verifies their presence and records the
commands that regenerate them. Those commands are deliberately not run by
default because the production studies use up to 20,000 replicates per scenario.
