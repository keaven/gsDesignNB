# Bootstrap-resample subjects within strata (internal)

Performs stratified cluster (subject-level) bootstrap resampling. When a
subject is drawn more than once its rows are duplicated and the
`subject_col` values are replaced with unique pseudo-IDs so that the
GLMM treats them as separate individuals (mirroring SAS
`PROC SURVEYSELECT` with `method=urs cluster=USUBJID`).

## Usage

``` r
.bootstrap_by_cluster(data, subject_col, strata_cols, n_boot)
```

## Arguments

- data:

  Data frame to resample.

- subject_col:

  Character. Column that identifies subjects (clusters).

- strata_cols:

  Character vector. Stratification column names. May be `character(0)`
  for no stratification.

- n_boot:

  Integer. Number of bootstrap replicates.

## Value

List of length `n_boot`. Each element is a copy of `data` with an added
`replicate` integer column and (when `n_boot > 1`) potentially modified
`subject_col` values for duplicated subjects.
