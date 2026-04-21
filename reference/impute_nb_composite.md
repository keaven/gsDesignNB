# Apply the composite ICE strategy: replace post-ICE outcomes with baseline

For subjects whose missingness flag matches `composite_value`, all
missing post-ICE count observations are set to the subject's baseline
count. This implements the composite estimand strategy for intercurrent
events such as death or treatment discontinuation due to disease
worsening, where the event itself is incorporated into the outcome
(e.g., baseline carried forward as a "worst case" placeholder).

## Usage

``` r
impute_nb_composite(
  data,
  outcome_col,
  imputed_value_col = "imputed_value",
  miss_flag_col,
  composite_value = "Comp",
  baseline_col
)
```

## Arguments

- data:

  Data frame.

- outcome_col:

  Character. Column with the original count outcome; used to identify
  which rows are missing (`NA`).

- imputed_value_col:

  Character. Column to update. If absent, it is created as a copy of
  `outcome_col` before the composite fill is applied. Default
  `"imputed_value"`.

- miss_flag_col:

  Character. Column with the missingness flag.

- composite_value:

  Character. Flag value triggering the composite strategy. Default
  `"Comp"`.

- baseline_col:

  Character. Column with the baseline count used as the fill value.

## Value

Data frame with `imputed_value_col` updated for composite rows.

## Details

The function is intentionally simple and requires no model. It can be
applied to a dataset already containing `imputed_value` from a prior MAR
or MNAR imputation step, or directly to the original data.

## Examples

``` r
df <- data.frame(
  count         = c(3L, NA,   NA,    5L),
  imputed_value = c(3L,  7L,  NA,    5L),
  miss_flag     = c(NA, "MAR", "Comp", NA),
  baseline      = c(4L,  4L,   4L,    6L)
)
impute_nb_composite(
  df,
  outcome_col     = "count",
  miss_flag_col   = "miss_flag",
  composite_value = "Comp",
  baseline_col    = "baseline"
)
#>   count imputed_value miss_flag baseline
#> 1     3             3      <NA>        4
#> 2    NA             7       MAR        4
#> 3    NA             4      Comp        4
#> 4     5             5      <NA>        6
```
