# Multiple imputation for longitudinal negative binomial counts

## Introduction

Longitudinal clinical trials with recurrent-event count endpoints (e.g.,
MRI lesion counts in multiple sclerosis, exacerbation counts in COPD)
routinely encounter missing data due to intercurrent events (ICEs) such
as premature treatment discontinuation, adverse events requiring rescue
therapy, or death. Proper handling of these missing values requires a
principled strategy that reflects the scientific question—the
*estimand*—as defined by the ICH E9(R1) addendum.

The `gsDesignNB` package provides three imputation strategies for
longitudinal negative binomial count data:

| Mechanism                                 | Flag value | Strategy                                    |
|:------------------------------------------|:-----------|:--------------------------------------------|
| Missing at Random                         | `"MAR"`    | GLMM-based imputation with subject BLUPs    |
| Missing Not at Random (non-reference arm) | `"MNAR"`   | Reference-based (copy-reference) imputation |
| Intercurrent event — composite            | `"Comp"`   | Baseline count carried forward              |

These mirror the SAS implementation using `PROC GLIMMIX` with
`dist=negbin link=log` and `PROC PLM` for counterfactual predictions.

## Statistical background

### Negative binomial count model

For subject i at visit t, let Y\_{it} denote the count endpoint. The
negative binomial (NB) model used throughout **gsDesignNB** is the
Gamma–Poisson mixture:

\Lambda_i \sim \text{Gamma}\\\left(\frac{1}{k},\\ \mu_i k\right), \qquad
Y\_{it} \mid \Lambda_i \sim \text{Poisson}(\Lambda_i).

Marginally, Y\_{it} has mean \mu_i and variance \text{Var}(Y\_{it}) =
\mu_i + k\\\mu_i^2, where k \geq 0 is the dispersion parameter. In the
longitudinal GLMM, the mean is modelled as

\log(\mu\_{it}) = \mathbf{x}\_{it}^{\top}\boldsymbol{\beta} + b_i,

where \mathbf{x}\_{it} contains fixed-effect covariates (treatment,
visit, baseline, stratification factors) and b_i is a subject-level
random intercept. The model is fitted via \[glmmTMB::glmmTMB()\] with
`family = nbinom2(link = "log")`.

### Imputation draw

All three model-based strategies draw imputed counts from the same
Gamma–Poisson compound:

\lambda_i^{(m)} \sim \text{Gamma}\\\left(\frac{1}{k},\\ \hat\mu_i \cdot
k\right), \qquad Y_i^{(m)} \mid \lambda_i^{(m)} \sim
\text{Poisson}\\\left(\lambda_i^{(m)}\right).

The strategies differ only in how \hat\mu_i is determined:

- **MAR**: \hat\mu_i = \hat\mu_i^{\text{BLUP}} — predicted mean
  including subject random effects.
- **MNAR reference-based**: \hat\mu_i = \hat\mu_i^{\text{FE, ref}}
  \times \dfrac{\hat\mu_i^{\text{BLUP}}}{\hat\mu_i^{\text{FE}}} —
  counterfactual mean under the reference treatment, adjusted by each
  subject’s estimated random effect. The ratio \hat\mu_i^{\text{BLUP}} /
  \hat\mu_i^{\text{FE}} is the multiplicative random-effect “BLUP
  multiplier” computed from the original treatment assignment.
- **Composite**: Y_i^{(m)} = y\_{i,0} (baseline value), no model draw.

### Bootstrap–MI combination

When `n_boot > 1`, subjects are resampled with replacement within strata
before each GLMM fit. This *boot-MI* approach propagates both imputation
uncertainty and model-fitting uncertainty into the final variance
estimate without requiring Rubin’s combining rules. When `n_boot = 1`,
standard MI is performed and Rubin’s rules should be applied when
pooling estimates across the `n_imp` imputed datasets.

## Worked example

### Simulating longitudinal count data

We simulate a 3-arm trial (placebo and two active doses) with 4
post-baseline visits, negative binomial counts (dispersion k = 0.5), and
a single binary stratification factor.

``` r
set.seed(2025)

n_subj   <- 90L   # 30 per arm
n_visit  <- 4L
k_true   <- 0.5   # dispersion
lambda   <- c(placebo = 1.8, low = 1.1, high = 0.7)   # mean counts/visit

trt_labels <- rep(c("placebo", "low", "high"), each = n_subj / 3L)
strat      <- sample(c(0L, 1L), n_subj, replace = TRUE)
baseline   <- rnbinom(n_subj, mu = 2.5, size = 1 / k_true)
id         <- seq_len(n_subj)

# Build long-format data (all visits complete at first)
long_data <- do.call(rbind, lapply(seq_len(n_subj), function(i) {
  mu_i <- lambda[trt_labels[i]] * exp(0.15 * strat[i] + 0.05 * baseline[i])
  data.frame(
    id       = i,
    visit    = seq_len(n_visit),
    trt      = trt_labels[i],
    strat1   = strat[i],
    baseline = baseline[i],
    count    = rnbinom(n_visit, mu = mu_i, size = 1 / k_true),
    stringsAsFactors = FALSE
  )
}))

# Sort
long_data <- long_data[order(long_data$id, long_data$visit), ]
rownames(long_data) <- NULL
str(long_data)
#> 'data.frame':    360 obs. of  6 variables:
#>  $ id      : int  1 1 1 1 2 2 2 2 3 3 ...
#>  $ visit   : int  1 2 3 4 1 2 3 4 1 2 ...
#>  $ trt     : chr  "placebo" "placebo" "placebo" "placebo" ...
#>  $ strat1  : int  0 0 0 0 1 1 1 1 1 1 ...
#>  $ baseline: num  3 3 3 3 5 5 5 5 4 4 ...
#>  $ count   : num  1 2 3 2 1 6 4 4 3 1 ...
```

### Introducing missing data by mechanism

We introduce three types of missing data to mirror a realistic trial:

- **MAR**: A random subset of subjects miss visits 3–4 due to
  non-disease reasons (e.g., logistical dropout). The probability of
  dropout is unrelated to unobserved counts given the observed history —
  a plausible MAR assumption.
- **MNAR**: A subset of active-arm subjects discontinue at visit 2 due
  to adverse events likely linked to drug exposure. Their later counts
  are expected to be higher (on the reference arm trajectory) than for
  subjects who remained on treatment — a MNAR mechanism for which
  reference-based imputation is appropriate.
- **Composite ICE**: A small number of subjects in any arm experience a
  severe disease worsening event (e.g., hospitalisation) that leads to
  treatment discontinuation. The composite estimand treats this event as
  part of the outcome; missing post-event counts are filled with the
  baseline count.

``` r
long_data$miss_flag <- NA_character_

# --- MAR: ~15% of subjects drop out from visit 3 onward ---
mar_subjects <- sample(id, size = round(0.15 * n_subj))
long_data$miss_flag[
  long_data$id %in% mar_subjects & long_data$visit >= 3
] <- "MAR"

# --- MNAR: ~10% of active-arm subjects withdraw at visit 2 ---
active_ids  <- id[trt_labels != "placebo"]
mnar_subjects <- sample(active_ids, size = round(0.10 * n_subj))
long_data$miss_flag[
  long_data$id %in% mnar_subjects & long_data$visit >= 2
] <- "MNAR"

# --- Composite ICE: ~5% of all subjects, post-event visits ---
comp_subjects <- sample(setdiff(id, c(mar_subjects, mnar_subjects)),
                        size = round(0.05 * n_subj))
comp_ice_visit <- sample(2L:4L, length(comp_subjects), replace = TRUE)
for (ci in seq_along(comp_subjects)) {
  long_data$miss_flag[
    long_data$id == comp_subjects[ci] &
    long_data$visit >= comp_ice_visit[ci]
  ] <- "Comp"
}

# Set the outcome to NA for flagged rows (simulate non-collection)
long_data$count[!is.na(long_data$miss_flag)] <- NA_integer_

cat("Missing rows by mechanism:\n")
#> Missing rows by mechanism:
print(table(long_data$miss_flag, useNA = "ifany"))
#> 
#> Comp  MAR MNAR <NA> 
#>    8   26   27  299
```

### Running `impute_nb()`

We now run the full MI pipeline. For this illustration we use
`n_boot = 1` (no bootstrap) and `n_imp = 5` imputations, fitting a
random-intercept NB GLMM.

``` r
library(gsDesignNB)

# The formula mirrors the fixed-effects structure from the SAS code.
# A random intercept per subject captures between-subject heterogeneity.
mi_formula <- count ~ baseline + strat1 + trt + visit + (1 | id)

imp_result <- impute_nb(
  data            = long_data,
  formula         = mi_formula,
  outcome_col     = "count",
  miss_flag_col   = "miss_flag",
  baseline_col    = "baseline",
  trt_col         = "trt",
  reference_trt   = "placebo",
  subject_col     = "id",
  strata_cols     = c("trt", "strat1"),
  mar_values      = "MAR",
  mnar_value      = "MNAR",
  composite_value = "Comp",
  n_imp           = 5L,
  n_boot          = 1L,
  seed            = 42L
)

# Dimensions: n_subj * n_visit * n_imp rows
dim(imp_result)
#> [1] 1800   10
names(imp_result)
#>  [1] "id"            "visit"         "trt"           "strat1"       
#>  [5] "baseline"      "count"         "miss_flag"     "replicate"    
#>  [9] "imputation"    "imputed_value"
```

### Inspecting imputed values

The `imputed_value` column contains:

- The **observed value** for non-missing rows.
- An **imputed draw** for missing rows.

``` r
# Show imputed rows only, first 2 imputations
imp_missing <- imp_result[
  !is.na(imp_result$miss_flag) & imp_result$imputation <= 2,
  c("id", "visit", "trt", "miss_flag", "baseline", "count", "imputation",
    "imputed_value")
]
head(imp_missing[order(imp_missing$miss_flag, imp_missing$id,
                        imp_missing$imputation), ], 20)
#>     id visit     trt miss_flag baseline count imputation imputed_value
#> 2    1     2 placebo      Comp        3    NA          1             3
#> 3    1     3 placebo      Comp        3    NA          1             3
#> 4    1     4 placebo      Comp        3    NA          1             3
#> 362  1     2 placebo      Comp        3    NA          2             3
#> 363  1     3 placebo      Comp        3    NA          2             3
#> 364  1     4 placebo      Comp        3    NA          2             3
#> 120 30     4 placebo      Comp        7    NA          1             7
#> 480 30     4 placebo      Comp        7    NA          2             7
#> 216 54     4     low      Comp        2    NA          1             2
#> 576 54     4     low      Comp        2    NA          2             2
#> 354 89     2    high      Comp        0    NA          1             0
#> 355 89     3    high      Comp        0    NA          1             0
#> 356 89     4    high      Comp        0    NA          1             0
#> 714 89     2    high      Comp        0    NA          2             0
#> 715 89     3    high      Comp        0    NA          2             0
#> 716 89     4    high      Comp        0    NA          2             0
#> 15   4     3 placebo       MAR        0    NA          1            13
#> 16   4     4 placebo       MAR        0    NA          1             7
#> 375  4     3 placebo       MAR        0    NA          2             0
#> 376  4     4 placebo       MAR        0    NA          2             0
```

### Comparing imputed means by strategy

We can verify that the imputation behaved as expected. Under MAR,
imputed values should centre around the model prediction for each arm.
Under reference-based MNAR, imputed values for non-placebo subjects
should look more like placebo outcomes (higher counts).

``` r
# Mean imputed value per mechanism and treatment arm
agg <- aggregate(
  imputed_value ~ miss_flag + trt,
  data   = imp_result[!is.na(imp_result$miss_flag), ],
  FUN    = mean
)
agg <- agg[order(agg$miss_flag, agg$trt), ]
print(agg, row.names = FALSE)
#>  miss_flag     trt imputed_value
#>       Comp    high      0.000000
#>       Comp     low      2.000000
#>       Comp placebo      4.000000
#>        MAR    high      0.600000
#>        MAR     low      1.833333
#>        MAR placebo      1.628571
#>       MNAR    high      1.613333
#>       MNAR     low      3.050000
```

The MNAR reference-based imputed means for the `low` and `high` arms
should be closer to the `placebo` mean than the MAR-imputed means,
reflecting the “what if this subject had been on placebo?”
counterfactual.

### Pooled analysis with Rubin’s rules (standard MI)

After standard MI (`n_boot = 1`), pool estimates across the `n_imp`
imputed datasets using Rubin’s rules. Here we compute the mean total
count per arm at visit 4 as a simple illustration.

``` r
# Per-imputation mean at visit 4
v4 <- imp_result[imp_result$visit == 4, ]
per_imp <- lapply(split(v4, v4$imputation), function(d) {
  tapply(d$imputed_value, d$trt, mean, na.rm = TRUE)
})

# Point estimate: average over imputations
Q_bar <- Reduce("+", per_imp) / length(per_imp)

# Within-imputation variance (using single-imputation SE ~ mean/n, illustrative)
# In practice use regression SE from each imputed dataset
cat("Pooled mean imputed count at visit 4 by treatment arm:\n")
#> Pooled mean imputed count at visit 4 by treatment arm:
print(round(Q_bar, 2))
#>    high     low placebo 
#>    1.05    1.46    2.09
```

### Bootstrap–MI for variance estimation

When `n_boot > 1`, no Rubin’s rules are needed. The variance across all
`n_boot * n_imp` completed datasets directly incorporates both sources
of uncertainty.

``` r
imp_boot <- impute_nb(
  data            = long_data,
  formula         = mi_formula,
  outcome_col     = "count",
  miss_flag_col   = "miss_flag",
  baseline_col    = "baseline",
  trt_col         = "trt",
  reference_trt   = "placebo",
  subject_col     = "id",
  strata_cols     = c("trt", "strat1"),
  n_imp           = 5L,
  n_boot          = 3L,   # use larger values in practice (e.g., 100–500)
  seed            = 99L
)

# Each replicate × imputation combination is an independent completed dataset
cat("Total rows:", nrow(imp_boot),
    "= n_subj * n_visit * n_boot * n_imp =",
    n_subj * n_visit * 3L * 5L, "\n")
#> Total rows: 5400 = n_subj * n_visit * n_boot * n_imp = 5400

# Variance of the arm-level mean at visit 4 across all completed datasets
v4b <- imp_boot[imp_boot$visit == 4, ]
dataset_means <- tapply(
  v4b$imputed_value,
  list(paste(v4b$replicate, v4b$imputation, sep = "_"), v4b$trt),
  mean
)
cat("\nVariance of dataset-level means across boot×imp datasets:\n")
#> 
#> Variance of dataset-level means across boot×imp datasets:
print(round(apply(dataset_means, 2, var), 4))
#>    high     low placebo 
#>  0.1420  0.1384  0.1264
```

### Using the building-block functions directly

Advanced users can call the individual strategy functions for more
control.

#### Step 1 — Fit the GLMM

``` r
obs_data <- long_data[!is.na(long_data$count), ]
fits <- fit_nb_glmm(
  data    = obs_data,
  formula = count ~ baseline + strat1 + trt + visit + (1 | id)
)
cat("Estimated dispersion k:", round(fits[["1"]]$k, 3), "\n")
#> Estimated dispersion k: 2.648
```

#### Step 2 — MAR imputation

``` r
imp_mar <- impute_nb_mar(
  data          = long_data,
  fits          = fits,
  outcome_col   = "count",
  miss_flag_col = "miss_flag",
  mar_values    = c("MAR", "MNAR"),  # reference-arm MNAR treated as MAR
  n_imp         = 3L
)
cat("MAR imputed rows:", sum(!is.na(imp_mar$imputed_value[
  imp_mar$imputation == 1 & !is.na(imp_mar$miss_flag)
])), "\n")
#> MAR imputed rows: 53
```

#### Step 3 — MNAR reference-based imputation

``` r
imp_mnar <- impute_nb_mnar_ref(
  data          = long_data,
  fits          = fits,
  outcome_col   = "count",
  miss_flag_col = "miss_flag",
  mnar_value    = "MNAR",
  trt_col       = "trt",
  reference_trt = "placebo",
  n_imp         = 3L
)
```

#### Step 4 — Composite strategy (no model required)

``` r
library(gsDesignNB)

# Composite fill can be applied to any data frame, no glmmTMB needed
df_comp_example <- data.frame(
  count         = c(3L,   NA,     NA,      5L),
  imputed_value = c(3L,    7L,    NA,      5L),
  miss_flag     = c(NA,  "MAR", "Comp",    NA),
  baseline      = c(4L,    4L,    4L,      6L)
)
impute_nb_composite(
  df_comp_example,
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

The third row (`miss_flag = "Comp"`) now has `imputed_value = 4` (the
baseline value); the second row (`miss_flag = "MAR"`) retains its
previously imputed value of 7.

## Key considerations

**Choice of random-effect structure.** The original SAS model uses an
unstructured residual covariance across visits within each subject ×
endpoint combination. In `glmmTMB` this can be approximated with
`(0 + visit_factor | id:param)`, but convergence may be slow for small
datasets. Start with a random intercept `(1 | id)` and add complexity
incrementally.

**Treatment variable type.** If `trt_col` is a factor, ensure that
`reference_trt` matches a level in the factor. For counterfactual
prediction the reference level is inserted into the newdata frame; if
new levels arise, `allow.new.levels = TRUE` handles them gracefully.

**Composite and MNAR rows that are observed.** The composite strategy
only overwrites rows where `is.na(outcome_col)`. Observed rows with
`miss_flag = "Comp"` (e.g., the actual ICE visit itself, if observed)
are left unchanged.

**Bootstrap sample size.** For the boot-MI approach,
`n_boot = 100`–`500` replicates is typical in practice. The smaller
values used in this vignette are for illustration only.

**Pooling.** With standard MI (`n_boot = 1`), pool regression-model
estimates across imputations using Rubin’s rules. With boot-MI
(`n_boot > 1`), the empirical variance across all `n_boot × n_imp`
datasets is a valid variance estimator without additional combining
rules.

## Session info

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] gsDesignNB_0.3.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6        TMB_1.9.21          xfun_0.57          
#>  [4] bslib_0.10.0        ggplot2_4.0.2       htmlwidgets_1.6.4  
#>  [7] lattice_0.22-9      numDeriv_2016.8-1.1 vctrs_0.7.3        
#> [10] tools_4.5.3         Rdpack_2.6.6        generics_0.1.4     
#> [13] parallel_4.5.3      sandwich_3.1-1      tibble_3.3.1       
#> [16] pkgconfig_2.0.3     Matrix_1.7-4        data.table_1.18.2.1
#> [19] RColorBrewer_1.1-3  S7_0.2.1-1          desc_1.4.3         
#> [22] gt_1.3.0            lifecycle_1.0.5     doFuture_1.2.1     
#> [25] compiler_4.5.3      farver_2.1.2        textshaping_1.0.5  
#> [28] codetools_0.2-20    htmltools_0.5.9     sass_0.4.10        
#> [31] yaml_2.3.12         pkgdown_2.2.0       pillar_1.11.1      
#> [34] nloptr_2.2.1        gsDesign_3.9.0      jquerylib_0.1.4    
#> [37] tidyr_1.3.2         MASS_7.3-65         cachem_1.1.0       
#> [40] iterators_1.0.14    reformulas_0.4.4    foreach_1.5.2      
#> [43] boot_1.3-32         parallelly_1.47.0   nlme_3.1-168       
#> [46] r2rtf_1.3.0         tidyselect_1.2.1    digest_0.6.39      
#> [49] mvtnorm_1.3-7       future_1.70.0       listenv_0.10.1     
#> [52] dplyr_1.2.1         purrr_1.2.2         splines_4.5.3      
#> [55] fastmap_1.2.0       grid_4.5.3          cli_3.6.6          
#> [58] magrittr_2.0.5      survival_3.8-6      future.apply_1.20.2
#> [61] scales_1.4.0        simtrial_1.0.2      rmarkdown_2.31     
#> [64] globals_0.19.1      otel_0.2.0          lme4_2.0-1         
#> [67] ragg_1.5.2          zoo_1.8-15          evaluate_1.0.5     
#> [70] knitr_1.51          glmmTMB_1.1.14      rbibutils_2.4.1    
#> [73] mgcv_1.9-4          rlang_1.2.0         Rcpp_1.1.1-1       
#> [76] xtable_1.8-8        glue_1.8.1          xml2_1.5.2         
#> [79] minqa_1.2.8         jsonlite_2.0.0      R6_2.6.1           
#> [82] systemfonts_1.3.2   fs_2.1.0
```
