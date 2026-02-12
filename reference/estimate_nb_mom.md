# Method of Moments Estimation for Negative Binomial Parameters

Estimates the event rate(s) and common dispersion parameter (k) for
negative binomial count data using the method of moments. This is a
robust alternative to Maximum Likelihood Estimation (MLE), especially
when MLE fails to converge or produces boundary estimates.

## Usage

``` r
estimate_nb_mom(data, group = NULL)
```

## Arguments

- data:

  A data frame containing the data. Must include columns `events`
  (number of events) and `tte` (total exposure/follow-up time).

- group:

  Optional character string specifying the grouping column name (e.g.,
  "treatment"). If provided, rates are estimated separately for each
  group, while a common dispersion parameter is estimated across groups.
  If NULL (default), a single rate and dispersion are estimated (blinded
  case).

## Value

A list containing:

- lambda:

  Estimated event rate(s). A single numeric value if `group` is NULL, or
  a named vector if `group` is provided.

- dispersion:

  Estimated common dispersion parameter (k).

## Details

The method of moments estimator for the dispersion parameter \\k\\ is
derived by equating the theoretical variance to the observed second
central moment, accounting for varying exposure times.

For a given group with rate \\\lambda\\, the expected count for subject
\\i\\ is \\\mu_i = \lambda t_i\\. The variance is \\V_i = \mu_i + k
\mu_i^2\\. The estimator is calculated as: \$\$\hat{k} = \max\left(0,
\frac{\sum (y_i - \hat{\mu}\_i)^2 - \sum y_i}{\sum
\hat{\mu}\_i^2}\right)\$\$ where \\y_i\\ is the number of events,
\\t_i\\ is the exposure time, and \\\hat{\mu}\_i = \hat{\lambda} t_i\\
is the estimated expected count.

When multiple groups are present, the numerator and denominator are
summed across all groups to estimate a common \\k\\.

## Examples

``` r
# Blinded estimation (single group)
df <- data.frame(events = c(1, 2, 0, 3), tte = c(1, 1.2, 0.5, 1.5))
estimate_nb_mom(df)
#> $lambda
#> [1] 1.428571
#> 
#> $dispersion
#> [1] 0
#> 

# Unblinded estimation (two groups)
df_group <- df
df_group$group <- c("A", "A", "B", "B")
estimate_nb_mom(df_group, group = "group")
#> $lambda
#>        A        B 
#> 1.363636 1.500000 
#> 
#> $dispersion
#> [1] 0
#> 
```
