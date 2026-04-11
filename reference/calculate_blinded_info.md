# Calculate blinded statistical information

Estimates the blinded dispersion \\k\\ and event rate \\\lambda\\ from
pooled (blinded) interim data and calculates the observed statistical
information \\\mathcal{I}\\ for the log rate ratio. The estimation
follows the approach of Friede & Schmidli (2010): a single negative
binomial model is fit to the pooled data, then the estimated overall
rate is split into group-specific rates using the *planned* rate ratio.

## Usage

``` r
calculate_blinded_info(
  data,
  ratio = 1,
  lambda1_planning,
  lambda2_planning,
  event_gap = NULL
)
```

## Arguments

- data:

  A data frame containing the blinded interim data. Must include columns
  `events` (number of events) and `tte` (total exposure/follow-up time
  per subject).

- ratio:

  Planned allocation ratio \\r = n_2/n_1\\. Default is 1.

- lambda1_planning:

  Planned event rate \\\lambda_1\\ for the control group (used to
  determine the rate split).

- lambda2_planning:

  Planned event rate \\\lambda_2\\ for the experimental group (used to
  determine the rate split).

- event_gap:

  Optional gap duration (numeric). If provided, planning rates are
  adjusted to effective rates \\\lambda\_{\mathrm{eff}} = \lambda / (1 +
  \lambda \cdot \mathrm{gap})\\ before computing the rate split.

## Value

A list containing:

- blinded_info:

  Estimated statistical information \\\mathcal{I}\\.

- dispersion_blinded:

  Estimated dispersion parameter \\k\\.

- lambda_blinded:

  Estimated overall (pooled) event rate.

- lambda1_adjusted:

  Re-estimated control rate \\\hat\lambda_1\\.

- lambda2_adjusted:

  Re-estimated experimental rate \\\hat\lambda_2\\.

## Details

The statistical information is computed as: \$\$\mathcal{I} =
\frac{1}{1/W_1 + 1/W_2}\$\$ where \\W_g = p_g \sum_i \mu\_{g,i} / (1 +
k\\\mu\_{g,i})\\ and \\\mu\_{g,i} = \lambda_g t_i\\ is the expected
count for subject \\i\\ if they belonged to group \\g\\.

## References

Friede, T., & Schmidli, H. (2010). Blinded sample size reestimation with
negative binomial counts in superiority and non-inferiority trials.
*Methods of Information in Medicine*, 49(06), 618–624.
[doi:10.3414/ME09-02-0060](https://doi.org/10.3414/ME09-02-0060)

## See also

[`blinded_ssr()`](https://keaven.github.io/gsDesignNB/reference/blinded_ssr.md)
for blinded sample size reestimation;
[`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md)
for the underlying sample size formula.

## Examples

``` r
interim <- data.frame(events = c(1, 2, 1, 3), tte = c(0.8, 1.0, 1.2, 0.9))
calculate_blinded_info(
  interim,
  ratio = 1,
  lambda1_planning = 0.5,
  lambda2_planning = 0.3
)
#> $blinded_info
#> [1] 1.640582
#> 
#> $dispersion_blinded
#> [1] 1.617394e-05
#> 
#> $lambda_blinded
#> (Intercept) 
#>    1.794874 
#> 
#> $lambda1_adjusted
#> (Intercept) 
#>    2.243592 
#> 
#> $lambda2_adjusted
#> (Intercept) 
#>    1.346155 
#> 
```
