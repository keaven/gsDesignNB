# Seasonal event simulation

``` r
library(gsDesignNB)
library(data.table)
library(dplyr)
```

This vignette demonstrates how to simulate recurrent events when event
rates vary by season. This is useful for conditions with seasonal
patterns, such as respiratory infections.

## Simulation setup

We design a small trial with the following characteristics:

- **Randomization start:** January 1, 2024.
- **Sample size:** 20 subjects.
- **Enrollment:** 6 months duration.
- **Follow-up:** 1 year.
- **Seasonality:** Higher event rates in Winter/Fall, lower in
  Spring/Summer.

### Define parameters

``` r
# Randomization starts in Winter
rand_start <- as.Date("2024-01-01")

# Enrollment: 20 subjects over 6 months
enroll_rate <- data.frame(
  rate = 20 / 6,
  duration = 6
) # time unit: months?
# Note: In nb_sim_seasonal, we assumed rate units match max_followup units.
# Let's use YEARS as the time unit for consistency with the package conventions often used.
# If max_followup = 1 (year), then rates are per year.
# Enrollment duration = 0.5 years.

enroll_rate <- data.frame(
  rate = 20 / 0.5,
  duration = 0.5
)

# Seasonal Failure Rates (per year)
# Winter: High
# Spring: Low
# Summer: Low
# Fall: Medium-High

fail_rate <- data.frame(
  treatment = rep(c("Control", "Experimental"), each = 4),
  season = rep(c("Winter", "Spring", "Summer", "Fall"), 2),
  rate = c(
    # Control
    2.0, # Winter
    0.5, # Spring
    0.2, # Summer
    1.5, # Fall
    # Experimental (assume 30% reduction)
    2.0 * 0.7,
    0.5 * 0.7,
    0.2 * 0.7,
    1.5 * 0.7
  ),
  dispersion = 0.5 # Constant dispersion
)

# Dropout (5% per year)
dropout_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0.05, 0.05),
  duration = c(100, 100)
)
```

## Run simulation

We use
[`nb_sim_seasonal()`](https://keaven.github.io/gsDesignNB/reference/nb_sim_seasonal.md)
to simulate the trial.

``` r
set.seed(123)

sim_data <- nb_sim_seasonal(
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  dropout_rate = dropout_rate,
  max_followup = 1, # 1 year
  randomization_start_date = rand_start,
  n = 20
)

# View structure of the output
head(sim_data)
#>   id treatment enroll_time season     start       end event  cal_start
#> 1  1   Control  0.02108643 Winter 0.0000000 0.1431846     0 2024-01-08
#> 2  1   Control  0.02108643 Spring 0.1431846 0.3950669     0 2024-03-01
#> 3  1   Control  0.02108643 Summer 0.3950669 0.6469492     0 2024-06-01
#> 4  1   Control  0.02108643   Fall 0.6469492 0.8446997     1 2024-09-01
#> 5  1   Control  0.02108643   Fall 0.8446997 0.8960936     0 2024-11-12
#> 6  1   Control  0.02108643 Winter 0.8960936 1.0000000     0 2024-12-01
#>      cal_end
#> 1 2024-03-01
#> 2 2024-06-01
#> 3 2024-09-01
#> 4 2024-11-12
#> 5 2024-12-01
#> 6 2025-01-07
```

The output contains multiple rows per subject, split by season
intervals. `event` is 1 if an event occurred at the end of the interval,
and 0 otherwise.

## Analysis at a cut date

We can cut the data at a specific calendar date using
[`cut_data_by_date()`](https://keaven.github.io/gsDesignNB/reference/cut_data_by_date.md).
This function aggregates the seasonal intervals up to the cut date,
subtracting any event gaps if specified.

``` r
# Cut data at 9 months (0.75 years)
cut_date <- 0.75

# Use a small gap (e.g., 7 days)
gap_days <- 7 / 365.25

analysis_data <- cut_data_by_date(sim_data, cut_date = cut_date, event_gap = gap_days)

head(analysis_data)
#>   id    treatment enroll_time season events        tte  tte_total
#> 1  1      Control  0.02108643   Fall      0 0.08196441 0.08196441
#> 2  1      Control  0.02108643 Spring      0 0.25188227 0.25188227
#> 3  1      Control  0.02108643 Summer      0 0.25188227 0.25188227
#> 4  1      Control  0.02108643 Winter      0 0.14318462 0.14318462
#> 5  2 Experimental  0.03550169   Fall      0 0.08196441 0.08196441
#> 6  2 Experimental  0.03550169 Spring      0 0.25188227 0.25188227
```

The `analysis_data` is aggregated by subject and season. This allows for
seasonal adjustments in the analysis if desired.

``` r
# Summarize events by season
library(dplyr)
analysis_data %>%
  group_by(season, treatment) %>%
  summarize(
    Subjects = n_distinct(id),
    TotalExposure = sum(tte),
    TotalEvents = sum(events),
    Rate = sum(events) / sum(tte)
  )
#> `summarise()` has grouped output by 'season'. You can override using the
#> `.groups` argument.
#> # A tibble: 8 × 6
#> # Groups:   season [4]
#>   season treatment    Subjects TotalExposure TotalEvents  Rate
#>   <chr>  <chr>           <int>         <dbl>       <dbl> <dbl>
#> 1 Fall   Control            10         0.820           0 0    
#> 2 Fall   Experimental       10         0.781           2 2.56 
#> 3 Spring Control            10         2.02            1 0.495
#> 4 Spring Experimental       10         2.06            3 1.46 
#> 5 Summer Control            10         2.52            0 0    
#> 6 Summer Experimental       10         2.52            0 0    
#> 7 Winter Control             5         0.395           0 0    
#> 8 Winter Experimental        5         0.408           0 0
```

## Seasonal and treatment effect estimation

We can estimate the seasonal effects and the treatment effect using a
negative binomial generalized linear model. We use the logarithm of the
exposure time as an offset.

``` r
# Fit negative binomial GLM
# We include season and treatment in the model
fit <- MASS::glm.nb(
  events ~ treatment + season + offset(log(tte)),
  data = analysis_data[analysis_data$tte > 0, ]
)

# Summary of the model
summary(fit)
#> 
#> Call:
#> MASS::glm.nb(formula = events ~ treatment + season + offset(log(tte)), 
#>     data = analysis_data[analysis_data$tte > 0, ], init.theta = 0.7046293947, 
#>     link = log)
#> 
#> Coefficients:
#>                         Estimate Std. Error z value Pr(>|z|)
#> (Intercept)           -9.344e-01  1.230e+00  -0.760    0.447
#> treatmentExperimental  1.712e+00  1.172e+00   1.460    0.144
#> seasonSpring          -1.881e-01  9.655e-01  -0.195    0.845
#> seasonSummer          -3.559e+01  1.501e+07   0.000    1.000
#> seasonWinter          -3.409e+01  1.850e+07   0.000    1.000
#> 
#> (Dispersion parameter for Negative Binomial(0.7046) family taken to be 1)
#> 
#>     Null deviance: 28.268  on 69  degrees of freedom
#> Residual deviance: 17.594  on 65  degrees of freedom
#> AIC: 45.764
#> 
#> Number of Fisher Scoring iterations: 1
#> 
#> 
#>               Theta:  0.70 
#>           Std. Err.:  1.24 
#> 
#>  2 x log-likelihood:  -33.764

# Extract estimates
coef_summary <- summary(fit)$coefficients

# Seasonal effects (relative to reference season, likely Fall based on alphabetical order)
# Treatment effect (Experimental vs Control)
print(coef_summary)
#>                          Estimate   Std. Error       z value  Pr(>|z|)
#> (Intercept)            -0.9344423 1.229766e+00 -7.598537e-01 0.4473421
#> treatmentExperimental   1.7117150 1.172198e+00  1.460261e+00 0.1442185
#> seasonSpring           -0.1881428 9.654901e-01 -1.948676e-01 0.8454966
#> seasonSummer          -35.5862751 1.500600e+07 -2.371470e-06 0.9999981
#> seasonWinter          -34.0933221 1.850416e+07 -1.842468e-06 0.9999985

# Estimated Rate Ratio (Experimental / Control)
rr <- exp(coef(fit)["treatmentExperimental"])
message("Estimated Rate Ratio (Experimental / Control): ", round(rr, 3))
#> Estimated Rate Ratio (Experimental / Control): 5.538
message("True Design Rate Ratio: ", 0.7)
#> True Design Rate Ratio: 0.7
```

### Variance and information

The variance of the treatment effect estimate can be extracted from the
covariance matrix of the model. The statistical information is the
inverse of this variance.

``` r
# Variance of the treatment effect coefficient
var_beta <- vcov(fit)["treatmentExperimental", "treatmentExperimental"]
message("Variance of Treatment Effect (log-scale): ", var_beta)
#> Variance of Treatment Effect (log-scale): 1.37404900110776

# Statistical Information
info <- 1 / var_beta
message("Statistical Information: ", info)
#> Statistical Information: 0.727776083090051
```
