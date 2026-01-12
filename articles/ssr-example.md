# Sample size re-estimation example

``` r
library(gsDesignNB)
library(gsDesign)
library(data.table)
library(MASS)
library(gt)
library(gt)
```

## Introduction

This vignette demonstrates a group sequential trial with a single
interim analysis where unblinded sample size re-estimation (SSR) is
performed.

We simulate a scenario where the **control event rate is lower than
assumed** and the **dispersion is higher than assumed**. Both factors
lead to slower information accumulation (or higher variance) than
planned. We will show how to:

1.  Monitor information accumulation at an interim analysis.
2.  Re-estimate the required sample size (or duration) to achieve the
    target power.
3.  Adjust the final analysis bounds if the final information differs
    from the plan.

## Trial setup and initial design

**Planned parameters:**

- Control rate (\\\lambda_1\\): 0.1 events/month
- Experimental rate (\\\lambda_2\\): 0.075 events/month (Hazard Ratio =
  0.75)
- Dispersion (\\k\\): 0.5
- Power: 90%
- One-sided Type I error (\\\alpha\\): 0.025
- Enrollment: 20 patients/month for 12 months (Total N = 240)
- Study duration: 24 months

**Actual parameters (simulation truth):**

- Control rate (\\\lambda_1\\): 0.08 events/month (Lower than planned)
- Experimental rate (\\\lambda_2\\): 0.06 events/month (HR = 0.75
  maintained)
- Dispersion (\\k\\): 0.65 (Higher than planned)

### Initial sample size calculation

``` r
# Planned parameters
lambda1_plan <- 0.1
lambda2_plan <- 0.075
k_plan <- 0.5
power_plan <- 0.9
alpha_plan <- 0.025
accrual_rate_plan <- 20
accrual_dur_plan <- 12
trial_dur_plan <- 24

# Calculate sample size
design_plan <- sample_size_nbinom(
  lambda1 = lambda1_plan,
  lambda2 = lambda2_plan,
  dispersion = k_plan,
  power = power_plan,
  alpha = alpha_plan,
  accrual_rate = accrual_rate_plan,
  accrual_duration = accrual_dur_plan,
  trial_duration = trial_dur_plan,
  max_followup = 12
)

# Convert to group sequential design
gs_plan <- design_plan |> 
   gsNBCalendar(
     k = 2,
     test.type = 4, # Non-binding futility
     alpha = alpha_plan,
     sfu = sfHSD,
     sfupar = -2, # Moderately aggressive alpha-spending
     sfl = sfHSD,
     sflpar = 1, # Pocock-like futility spending
     analysis_times = c(accrual_dur_plan - 2, trial_dur_plan)
   ) |> gsDesignNB::toInteger() # Round to integer sample size
```

``` r
summary(gs_plan)
#> Asymmetric two-sided with non-binding futility bound group sequential design
#> for negative binomial outcomes, 2 analyses, total sample size 882.0, 90 percent
#> power, 2.5 percent (1-sided) Type I error. Control rate 0.1000, treatment rate
#> 0.0750, risk ratio 0.7500, dispersion 0.5000. Accrual duration 12.0, trial
#> duration 24.0, max follow-up 12.0, average exposure 12.00. Randomization ratio
#> 1:1. Upper spending: Hwang-Shih-DeCani (gamma = -2) Lower spending:
#> Hwang-Shih-DeCani (gamma = 1)
#> Asymmetric two-sided with non-binding futility bound group sequential design
#> for negative binomial outcomes, 2 analyses, total sample size 882.0, 90 percent
#> power, 2.5 percent (1-sided) Type I error. Control rate 0.1000, treatment rate
#> 0.0750, risk ratio 0.7500, dispersion 0.5000. Accrual duration 12.0, trial
#> duration 24.0, max follow-up 12.0, average exposure 12.00. Randomization ratio
#> 1:1. Upper spending: Hwang-Shih-DeCani (gamma = -2) Lower spending:
#> Hwang-Shih-DeCani (gamma = 1)
```

``` r
gsBoundSummary(gs_plan,
    deltaname = "RR",
    logdelta = TRUE,
    Nname = "Information",
    timename = "Month",
    digits = 4,
    ddigits = 2) |> gt() |>
  tab_header(
    title = "Group Sequential Design Bounds for Negative Binomial Outcome",
    subtitle = paste0(
      "N = ", ceiling(gs_plan$n_total[gs_plan$k]),
      ", Expected events = ", round(gs_plan$nb_design$total_events, 1)
    )
  )
```

| Group Sequential Design Bounds for Negative Binomial Outcome |                     |          |          |
|--------------------------------------------------------------|---------------------|----------|----------|
| N = 882, Expected events = 785.4                             |                     |          |          |
| Analysis                                                     | Value               | Efficacy | Futility |
| IA 1: 43%                                                    | Z                   | 2.5498   | 0.7234   |
| Information: 64.85                                           | p (1-sided)         | 0.0054   | 0.2347   |
| Month: 10                                                    | ~RR at bound        | 0.7286   | 0.9141   |
|                                                              | P(Cross) if RR=1    | 0.0054   | 0.7653   |
|                                                              | P(Cross) if RR=0.75 | 0.4077   | 0.0556   |
| Final                                                        | Z                   | 2.0152   | 2.0152   |
| Information: 149.77                                          | p (1-sided)         | 0.0219   | 0.0219   |
| Month: 24                                                    | ~RR at bound        | 0.8481   | 0.8481   |
|                                                              | P(Cross) if RR=1    | 0.0219   | 0.9781   |
|                                                              | P(Cross) if RR=0.75 | 0.9027   | 0.0973   |

## Simulation

We simulate the trial with the same rate ratio, but lower actual rates
and a higher dispersion. We will limit the maximum sample size to 40%
more than the planned size (i.e., max N = ) to avoid unbounded
increases.

``` r
set.seed(1234)

# Actual parameters
lambda1_true <- 0.08 # Lower event rates, same rr
lambda2_true <- 0.06
k_true <- 0.65 # Higher dispersion

# Enrollment and rates for simulation
# We simulate a larger pool to allow for potential sample size increase
enroll_rate <- data.frame(rate = accrual_rate_plan, duration = accrual_dur_plan * 2) 
fail_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(lambda1_true, lambda2_true),
  dispersion = k_true
)
dropout_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0, 0),
  duration = c(100, 100)
)

sim_data <- nb_sim(
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  dropout_rate = dropout_rate,
  max_followup = trial_dur_plan, 
  n = 600 
)

# Limit to planned enrollment for the initial cut
# We will "open" more enrollment if needed
sim_data_planned <- sim_data[1:ceiling(design_plan$n_total), ]
```

## Interim analysis

We perform an interim analysis at **Month 10**, which is 2 months prior
to the end of planned enrollment (Month 12).

``` r
interim_time <- 10
interim_data <- cut_data_by_date(sim_data_planned, cut_date = interim_time)

# Summary
table(interim_data$treatment)
#> 
#>      Control Experimental 
#>           99           98
sum(interim_data$events)
#> [1] 69
mean(interim_data$tte)
#> [1] 5.086404
```

### Information computation

We calculate the statistical information accumulated so far.

**Blinded Information:** Calculated assuming a common rate and
dispersion, often used for monitoring without unblinding.

**Unblinded Information:** Calculated using the observed rates and
dispersion in each group. This is the true Fisher information for the
log rate ratio.

\\ \mathcal{I} = \left( \frac{1}{\text{Var}(\hat{\beta}\_{trt})} \right)
\\

``` r
# Blinded SSR
blinded_res <- blinded_ssr(
  data = interim_data,
  ratio = 1,
  lambda1_planning = lambda1_plan,
  lambda2_planning = lambda2_plan,
  power = power_plan,
  alpha = alpha_plan,
  accrual_rate = accrual_rate_plan,
  accrual_duration = accrual_dur_plan,
  trial_duration = trial_dur_plan
)

# Unblinded SSR
ssr_res <- unblinded_ssr(
  data = interim_data,
  ratio = 1,
  lambda1_planning = lambda1_plan,
  lambda2_planning = lambda2_plan,
  power = power_plan,
  alpha = alpha_plan,
  accrual_rate = accrual_rate_plan,
  accrual_duration = accrual_dur_plan,
  trial_duration = trial_dur_plan
)

cat("Blinded SSR N:", ceiling(blinded_res$n_total_blinded), "\n")
#> Blinded SSR N: 914
cat("Unblinded SSR N:", ceiling(ssr_res$n_total_unblinded), "\n")
#> Unblinded SSR N: 888
print(ssr_res)
#> $n_total_unblinded
#> (Intercept) 
#>         888 
#> 
#> $dispersion_unblinded
#> [1] 0.8891071
#> 
#> $lambda1_unblinded
#> (Intercept) 
#>  0.07878905 
#> 
#> $lambda2_unblinded
#> (Intercept) 
#>  0.05726336 
#> 
#> $info_fraction
#> [1] 0.0950141
#> 
#> $unblinded_info
#> [1] 12.06309
#> 
#> $target_info
#> [1] 126.9611
```

The estimated control rate is 0.079, which is lower than the planned
0.1. The information fraction is 0.095.

### Sample size re-estimation

Since the control rate is lower and dispersion is higher, the
information accumulates slower than planned. To maintain power, we need
to increase the sample size or follow-up.

The `unblinded_ssr` function estimates the required total sample size to
be 888. We ensure the sample size is at least as large as the planned
sample size.

``` r
n_planned <- ceiling(design_plan$n_total)
n_estimated <- ceiling(ssr_res$n_total_unblinded)

# Ensure we don't decrease sample size
n_final <- max(n_planned, n_estimated)

if (n_final > n_planned) {
  cat("Increasing sample size from", n_planned, "to", n_final, "\n")
  
  # Calculate new durations based on constant accrual rate
  # We extend enrollment to reach the new target N
  new_accrual_dur <- n_final / accrual_rate_plan
  
  # We maintain the planned max follow-up of 12 months
  # So the trial duration extends by the same amount as the accrual duration
  new_trial_dur <- new_accrual_dur + 12 
  
  cat("New Accrual Duration:", round(new_accrual_dur, 2), "months\n")
  cat("New Trial Duration:", round(new_trial_dur, 2), "months\n")
} else {
  cat("No sample size increase needed.\n")
  new_accrual_dur <- accrual_dur_plan
  new_trial_dur <- trial_dur_plan
  n_final <- n_planned
}
#> Increasing sample size from 748 to 888 
#> New Accrual Duration: 44.4 months
#> New Trial Duration: 56.4 months
```

## Final analysis

We proceed to the final analysis at Month 56.4. We include the
additional patients if the sample size was increased.

``` r
# Update the dataset to include the additional patients
sim_data_final <- sim_data[1:n_final, ]

# Cut at final analysis time
final_data <- cut_data_by_date(sim_data_final, cut_date = new_trial_dur)

# Average exposure
mean(final_data$tte)
#> [1] 23.95368

# Fit final model
fit_final <- glm.nb(events ~ treatment + offset(log(tte)), data = final_data)
summary(fit_final)
#> 
#> Call:
#> glm.nb(formula = events ~ treatment + offset(log(tte)), data = final_data, 
#>     init.theta = 1.782870921, link = log)
#> 
#> Coefficients:
#>                       Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)           -2.58248    0.08222 -31.411   <2e-16 ***
#> treatmentExperimental -0.14688    0.11826  -1.242    0.214    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for Negative Binomial(1.7829) family taken to be 1)
#> 
#>     Null deviance: 364.75  on 330  degrees of freedom
#> Residual deviance: 363.20  on 329  degrees of freedom
#> AIC: 1166
#> 
#> Number of Fisher Scoring iterations: 1
#> 
#> 
#>               Theta:  1.783 
#>           Std. Err.:  0.316 
#> 
#>  2 x log-likelihood:  -1159.959

# Calculate final information
var_beta <- vcov(fit_final)[2, 2]
final_info <- 1 / var_beta
cat("Final Information:", final_info, "\n")
#> Final Information: 71.50829
cat("Target Information:", ssr_res$target_info, "\n")
#> Target Information: 126.9611
```

### Updating bounds

If the final information differs from the target (e.g., it is lower), we
need to adjust the critical values to ensure we spend exactly
\\\alpha\\. We use the `gsDesign` package with `usTime` (user-specified
information fraction) to update the bounds.

``` r
# Information fraction at final analysis
# If final_info < target_info, fraction < 1. 
# However, for the final analysis, we typically want to spend all remaining alpha.
# We can treat the current info as the "maximum" info for the spending function.

# Example group sequential design (O'Brien-Fleming spending)
# We assume the interim was at the observed information fraction
info_fractions <- c(ssr_res$info_fraction, final_info / ssr_res$target_info)

# If final fraction < 1, we can force it to 1 to spend all alpha, 
# effectively accepting the lower power.
# Or we can use the actual information fractions in a spending function approach.

# Let's assume we want to spend all alpha at the final analysis regardless of info.
# We set the final timing to 1.
info_fractions_spending <- c(ssr_res$info_fraction, 1)

# Update bounds
gs_update <- gsDesign(
  k = 2,
  test.type = 1,
  alpha = alpha_plan,
  sfu = sfLDOF,
  usTime = info_fractions_spending,
  n.I = c(ssr_res$unblinded_info, final_info)
)

gsBoundSummary(gs_update,
    deltaname = "RR",
    logdelta = TRUE,
    Nname = "Information",
    timename = "Month",
    digits = 4,
    ddigits = 2) |> gt() |>
  tab_header(
    title = "Updated Group Sequential Design Bounds",
    subtitle = paste0(
      "Final Information = ", round(final_info, 2)
    )
  )
```

| Updated Group Sequential Design Bounds |                     |          |
|----------------------------------------|---------------------|----------|
| Final Information = 71.51              |                     |          |
| Analysis                               | Value               | Efficacy |
| IA 1: 17%                              | Z                   | 7.1773   |
| Information: 12.06                     | p (1-sided)         | 0.0000   |
|                                        | ~RR at bound        | 1.8918   |
|                                        | P(Cross) if RR=1    | 0.0000   |
|                                        | P(Cross) if RR=2.72 | 1.0000   |
| Final                                  | Z                   | 1.9600   |
| Information: 71.51                     | p (1-sided)         | 0.0250   |
|                                        | ~RR at bound        | 1.0741   |
|                                        | P(Cross) if RR=1    | 0.0250   |
|                                        | P(Cross) if RR=2.72 | 1.0000   |

``` r

# Test statistic
z_score <- coef(fit_final)["treatmentExperimental"] / sqrt(var_beta)
# Note: gsDesign assumes Z > 0 for efficacy. 
# Our model gives log(RR) < 0 for efficacy.
# So Z_stat = - (log(RR)) / SE
z_stat <- -coef(fit_final)["treatmentExperimental"] / sqrt(var_beta)

cat("Z-statistic:", z_stat, "\n")
#> Z-statistic: 1.24208
cat("Final Bound:", gs_update$upper$bound[2], "\n")
#> Final Bound: 1.959964

if (z_stat > gs_update$upper$bound[2]) {
  cat("Reject Null Hypothesis\n")
} else {
  cat("Fail to Reject Null Hypothesis\n")
}
#> Fail to Reject Null Hypothesis
```

## References
