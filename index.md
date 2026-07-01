# gsDesignNB

gsDesignNB provides design, simulation, and interim monitoring tools for
recurrent-event trials analyzed with negative binomial rate models, with
Poisson methods available as the limiting special case when dispersion
is negligible.

The package is NB-first: plan designs with
[`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md),
simulate recurrent-event data with
[`nb_sim()`](https://keaven.github.io/gsDesignNB/reference/nb_sim.md) or
[`nb_sim_seasonal()`](https://keaven.github.io/gsDesignNB/reference/nb_sim_seasonal.md),
and evaluate group sequential monitoring or sample size re-estimation
with
[`sim_gs_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_gs_nbinom.md)
and
[`sim_ssr_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_ssr_nbinom.md).
Planning and simulation can use either Wald or score-test inference for
rate ratios. The gsDesign package supplies the underlying
spending-function and boundary calculations used by those workflows.

## Start Here

- `sample-size-nbinom` for fixed-design planning
- `score-vs-wald-simulation` for Wald/score sizing and Type I error
  guidance
- `simulation-example` and `seasonal-simulation` for recurrent-event
  data generation
- `ssr-example` and `ssr-simulation-study` for negative binomial SSR
  workflows

## Installation

You can install gsDesignNB from CRAN with:

``` r

install.packages("gsDesignNB")
```

Or install the development version from GitHub with:

``` r

remotes::install_github("keaven/gsDesignNB")
```

## Code style

This package follows the [tidyverse style
guide](https://style.tidyverse.org/).
