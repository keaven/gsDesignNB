# gsDesignNB

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/gsDesignNB)](https://CRAN.R-project.org/package=gsDesignNB)
[![Codecov test coverage](https://codecov.io/gh/keaven/gsDesignNB/graph/badge.svg)](https://app.codecov.io/gh/keaven/gsDesignNB)
[![pkgdown](https://github.com/keaven/gsDesignNB/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/keaven/gsDesignNB/actions/workflows/pkgdown.yaml)
[![R-CMD-check](https://github.com/keaven/gsDesignNB/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/keaven/gsDesignNB/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

gsDesignNB provides design, simulation, and interim monitoring tools for
recurrent-event trials analyzed with negative binomial rate models, with
Poisson methods available as the limiting special case when dispersion is
negligible.

The package is NB-first: plan designs with `sample_size_nbinom()`, simulate
recurrent-event data with `nb_sim()` or `nb_sim_seasonal()`, and evaluate
group sequential monitoring or sample size re-estimation with
`sim_gs_nbinom()` and `sim_ssr_nbinom()`. Planning and simulation can use either
Wald or score-test inference for rate ratios. The `gsDesign` package supplies
the underlying spending-function and boundary calculations used by those
workflows.

## Start Here

- `sample-size-nbinom` for fixed-design planning
- `score-vs-wald-simulation` for Wald/score sizing and Type I error guidance
- `simulation-example` and `seasonal-simulation` for recurrent-event data generation
- `ssr-example` and `ssr-simulation-study` for negative binomial SSR workflows

## Installation

You can install gsDesignNB from CRAN with:

```r
install.packages("gsDesignNB")
```

Or install the development version from GitHub with:

```r
remotes::install_github("keaven/gsDesignNB")
```

## Local pkgdown site

After `pkgdown::build_site()`, prefer **not** using `pkgdown::preview_site()` alone:
it opens a `file://` URL, and browsers often load little or no CSS/JS for local
files. From the package root, run:

```r
gsDesignNB::preview_pkgdown_site()
```

That serves `docs/` at `http://127.0.0.1:8787` so the site matches GitHub Pages.

## Code style

This package follows the [tidyverse style guide](https://style.tidyverse.org/).
