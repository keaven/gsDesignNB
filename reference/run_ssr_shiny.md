# Launch the SSR Shiny prototype

Opens a lightweight Shiny interface for interactive exploration of
adaptive sample size re-estimation (SSR) scenarios. The app is a thin
wrapper over
[`sample_size_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sample_size_nbinom.md),
[`gsNBCalendar()`](https://keaven.github.io/gsDesignNB/reference/gsNBCalendar.md),
[`sim_ssr_nbinom()`](https://keaven.github.io/gsDesignNB/reference/sim_ssr_nbinom.md),
and
[`summarize_ssr_sim()`](https://keaven.github.io/gsDesignNB/reference/summarize_ssr_sim.md),
so the statistical computations remain in the package.

## Usage

``` r
run_ssr_shiny(
  display.mode = c("normal", "showcase"),
  launch.browser = interactive()
)
```

## Arguments

- display.mode:

  Character; passed to
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html).
  Defaults to `"normal"`.

- launch.browser:

  Logical; passed to
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html).
  Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html).

## Value

Invisibly returns the Shiny app object.

## Examples

``` r
if (FALSE) { # \dontrun{
run_ssr_shiny()
} # }
```
