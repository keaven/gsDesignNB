# Preview built pkgdown site in the browser

[`pkgdown::preview_site()`](https://pkgdown.r-lib.org/reference/preview_site.html)
opens `docs/index.html` as a `file://` URL. Many browsers restrict
loading stylesheets, scripts, and fonts from local files, so the site
can appear almost unstyled. This function serves `docs/` over HTTP on
the loopback interface so the site matches the published GitHub Pages
appearance.

## Usage

``` r
preview_pkgdown_site(port = 8787L)
```

## Arguments

- port:

  Integer; TCP port for the static server (default `8787`).

## Value

Invisibly, `NULL` (called for its side effect of running the server).

## Details

Call from the **package root** after
[`pkgdown::build_site()`](https://pkgdown.r-lib.org/reference/build_site.html)
(or any build that populates `docs/`).

The server runs until you interrupt the R session (e.g. **Esc** in
RStudio or Ctrl+C in the terminal).

## Examples

``` r
if (FALSE) { # \dontrun{
pkgdown::build_site()
preview_pkgdown_site()
} # }
```
