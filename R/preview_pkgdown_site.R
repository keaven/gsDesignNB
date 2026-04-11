#' Preview built pkgdown site in the browser
#'
#' `pkgdown::preview_site()` opens `docs/index.html` as a `file://` URL. Many
#' browsers restrict loading stylesheets, scripts, and fonts from local files, so
#' the site can appear almost unstyled. This function serves `docs/` over HTTP
#' on the loopback interface so the site matches the published GitHub Pages
#' appearance.
#'
#' Call from the **package root** after [pkgdown::build_site()] (or any build
#' that populates `docs/`).
#'
#' The server runs until you interrupt the R session (e.g. **Esc** in RStudio or
#' Ctrl+C in the terminal).
#'
#' @param port Integer; TCP port for the static server (default `8787`).
#'
#' @return Invisibly, `NULL` (called for its side effect of running the server).
#'
#' @examples
#' \dontrun{
#' pkgdown::build_site()
#' preview_pkgdown_site()
#' }
#'
#' @export
preview_pkgdown_site <- function(port = 8787L) { # nocov start
  docs <- file.path(getwd(), "docs")
  if (!dir.exists(docs)) {
    stop(
      "No docs/ directory in the working directory. Run pkgdown::build_site() from the package root first.",
      call. = FALSE
    )
  }
  if (!requireNamespace("httpuv", quietly = TRUE)) {
    stop("Install the httpuv package: install.packages(\"httpuv\")", call. = FALSE)
  }
  port <- as.integer(port)[1L]
  httpuv::runStaticServer(
    normalizePath(docs, winslash = "/", mustWork = TRUE),
    host = "127.0.0.1",
    port = port,
    browse = interactive(),
    background = FALSE
  )
  invisible(NULL)
} # nocov end
