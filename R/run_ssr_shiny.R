#' Launch the SSR Shiny prototype
#'
#' Opens a lightweight Shiny interface for interactive exploration of adaptive
#' sample size re-estimation (SSR) scenarios. The app is a thin wrapper over
#' [sample_size_nbinom()], [gsNBCalendar()], [sim_ssr_nbinom()], and
#' [summarize_ssr_sim()], so the statistical computations remain in the package.
#'
#' @param display.mode Character; passed to [shiny::runApp()]. Defaults to
#'   `"normal"`.
#' @param launch.browser Logical; passed to [shiny::runApp()]. Defaults to
#'   `interactive()`.
#'
#' @return Invisibly returns the Shiny app object.
#'
#' @examples
#' \dontrun{
#' run_ssr_shiny()
#' }
#'
#' @export
run_ssr_shiny <- function( # nocov start
  display.mode = c("normal", "showcase"),
  launch.browser = interactive()
) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Install the shiny package: install.packages(\"shiny\")", call. = FALSE)
  }

  display.mode <- match.arg(display.mode)

  app_dir <- system.file("shiny-examples", "ssr-prototype", package = "gsDesignNB")
  if (!nzchar(app_dir)) {
    candidates <- c(
      file.path(getwd(), "inst", "shiny-examples", "ssr-prototype"),
      file.path(getwd(), "..", "inst", "shiny-examples", "ssr-prototype")
    )
    existing <- candidates[dir.exists(candidates)]
    if (length(existing) > 0) {
      app_dir <- normalizePath(existing[1], winslash = "/", mustWork = TRUE)
    }
  }

  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    stop(
      "Could not locate the SSR Shiny prototype in inst/shiny-examples/ssr-prototype.",
      call. = FALSE
    )
  }

  shiny::runApp(
    appDir = app_dir,
    display.mode = display.mode,
    launch.browser = launch.browser
  )
} # nocov end
