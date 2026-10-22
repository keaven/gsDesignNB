if (requireNamespace("shiny", quietly = TRUE)) assert("SSR Shiny prototype sources successfully", {
  app_path <- system.file(
    "shiny-examples",
    "ssr-prototype",
    "app.R",
    package = "gsDesignNB"
  )
  if (!nzchar(app_path)) {
    app_path <- file.path("inst", "shiny-examples", "ssr-prototype", "app.R")
  }

  (file.exists(app_path))

  app_obj <- source(app_path, local = new.env(parent = baseenv()))$value
  (inherits(app_obj, "shiny.appobj"))
})
