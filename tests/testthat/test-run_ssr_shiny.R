test_that("SSR Shiny prototype sources successfully", {
  skip_if_not_installed("shiny")

  app_path <- system.file(
    "shiny-examples",
    "ssr-prototype",
    "app.R",
    package = "gsDesignNB"
  )
  if (!nzchar(app_path)) {
    app_path <- file.path("inst", "shiny-examples", "ssr-prototype", "app.R")
  }

  expect_true(file.exists(app_path))

  app_obj <- source(app_path, local = new.env(parent = baseenv()))$value
  expect_s3_class(app_obj, "shiny.appobj")
})
