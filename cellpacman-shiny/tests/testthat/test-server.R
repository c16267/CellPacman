test_that("app server loads feature and metadata uploads through package API", {
  old_plan <- future::plan(future::multisession, workers = I(1L))
  on.exit(future::plan(old_plan), add = TRUE)
  feature_path <- tempfile(fileext = ".csv")
  metadata_path <- tempfile(fileext = ".csv")
  readr::write_csv(feature_fixture(), feature_path)
  readr::write_csv(metadata_fixture(), metadata_path)

  shiny::testServer(app_server, {
    session$setInputs(
      feature_file = data.frame(
        name = "features.csv",
        datapath = feature_path,
        stringsAsFactors = FALSE
      ),
      metadata_file = data.frame(
        name = "metadata.csv",
        datapath = metadata_path,
        stringsAsFactors = FALSE
      )
    )

    wait_for_app(data_table_ready, session)
    expect_s3_class(data(), "cellpainting_data")
    expect_equal(nrow(display_data()), nrow(feature_fixture()))
  })
})

test_that("the app object is built with validated settings", {
  expect_s3_class(cellpacman_app(use_example_data = TRUE), "shiny.appobj")
  expect_error(cellpacman_app(use_example_data = "yes"), "TRUE or FALSE")
  expect_error(cellpacman_app(workers = 0), "positive integer")
})

test_that("runtime configuration is read from the environment", {
  env <- function(x, unset) switch(x, CELLPACMAN_EXAMPLE_DATA = "TRUE", CELLPACMAN_WORKERS = "3", unset)
  expect_identical(app_config(env), list(use_example_data = TRUE, workers = 3L))
  env <- function(x, unset) switch(x, CELLPACMAN_WORKERS = "not-a-number", unset)
  expect_identical(app_config(env), list(use_example_data = FALSE, workers = 1L))
})

test_that("the UI attaches the app's static assets", {
  dep <- app_html_dependency()
  expect_s3_class(dep, "html_dependency")
  expect_true(file.exists(file.path(dep$src$file, dep$stylesheet)))
  expect_true(file.exists(file.path(dep$src$file, dep$script)))
  names <- vapply(htmltools::findDependencies(app_ui()), function(d) d$name, character(1))
  expect_true("cellpacman-app" %in% names)
})

test_that("example tables are read through the package API", {
  features <- app_read_table(NULL, "features")
  metadata <- app_read_table(NULL, "metadata")
  expect_identical(names(features)[1], "WellId")
  expect_true(all(c("WellId", "Plate", "Compound", "Concentration") %in% names(metadata)))
  expect_equal(nrow(features), nrow(metadata))
})
