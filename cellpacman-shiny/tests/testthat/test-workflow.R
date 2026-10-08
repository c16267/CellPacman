test_that("background workflow preserves analysis results and resets after uploads", {
  old_plan <- future::plan(future::multisession, workers = I(1L))
  on.exit(future::plan(old_plan), add = TRUE)
  features <- tibble::tibble(
    WellId = paste0("W", seq_len(36)),
    FeatureA = sin(seq_len(36)), FeatureB = cos(seq_len(36)),
    FeatureC = seq_len(36)
  )
  metadata <- tibble::tibble(
    WellId = features$WellId, Plate = "P1",
    Compound = rep(c("DrugA", "DrugB"), each = 18),
    Concentration = rep(c(0.1, 1, 10), 12)
  )
  feature_path <- tempfile(fileext = ".csv")
  metadata_path <- tempfile(fileext = ".csv")
  on.exit(unlink(c(feature_path, metadata_path)), add = TRUE)
  readr::write_csv(features, feature_path)
  readr::write_csv(metadata, metadata_path)
  expected <- app_dimred_analysis(loadCellPainting(feature_path, metadata_path))
  expect_named(dplyr::select(expected$dimred$reduced_data, x, y), c("x", "y"))
  expect_identical(expected$dimred$reduction$n_dim, 2L)

  shiny::testServer(function(input, output, session) {
    workflow <- app_workflow_server(input, session, FALSE)
  }, {
    session$setInputs(
      feature_file = data.frame(name = "features.csv", datapath = feature_path),
      metadata_file = data.frame(name = "metadata.csv", datapath = metadata_path),
      cluster_method = "dbscan", cluster_eps = 100, cluster_min_pts = 2,
      curve_method = "princurve", curve_smoother = "lowess"
    )
    wait_for_app(workflow$data_ready, session)
    session$setInputs(run_dimred_button = 1)
    wait_for_app(workflow$dimred_ready, session)
    expect_equal(workflow$dimred_results()$dimred$reduced_data, expected$dimred$reduced_data)
    session$setInputs(run_clustering_button = 1)
    wait_for_app(workflow$clustering_ready, session)
    expect_s3_class(workflow$clustering_results()$clusters, "cellpainting_clusters")
    session$setInputs(run_trajectory_button = 1)
    wait_for_app(workflow$trajectory_ready, session)
    expect_null(workflow$trajectory_results()$error)

    # The second request has identical inputs and is ready without another job.
    session$setInputs(run_dimred_button = 2)
    expect_true(workflow$dimred_ready())
    expect_false(workflow$clustering_ready())
    expect_false(workflow$trajectory_ready())

    session$setInputs(metadata_file = data.frame(name = "replacement.csv", datapath = metadata_path))
    expect_false(workflow$data_ready())
    expect_false(workflow$dimred_ready())
    expect_error(workflow$dimred_results(), class = "shiny.silent.error")
    wait_for_app(workflow$data_ready, session)
  })
})

test_that("upload validation failures are displayed and recover on replacement", {
  feature_path <- tempfile(fileext = ".csv")
  metadata_path <- tempfile(fileext = ".csv")
  on.exit(unlink(c(feature_path, metadata_path)), add = TRUE)
  readr::write_csv(feature_fixture(), feature_path)
  readr::write_csv(dplyr::select(metadata_fixture(), -Plate), metadata_path)
  shiny::testServer(function(input, output, session) {
    workflow <- app_workflow_server(input, session, FALSE)
  }, {
    session$setInputs(
      feature_file = data.frame(name = "features.csv", datapath = feature_path),
      metadata_file = data.frame(name = "metadata.csv", datapath = metadata_path)
    )
    wait_for_app(function() workflow$load_status() == "error", session)
    expect_false(workflow$data_ready())
    expect_error(workflow$load_message(), "Plate")
    readr::write_csv(metadata_fixture(), metadata_path)
    session$setInputs(metadata_file = data.frame(name = "fixed.csv", datapath = metadata_path))
    wait_for_app(workflow$data_ready, session)
    expect_identical(workflow$load_message(), "Data ready.")
  })
})
