test_that("uploads and analysis acknowledge work before it finishes in the browser", {
  skip_on_cran()
  skip_if_not_installed("shinytest2")
  app <- shinytest2::AppDriver$new(
    test_path("apps", "responsive"), name = "responsiveness",
    load_timeout = 30000, timeout = 30000
  )
  on.exit(app$stop(), add = TRUE)
  feature_path <- tempfile(fileext = ".csv")
  metadata_path <- tempfile(fileext = ".csv")
  on.exit(unlink(c(feature_path, metadata_path)), add = TRUE)
  readr::write_csv(tibble::tibble(
    WellId = paste0("W", seq_len(36)),
    FeatureA = sin(seq_len(36)), FeatureB = cos(seq_len(36)), FeatureC = seq_len(36)
  ), feature_path)
  readr::write_csv(tibble::tibble(
    WellId = paste0("W", seq_len(36)), Plate = "P1",
    Compound = rep(c("DrugA", "DrugB"), each = 18),
    Concentration = rep(c(0.1, 1, 10), 12)
  ), metadata_path)
  app$upload_file(feature_file = feature_path, wait_ = FALSE)
  app$upload_file(metadata_file = metadata_path, wait_ = FALSE)
  app$wait_for_js("document.querySelector('#load_status').textContent.includes('Loading')")
  app$wait_for_js("!document.querySelector('#feature_table').closest('.shiny-spinner-output-container').querySelector('.load-container').classList.contains('shiny-spinner-hidden')")
  app$run_js("Shiny.setInputValue('ping_request', 'during-upload')")
  app$wait_for_js("document.querySelector('#ping_response').textContent === 'during-upload'")
  expect_true(app$get_js("document.querySelector('#load_status').textContent.includes('Loading')"))
  app$wait_for_js("document.querySelector('#load_status').textContent === 'Data ready.'")

  for (stage in c("dimred", "clustering", "trajectory")) {
    tab <- switch(stage, dimred = "dimension_reduction", clustering = "clustering", trajectory = "trajectory")
    plot <- switch(stage, dimred = "umap_plot", clustering = "cluster_plot", trajectory = "trajectory_curve_plot")
    button <- paste0("run_", stage, "_button")
    app$run_js(sprintf("document.querySelector('#main_nav a[data-value=\"%s\"]').click()", tab))
    # Read the button in the same browser turn as the click: the busy state
    # must not require a server round-trip.
    expect_true(app$get_js(sprintf(
      "(() => {const b = document.getElementById('%s'); b.click(); return b.disabled;})()", button
    )))
    app$wait_for_js(sprintf(
      "!document.getElementById('%s').closest('.shiny-spinner-output-container').querySelector('.load-container').classList.contains('shiny-spinner-hidden')", plot
    ))
    app$run_js(sprintf("Shiny.setInputValue('ping_request', '%s')", stage))
    app$wait_for_js(sprintf("document.querySelector('#ping_response').textContent === '%s'", stage))
    expect_true(app$get_js(sprintf("document.getElementById('%s').disabled", button)))
    app$wait_for_js(sprintf("!document.getElementById('%s').disabled", button))
  }
})
