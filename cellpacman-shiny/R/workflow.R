# The workflow handles task execution and invalidation; renderers only consume
# completed results. Scientific operations remain in the package API.
app_workflow_server <- function(input, session, use_example_data) {
  feature_task <- app_task(app_read_table, session)
  metadata_task <- app_task(app_read_table, session)
  data_task <- app_task(cellpacman::loadCellPainting, session)
  dimred_task <- app_task(app_dimred_analysis, session)
  clustering_task <- app_task(app_cluster_analysis, session)
  trajectory_task <- app_task(app_trajectory_analysis, session)
  dataset_id <- 0L
  dimred_id <- 0L
  
  invalidate_analysis <- function() {
    dimred_task$invalidate()
    clustering_task$invalidate()
    trajectory_task$invalidate()
  }
  
  shiny::observeEvent(list(input$feature_file, input$metadata_file), {
    dataset_id <<- dataset_id + 1L
    data_task$invalidate()
    invalidate_analysis()
  }, ignoreNULL = FALSE, priority = 100)
  
  shiny::observeEvent(input$feature_file, {
    if (is.null(input$feature_file) && !use_example_data) {
      feature_task$invalidate()
    } else {
      feature_task$invoke(
        list(upload = input$feature_file, example = "features"),
        key = input$feature_file
      )
    }
  }, ignoreNULL = FALSE)
  
  shiny::observeEvent(input$metadata_file, {
    if (is.null(input$metadata_file) && !use_example_data) {
      metadata_task$invalidate()
    } else {
      metadata_task$invoke(
        list(upload = input$metadata_file, example = "metadata"),
        key = input$metadata_file
      )
    }
  }, ignoreNULL = FALSE)
  
  shiny::observeEvent(list(feature_task$result(), metadata_task$result()), {
    data_task$invoke(
      list(file.feature = feature_task$result(), file.metadata = metadata_task$result()),
      key = dataset_id
    )
  })
  
  # A task button acknowledges the click in the browser before any R work.
  # Explicit reset also handles validation failures and cached results.
  run_button <- function(id, task, ready, run) {
    shiny::observe({
      bslib::update_task_button(
        id, state = if (task$status() == "running") "busy" else "ready", session = session
      )
    })
    shiny::observeEvent(input[[id]], {
      on.exit({
        bslib::update_task_button(
          id, state = if (task$status() == "running") "busy" else "ready", session = session
        )
      })
      shiny::req(ready(), task$status() != "running")
      run()
    })
  }
  
  data_ready <- shiny::reactive(data_task$status() == "ready")
  dimred_ready <- shiny::reactive(dimred_task$status() == "ready")
  clustering_ready <- shiny::reactive(clustering_task$status() == "ready")
  trajectory_ready <- shiny::reactive(trajectory_task$status() == "ready")
  
  run_button("run_dimred_button", dimred_task, data_ready, function() {
    dimred_id <<- dimred_id + 1L
    clustering_task$invalidate()
    trajectory_task$invalidate()
    dimred_task$invoke(list(data = data_task$result()), key = dataset_id)
  })
  
  run_button("run_clustering_button", clustering_task, dimred_ready, function() {
    trajectory_task$invalidate()
    args <- list(
      dimred = dimred_task$result()$dimred,
      method = input$cluster_method, eps = input$cluster_eps,
      min_pts = as.integer(input$cluster_min_pts)
    )
    clustering_task$invoke(args, key = list(dimred_id, args$method, args$eps, args$min_pts))
  })
  
  run_button("run_trajectory_button", trajectory_task, clustering_ready, function() {
    args <- list(
      dimred = dimred_task$result()$dimred,
      data_norm = dimred_task$result()$data_norm,
      method = input$curve_method, eps = input$cluster_eps,
      min_pts = as.integer(input$cluster_min_pts), smoother = input$curve_smoother
    )
    trajectory_task$invoke(args, key = list(dimred_id, args$method, args$eps, args$min_pts, args$smoother))
  })
  
  # Loading/validation has no analysis button, so keep its status visible too.
  load_status <- shiny::reactive({
    statuses <- c(feature_task$status(), metadata_task$status(), data_task$status())
    if ("error" %in% statuses) return("error")
    if ("running" %in% statuses) return("running")
    if (data_ready()) return("ready")
    "empty"
  })
  
  list(
    feature_data = feature_task$result, 
    metadata_data = metadata_task$result,
    data = data_task$result, 
    dimred_results = dimred_task$result,
    clustering_results = clustering_task$result, 
    trajectory_results = trajectory_task$result,
    data_ready = data_ready, 
    dimred_ready = dimred_ready,
    clustering_ready = clustering_ready, 
    trajectory_ready = trajectory_ready,
    load_status = load_status,
    load_message = shiny::reactive({
      # Surface both parser and joint validation errors next to the uploads.
      if (load_status() == "error") {
        for (task in list(feature_task, metadata_task, data_task)) {
          if (task$status() == "error") task$result()
        }
      }
      switch(load_status(), running = "Loading and validating data...",
             ready = "Data ready.", empty = "", error = "Unable to load data.")
    })
  )
}
