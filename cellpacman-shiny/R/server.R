# Server: consumes the completed stage results exposed by app_workflow_server()
# and renders tables and plots. No analysis is performed here; the
# comparison panels call the cellpacman API directly on the normalized data.
app_server <- function(input, output, session, use_example_data = FALSE) {
  workflow <- app_workflow_server(input, session, use_example_data)
  data_table_ready <- workflow$data_ready
  dimred_ready <- workflow$dimred_ready
  clustering_ready <- workflow$clustering_ready
  trajectory_ready <- workflow$trajectory_ready
  feature_data <- workflow$feature_data
  metadata_data <- workflow$metadata_data
  data <- workflow$data
  dimred_results <- workflow$dimred_results
  clustering_results <- workflow$clustering_results
  trajectory_results <- workflow$trajectory_results
  
  output$load_status <- shiny::renderText(workflow$load_message())
  
  set_nav_disabled <- function(value, disabled) {
    session$sendCustomMessage(
      "setNavDisabled",
      list(
        navId = "main_nav",
        value = value,
        disabled = isTRUE(disabled)
      )
    )
  }
  
  shiny::observe({
    set_nav_disabled("dimension_reduction", !isTRUE(data_table_ready()))
    set_nav_disabled("clustering", !isTRUE(dimred_ready()))
    set_nav_disabled("trajectory", !isTRUE(clustering_ready()))
  })
  
  app_pubr_minimal_theme <- function() {
    if (requireNamespace("ggpubr", quietly = TRUE)) {
      exports <- getNamespaceExports("ggpubr")
      if ("pubr_minimal" %in% exports) {
        theme_fun <- utils::getFromNamespace("pubr_minimal", "ggpubr")
        return(tryCatch(theme_fun(base_size = 12), error = function(e) theme_fun()))
      }
      if ("theme_pubr" %in% exports) {
        return(ggpubr::theme_pubr(base_size = 12))
      }
    }
    
    ggplot2::theme_minimal(base_size = 12)
  }
  
  style_umap_plot <- function(plot, plot_limits, legend_position = "right") {
    plot$coordinates <- ggplot2::coord_cartesian(
      xlim = plot_limits$x,
      ylim = plot_limits$y,
      expand = FALSE
    )
    
    plot +
      app_pubr_minimal_theme() +
      ggplot2::theme(legend.position = legend_position)
  }
  
  display_data <- shiny::reactive({
    shiny::req(data())
    dplyr::left_join(data()$metadata, data()$feature_data, by = "WellId")
  })
  
  table_options <- list(
    wrap = FALSE,
    searchable = TRUE,
    pagination = TRUE,
    highlight = TRUE,
    striped = TRUE,
    bordered = TRUE
  )
  
  feature_comparison_table <- function(results) {
    numeric_col <- function(digits = 3, scientific = FALSE) {
      reactable::colDef(
        cell = function(value) {
          if (is.na(value)) {
            return(shiny::div(""))
          }
          
          full_value <- format(value, digits = 16, scientific = TRUE, trim = TRUE)
          displayed_value <- format(value, digits = digits, scientific = scientific, trim = TRUE)
          shiny::div(title = full_value, displayed_value)
        }
      )
    }
    
    columns <- list(
      feature = reactable::colDef(
        minWidth = 220,
        cell = function(value) {
          shiny::div(
            title = as.character(value),
            style = "max-width: 20rem; overflow: hidden; text-overflow: ellipsis; white-space: nowrap;",
            value
          )
        }
      )
    )
    
    for (column in intersect(c("p_value", "adjusted_p_value"), names(results))) {
      columns[[column]] <- numeric_col(digits = 3, scientific = TRUE)
    }
    for (column in intersect("median_difference", names(results))) {
      columns[[column]] <- numeric_col(digits = 3)
    }
    
    reactable::reactable(
      results,
      columns = columns,
      defaultColDef = reactable::colDef(minWidth = 90, resizable = TRUE),
      compact = TRUE,
      style = list(overflowX = "auto", fontSize = "0.8rem"),
      wrap = FALSE,
      searchable = TRUE,
      pagination = TRUE,
      highlight = TRUE,
      striped = TRUE,
      bordered = TRUE
    )
  }
  
  output$feature_table <- reactable::renderReactable({
    do.call(reactable::reactable, c(list(feature_data()), table_options))
  })
  
  output$metadata_table <- reactable::renderReactable({
    do.call(reactable::reactable, c(list(metadata_data()), table_options))
  })
  
  cluster_display_choices <- shiny::reactive({
    shiny::req(clustering_results())
    
    clusters <- clustering_results()$clusters
    choices <- character(0)
    
    if (!is.null(clusters$unknown_assignment) && nrow(clusters$unknown_assignment) > 0) {
      choices <- c("Unknowns" = "__unknown__")
    }
    
    compound_choices <- names(clusters$known_assignments)
    if (length(compound_choices) > 0) {
      choices <- c(choices, stats::setNames(compound_choices, compound_choices))
    }
    
    choices
  })
  
  output$umap_plot <- plotly::renderPlotly({
    results <- dimred_results()
    plot <- cellpacman::plotDimred(results$dimred) +
      ggplot2::labs(title = NULL, subtitle = NULL)
    plot |>
      style_umap_plot(results$plot_limits, legend_position = "top") |>
      plotly::ggplotly() |>
      plotly::config(responsive = TRUE)
  })
  
  compound_comparison_choices <- shiny::reactive({
    shiny::req(dimred_results())
    compound_values <- data()$metadata$Compound
    compounds <- sort(unique(as.character(compound_values[!is.na(compound_values)])))
    choices <- stats::setNames(compounds, compounds)
    if (any(is.na(compound_values))) {
      choices <- c(choices, "Unknown" = "__unknown__")
    }
    choices
  })
  
  output$compound_comparison_controls <- shiny::renderUI({
    if (!isTRUE(dimred_ready())) {
      return(NULL)
    }
    
    choices <- compound_comparison_choices()
    if (length(choices) < 2) {
      return(shiny::p("At least two compound groups are required for comparison."))
    }
    
    shiny::tagList(
      shiny::selectInput("first_compound", "First compound", choices = choices, selected = choices[[1]]),
      shiny::selectInput("second_compound", "Second compound", choices = choices, selected = choices[[2]])
    )
  })
  
  shiny::observeEvent(compound_comparison_choices(), {
    choices <- compound_comparison_choices()
    if (length(choices) < 2 || !is.null(input$first_compound)) {
      return()
    }
    
    session$onFlushed(function() {
      shiny::updateSelectInput(
        session,
        "first_compound",
        selected = unname(choices[[1]])
      )
      shiny::updateSelectInput(
        session,
        "second_compound",
        selected = unname(choices[[2]])
      )
    }, once = TRUE)
  }, ignoreInit = FALSE)
  
  compound_comparison <- shiny::reactive({
    shiny::req(dimred_results(), input$first_compound, input$second_compound)
    shiny::validate(shiny::need(
      !identical(input$first_compound, input$second_compound),
      "Select two different compounds to compare."
    ))
    cellpacman::compareCompoundFeatures(
      dimred_results()$data_norm,
      input$first_compound,
      input$second_compound
    )
  })
  
  shiny::observeEvent(input$first_compound, {
    choices <- compound_comparison_choices()
    available_second_choices <- choices[unname(choices) != input$first_compound]
    if (length(available_second_choices) == 0) {
      return()
    }
    
    selected <- input$second_compound
    if (is.null(selected) || selected == input$first_compound || !selected %in% unname(available_second_choices)) {
      selected <- unname(available_second_choices[[1]])
    }
    shiny::updateSelectInput(session, "second_compound", choices = available_second_choices, selected = selected)
  }, ignoreInit = TRUE)
  
  output$compound_comparison_table_panel <- shiny::renderUI({
    if (!isTRUE(dimred_ready())) {
      return(NULL)
    }
    
    shinycssloaders::withSpinner(
      reactable::reactableOutput("compound_comparison_table", height = "100%"),
      type = 1,
      color = "#18344d",
      proxy.height = "24rem"
    )
  })
  
  output$compound_comparison_table <- reactable::renderReactable({
    results <- compound_comparison()
    feature_comparison_table(dplyr::select(results, -significance))
  })
  
  output$compound_comparison_volcano_panel <- shiny::renderUI({
    if (!isTRUE(dimred_ready())) {
      return(NULL)
    }
    
    trajectory_plot_output("compound_comparison_volcano")
  })
  
  output$compound_comparison_volcano <- plotly::renderPlotly({
    results <- compound_comparison()
    shiny::validate(shiny::need(nrow(results) > 0, "No comparable feature values are available."))
    title <- paste(
      if (identical(input$second_compound, "__unknown__")) "Unknown" else input$second_compound,
      "versus",
      if (identical(input$first_compound, "__unknown__")) "Unknown" else input$first_compound
    )
    app_volcano_plot(results, "median_difference", "adjusted_p_value", title)
  })
  
  output$compound_comparison_heatmap_panel <- shiny::renderUI({
    if (!isTRUE(dimred_ready())) {
      return(NULL)
    }
    
    trajectory_plot_output("compound_comparison_heatmap")
  })
  
  output$compound_comparison_heatmap <- plotly::renderPlotly({
    results <- compound_comparison()
    significant_features <- results$feature[results$adjusted_p_value < 0.05]
    shiny::validate(shiny::need(
      length(significant_features) > 0,
      "No features meet the FDR < 0.05 threshold for this comparison."
    ))
    app_heatmap_plot(
      cellpacman::summarizeCompoundFeatures(dimred_results()$data_norm, significant_features),
      "compound", "Compound", "Median\nnormalized value", rotate_x = TRUE
    )
  })
  
  output$cluster_plot <- plotly::renderPlotly({
    clusters <- clustering_results()$clusters
    shiny::validate(shiny::need(length(cluster_display_choices()) > 0,
                                "No clusterable compounds detected. Try adjusting DBSCAN settings."))
    shiny::req(input$comparison_compound)
    cluster_plot <- cellpacman::plotKnownClusters(clusters, input$comparison_compound)
    
    shiny::validate(
      shiny::need(
        !is.null(cluster_plot),
        "No clusters were detected for the selected display."
      )
    )
    
    cluster_plot |>
      style_umap_plot(dimred_results()$plot_limits) |>
      plotly::ggplotly() |>
      plotly::config(responsive = TRUE)
  })
  
  cluster_comparison_compounds <- shiny::reactive({
    shiny::req(clustering_results())
    compounds <- names(clustering_results()$clusters$known_assignments)
    stats::setNames(compounds, compounds)
  })
  
  selected_cluster_assignments <- shiny::reactive({
    shiny::req(clustering_results(), input$comparison_compound)
    clustering_results()$clusters$known_assignments[[input$comparison_compound]]
  })
  
  cluster_comparison_choices <- shiny::reactive({
    assignments <- selected_cluster_assignments()
    clusters <- sort(unique(as.character(assignments$cluster[assignments$cluster != "0"])))
    stats::setNames(clusters, paste("Cluster", clusters))
  })
  
  output$cluster_comparison_controls <- shiny::renderUI({
    if (!isTRUE(clustering_ready())) {
      return(NULL)
    }
    
    compounds <- cluster_comparison_compounds()
    if (length(compounds) == 0) {
      return(shiny::p("No known compound clusters are available for comparison."))
    }
    
    selected_compound <- input$comparison_compound
    if (is.null(selected_compound) || !selected_compound %in% unname(compounds)) {
      selected_compound <- unname(compounds[[1]])
    }
    assignments <- clustering_results()$clusters$known_assignments[[selected_compound]]
    clusters <- sort(unique(as.character(assignments$cluster[assignments$cluster != "0"])))
    if (length(clusters) < 2) {
      return(shiny::tagList(
        shiny::selectInput("comparison_compound", "Compound", choices = compounds, selected = selected_compound),
        shiny::p("The selected compound needs at least two non-noise clusters.")
      ))
    }
    
    cluster_choices <- stats::setNames(clusters, paste("Cluster", clusters))
    shiny::tagList(
      shiny::selectInput("comparison_compound", "Compound", choices = compounds, selected = selected_compound),
      shiny::selectInput("first_cluster", "First cluster", choices = cluster_choices, selected = clusters[[1]]),
      shiny::selectInput("second_cluster", "Second cluster", choices = cluster_choices, selected = clusters[[2]])
    )
  })
  
  shiny::observeEvent(input$comparison_compound, {
    choices <- cluster_comparison_choices()
    if (length(choices) < 2) {
      return()
    }
    session$onFlushed(function() {
      shiny::updateSelectInput(session, "first_cluster", choices = choices, selected = unname(choices[[1]]))
      shiny::updateSelectInput(session, "second_cluster", choices = choices, selected = unname(choices[[2]]))
    }, once = TRUE)
  }, ignoreInit = TRUE)
  
  cluster_comparison <- shiny::reactive({
    shiny::req(dimred_results(), input$first_cluster, input$second_cluster)
    shiny::validate(shiny::need(!identical(input$first_cluster, input$second_cluster), "Select two different clusters to compare."))
    cellpacman::compareClusterFeatures(
      dimred_results()$data_norm,
      selected_cluster_assignments(),
      input$first_cluster,
      input$second_cluster
    )
  })
  
  output$cluster_comparison_table_panel <- shiny::renderUI({
    if (!isTRUE(clustering_ready())) {
      return(NULL)
    }
    shinycssloaders::withSpinner(
      reactable::reactableOutput("cluster_comparison_table", height = "100%"),
      type = 1,
      color = "#18344d",
      proxy.height = "24rem"
    )
  })
  
  output$cluster_comparison_table <- reactable::renderReactable({
    feature_comparison_table(dplyr::select(cluster_comparison(), -significance))
  })
  
  output$cluster_comparison_heatmap_panel <- shiny::renderUI({
    if (!isTRUE(clustering_ready())) {
      return(NULL)
    }
    trajectory_plot_output("cluster_comparison_heatmap")
  })
  
  output$cluster_comparison_heatmap <- plotly::renderPlotly({
    results <- cluster_comparison()
    features <- results$feature[results$adjusted_p_value < 0.05]
    shiny::validate(shiny::need(length(features) > 0, "No features meet the FDR < 0.05 threshold for this comparison."))
    app_heatmap_plot(
      cellpacman::summarizeClusterFeatures(dimred_results()$data_norm, selected_cluster_assignments(), features),
      "cluster", "Cluster", "Median\nnormalized value"
    )
  })
  
  output$cluster_comparison_volcano_panel <- shiny::renderUI({
    if (!isTRUE(clustering_ready())) {
      return(NULL)
    }
    trajectory_plot_output("cluster_comparison_volcano")
  })
  
  output$cluster_comparison_volcano <- plotly::renderPlotly({
    results <- cluster_comparison()
    shiny::validate(shiny::need(nrow(results) > 0, "No comparable feature values are available."))
    app_volcano_plot(
      results,
      "median_difference",
      "adjusted_p_value",
      paste("Cluster", input$second_cluster, "versus Cluster", input$first_cluster)
    )
  })
  
  trajectory_compound_choices <- shiny::reactive({
    shiny::req(trajectory_results())
    
    results <- trajectory_results()
    if (!is.null(results$error)) {
      return(character(0))
    }
    
    compounds <- names(results$trajectory_curve_result$curves)
    stats::setNames(compounds, compounds)
  })
  
  shiny::observeEvent(trajectory_results(), {
    choices <- trajectory_compound_choices()
    selected <- if (length(choices) > 0) unname(choices[[1]]) else character(0)
    
    shiny::updateSelectInput(
      session,
      "trajectory_compound",
      choices = choices,
      selected = selected
    )
  })
  
  trajectory_plot_output <- app_plot_output
  
  output$trajectory_controls <- shiny::renderUI({
    if (!isTRUE(trajectory_ready())) {
      return(NULL)
    }
    
    choices <- trajectory_compound_choices()
    results <- trajectory_results()
    
    if (!is.null(results$error)) {
      return(
        shiny::div(
          style = "margin-top: 1rem; color: #666;",
          shiny::p("Curve fitting failed.")
        )
      )
    }
    
    if (length(choices) == 0) {
      return(
        shiny::div(
          style = "margin-top: 1rem; color: #666;",
          shiny::p("No fitted concentration trajectory curves are available.")
        )
      )
    }
    
    shiny::selectInput(
      "trajectory_compound",
      "Compound",
      choices = choices,
      selected = unname(choices[[1]])
    )
  })
  
  output$trajectory_curve_plot <- plotly::renderPlotly({
    results <- trajectory_results()
    shiny::validate(
      shiny::need(is.null(results$error), paste("Curve fitting failed:", results$error)),
      shiny::need(length(trajectory_compound_choices()) > 0,
                  "No concentration trajectory curves fitted. No compounds had enough clustered points for principal curve fitting.")
    )
    shiny::req(input$trajectory_compound)
    cellpacman::plotKnownCurvesPlotly(
      results$trajectory_curve_result,
      input$trajectory_compound,
      plot_limits = dimred_results()$plot_limits
    )
  })
  
  output$trajectory_projection_plot <- plotly::renderPlotly({
    results <- trajectory_results()
    shiny::validate(
      shiny::need(is.null(results$error), "Unknown projections unavailable. Curve fitting did not complete."),
      shiny::need(length(trajectory_compound_choices()) > 0,
                  "Unknown projections unavailable. No fitted concentration trajectory curves are available."),
      shiny::need(!is.null(results$projections),
                  "No unknown projections. No unknown compound clusters were available for projection.")
    )
    shiny::req(input$trajectory_compound)
    
    cellpacman::plotUnknownProjectionsPlotly(
      results$trajectory_curve_result,
      results$projections,
      input$trajectory_compound,
      plot_limits = dimred_results()$plot_limits
    )
  })
  
  trajectory_feature_result <- shiny::reactive({
    shiny::req(trajectory_results(), input$trajectory_compound)
    results <- trajectory_results()
    feature_error <- if (is.null(results$feature_error)) "Feature analysis failed." else results$feature_error
    shiny::validate(
      shiny::need(is.null(results$feature_error), feature_error),
      shiny::need(!is.null(results$feature_results), "Feature analysis is unavailable."),
      shiny::need(input$trajectory_compound %in% names(results$feature_results), "No feature analysis is available for the selected compound.")
    )
    results$feature_results[[input$trajectory_compound]]
  })
  
  output$trajectory_feature_table_panel <- shiny::renderUI({
    if (!isTRUE(trajectory_ready())) {
      return(NULL)
    }
    shinycssloaders::withSpinner(
      reactable::reactableOutput("trajectory_feature_table", height = "100%"),
      type = 1,
      color = "#18344d",
      proxy.height = "24rem"
    )
  })
  
  output$trajectory_feature_table <- reactable::renderReactable({
    results <- trajectory_feature_result()
    feature_comparison_table(dplyr::select(results$stats, -significance))
  })
  
  output$trajectory_feature_heatmap_panel <- shiny::renderUI({
    if (!isTRUE(trajectory_ready())) {
      return(NULL)
    }
    trajectory_plot_output("trajectory_feature_heatmap")
  })
  
  output$trajectory_feature_heatmap <- plotly::renderPlotly({
    results <- trajectory_feature_result()
    heatmap_data <- results$heatmap_data
    feature_cols <- setdiff(names(heatmap_data), c("bin", "mean_concentration"))
    shiny::validate(shiny::need(length(feature_cols) > 0, "No feature values are available for the heatmap."))
    
    plot_data <- heatmap_data |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(feature_cols),
        names_to = "feature",
        values_to = "value",
        cols_vary = "slowest"
      ) |>
      dplyr::transmute(
        feature = feature,
        bin = bin,
        value = value,
        concentration = mean_concentration
      )
    
    heatmap_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = bin,
        y = feature,
        fill = value,
        text = paste(feature, "Bin", bin, round(value, 3))
      )
    ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "green", mid = "white", high = "red", midpoint = 0) +
      ggplot2::labs(x = "Trajectory bin", y = "Feature", fill = "Mean\nnormalized value") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
    
    plotly::ggplotly(heatmap_plot, tooltip = "text") |>
      plotly::config(responsive = TRUE)
  })
  
  output$trajectory_feature_volcano_panel <- shiny::renderUI({
    if (!isTRUE(trajectory_ready())) {
      return(NULL)
    }
    trajectory_plot_output("trajectory_feature_volcano")
  })
  
  output$trajectory_feature_volcano <- plotly::renderPlotly({
    results <- trajectory_feature_result()
    plot_data <- dplyr::mutate(
      results$stats,
      significance = dplyr::recode(significance, "Not Significant" = "Not significant")
    )
    shiny::validate(shiny::need(nrow(plot_data) > 0, "No comparable feature values are available."))
    app_volcano_plot(plot_data, "median_difference", "adjusted_p_value", input$trajectory_compound)
  })
}
