# User interface: one bslib navbar page with four tabs (Data, Dimension
# Reduction, Clustering, Trajectory). Each analysis tab is enabled by the
# server only once the stage it depends on has completed.

# Static assets (www/cellpacman.css and www/cellpacman.js) are attached as an
# HTML dependency rather than through Shiny's implicit www/ serving, so that
# any host of app_ui() -- including the test app under tests/testthat/apps --
# serves the same files.
app_html_dependency <- function() {
  htmltools::htmlDependency(
    name = "cellpacman-app",
    version = as.character(utils::packageVersion("cellpacman")),
    src = c(file = app_www_dir()),
    stylesheet = "cellpacman.css",
    script = "cellpacman.js"
  )
}

app_ui <- function(use_example_data = FALSE){
  bslib::page_navbar(
    title = "PACMAN (PAinted Cell and coMpound ANalysis)",
    id = "main_nav",
    theme = bslib::bs_theme(version = 5),
    navbar_options = bslib::navbar_options(collapsible = FALSE),
    header = app_html_dependency(),

    # -------------------- DATA TAB --------------------
    bslib::nav_panel(
      "Data",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          shiny::fileInput("feature_file", "Feature table", accept = c(".csv", ".tsv", ".txt", ".xlsx")),
          shiny::fileInput("metadata_file", "Metadata table", accept = c(".csv", ".tsv", ".txt", ".xlsx")),
          shiny::textOutput("load_status"),
          if (use_example_data) {
            shiny::tags$p(
              class = "text-muted",
              "Example feature and metadata tables are loaded. Upload files to replace either table."
            )
          },
          width = 300
        ),
        
          shiny::div(
            class = "app-panel-grid",
          
          bslib::card(
            bslib::card_header("Features"),
            app_upload_output("feature", "feature", use_example_data),
            fill = TRUE
          ),
          bslib::card(
            bslib::card_header("Metadata"),
            app_upload_output("metadata", "metadata", use_example_data),
            fill = TRUE
          )
        )
      ),
      value = "data"
    ),

    # -------------------- DIMENSION REDUCTION TAB --------------------
    bslib::nav_panel(
      "Dimension Reduction",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          bslib::input_task_button("run_dimred_button", "Run Dimension Reduction", type = "default", auto_reset = FALSE),
          width = 300
        ),
        
        shiny::div(
          class = "app-panel-grid",
          bslib::layout_columns(
            col_widths = c(6, 6),
            class = "app-panel-row",
            bslib::card(
              bslib::card_header("UMAP Distribution: Known and Unknown Compounds"),
               app_plot_output("umap_plot"),
              fill = TRUE,
              class = "plot-card"
            ),
          bslib::card(
            bslib::card_header("Compare Compounds"),
            shiny::uiOutput("compound_comparison_controls"),
            fill = TRUE
          ),
          ),
          bslib::layout_columns(
            col_widths = c(4, 4, 4),
            class = "app-panel-row",
          bslib::card(
            bslib::card_header("Feature Comparison"),
            shiny::uiOutput("compound_comparison_table_panel"),
            fill = TRUE
          ),
          bslib::card(
            bslib::card_header("Significant Feature Heatmap"),
            shiny::uiOutput("compound_comparison_heatmap_panel"),
            fill = TRUE,
            class = "plot-card"
          ),
          bslib::card(
            bslib::card_header("Volcano Plot"),
            shiny::uiOutput("compound_comparison_volcano_panel"),
            fill = TRUE,
            class = "plot-card"
          )
          )
        )
      ),
      value = "dimension_reduction"
    ),

    # -------------------- CLUSTERING TAB --------------------
    bslib::nav_panel(
      "Clustering",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          bslib::input_task_button("run_clustering_button", "Run Clustering", type = "default", auto_reset = FALSE),
          shiny::tags$details(
            class = "app-sidebar-parameters",
            shiny::tags$summary("Clustering Parameters"),
            shiny::div(
              class = "app-sidebar-parameters-content",
              shiny::selectInput(
                "cluster_method",
                "Method",
                choices = c("DBSCAN" = "dbscan"),
                selected = "dbscan"
              ),
              shiny::numericInput(
                "cluster_eps",
                "Epsilon",
                value = 0.5,
                min = 0,
                step = 0.1
              ),
              shiny::numericInput(
                "cluster_min_pts",
                "Minimum Points",
                value = 5,
                min = 1,
                step = 1
              )
            )
          ),
          width = 300
        ),
        
        shiny::div(
          class = "app-panel-grid",
          bslib::layout_columns(
            col_widths = c(6, 6),
            class = "app-panel-row",
             bslib::card(
               bslib::card_header("DBSCAN Clustering"),
                app_plot_output("cluster_plot"),
               fill = TRUE,
               class = "plot-card"
             ),
            bslib::card(
              bslib::card_header("Compare Clusters"),
              shiny::uiOutput("cluster_comparison_controls"),
              fill = TRUE
            )
          ),
          bslib::layout_columns(
            col_widths = c(4, 4, 4),
            class = "app-panel-row",
            bslib::card(
              bslib::card_header("Feature Comparison"),
              shiny::uiOutput("cluster_comparison_table_panel"),
              fill = TRUE
            ),
            bslib::card(
              bslib::card_header("Significant Feature Heatmap"),
              shiny::uiOutput("cluster_comparison_heatmap_panel"),
              fill = TRUE,
              class = "plot-card"
            ),
            bslib::card(
              bslib::card_header("Volcano Plot"),
              shiny::uiOutput("cluster_comparison_volcano_panel"),
              fill = TRUE,
              class = "plot-card"
            )
          )
        )
      ),
      value = "clustering"
    ),
    
    # -------------------- TRAJECTORY TAB --------------------
    bslib::nav_panel(
      "Trajectory",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          bslib::input_task_button("run_trajectory_button", "Run Trajectory Analysis", type = "default", auto_reset = FALSE),
          shiny::uiOutput("trajectory_controls"),
          shiny::tags$details(
            class = "app-sidebar-parameters",
            shiny::tags$summary("Curve Parameters"),
            shiny::div(
              class = "app-sidebar-parameters-content",
              shiny::selectInput(
                "curve_method",
                "Method",
                choices = c("Principal Curve" = "princurve"),
                selected = "princurve"
              ),
              shiny::selectInput(
                "curve_smoother",
                "Smoother",
                choices = c("Lowess" = "lowess", "Spline" = "smooth_spline"),
                selected = "lowess"
              )
            )
          ),
          width = 300
        ),
        
         shiny::div(
           class = "app-panel-grid",
           bslib::layout_columns(
             col_widths = c(6, 6),
             class = "app-panel-row",
             bslib::card(
                bslib::card_header("Concentration Trajectory Curve"),
                app_plot_output("trajectory_curve_plot"),
               fill = TRUE,
               class = "plot-card"
             ),
             bslib::card(
               bslib::card_header("Unknown Cluster Projections"),
                app_plot_output("trajectory_projection_plot"),
               fill = TRUE,
               class = "plot-card"
             )
           ),
           bslib::layout_columns(
             col_widths = c(4, 4, 4),
             class = "app-panel-row",
             bslib::card(
               bslib::card_header("Feature Comparison"),
               shiny::uiOutput("trajectory_feature_table_panel"),
               fill = TRUE
             ),
             bslib::card(
               bslib::card_header("Significant Feature Heatmap"),
               shiny::uiOutput("trajectory_feature_heatmap_panel"),
               fill = TRUE,
               class = "plot-card"
             ),
             bslib::card(
               bslib::card_header("Volcano Plot"),
               shiny::uiOutput("trajectory_feature_volcano_panel"),
               fill = TRUE,
               class = "plot-card"
             )
           ),
           fill = TRUE
         )
      ),
      value = "trajectory"
    )
  )
}
