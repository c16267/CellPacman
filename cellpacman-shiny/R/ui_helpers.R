app_plot_output <- function(id) {
  shiny::div(
    class = "plot-fill",
    style = "height: 100%; min-height: 0; overflow: hidden;",
    shinycssloaders::withSpinner(
      plotly::plotlyOutput(id, width = "100%", height = "100%"),
      type = 1, color = "#18344d", proxy.height = "100%"
    )
  )
}

app_upload_output <- function(id, label, use_example_data) {
  condition <- if (use_example_data) "true" else paste0("input.", id, "_file != null")
  shiny::tagList(
    shiny::conditionalPanel(
      condition = paste0("!(", condition, ")"),
      shiny::div(
        style = "text-align: center; padding: 50px; color: #666;",
        shiny::h4(paste("No", label, "table uploaded yet")),
        shiny::p(paste("Please upload a", label, "table to view it here."))
      )
    ),
    shiny::conditionalPanel(
      condition = condition,
      shinycssloaders::withSpinner(
        reactable::reactableOutput(paste0(id, "_table"), height = "100%"),
        type = 1, color = "#18344d", proxy.height = "24rem", caption = "loading data..."
      )
    )
  )
}
