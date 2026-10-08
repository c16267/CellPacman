# Server-side helpers: upload parsing, the analysis wrappers executed in the
# background workers, and the presentation-only plots (volcano, heatmap).
#
# Everything scientific is delegated to the cellpacman package through its
# exported API. Package calls are namespace-qualified because these functions
# are serialized to background R processes that have not attached the package.

read_uploaded_table <- function(upload) {
  ext <- tolower(tools::file_ext(upload$name))
  tryCatch(
    cellpacman::readCellPainting(upload$datapath, ext = ext),
    error = function(e) stop(paste("Error reading file:", e$message))
  )
}

# `example` is "features" or "metadata" (see cellpacman::exampleDataPath()).
app_read_table <- function(upload, example) {
  if (!is.null(upload)) return(read_uploaded_table(upload))
  cellpacman::readCellPainting(cellpacman::exampleDataPath(example))
}

app_cluster_analysis <- function(dimred, method, eps, min_pts) {
  list(clusters = cellpacman::compoundCluster(dimred, method = method, eps = eps, minPts = min_pts))
}

app_dimred_analysis <- function(data) {
  data_norm <- cellpacman::normalize(data)
  dimred <- cellpacman::dimReduce(
    data_norm,
    random_seed = digest::digest2int("Random Seed")
  )
  list(
    data_norm = data_norm, dimred = dimred,
    plot_limits = list(
      x = range(dimred$reduced_data$x, na.rm = TRUE, finite = TRUE),
      y = range(dimred$reduced_data$y, na.rm = TRUE, finite = TRUE)
    )
  )
}

app_trajectory_analysis <- function(dimred, data_norm, method, eps, min_pts, smoother) {
  trajectory_curve_result <- tryCatch(
    cellpacman::curveEstimate(dimred, method = method, eps = eps, minPts = min_pts, smoother = smoother),
    error = function(e) e
  )
  if (inherits(trajectory_curve_result, "error")) {
    return(list(
      trajectory_curve_result = NULL,
      projections = NULL,
      feature_results = NULL,
      feature_error = NULL,
      error = conditionMessage(trajectory_curve_result)
    ))
  }
  
  projections <- tryCatch(cellpacman::curveProject(trajectory_curve_result), error = function(e) NULL)
  feature_results <- tryCatch(cellpacman::selectFeatureCurve(trajectory_curve_result, data_norm), error = function(e) e)
  feature_error <- NULL
  if (inherits(feature_results, "error")) {
    feature_error <- conditionMessage(feature_results)
    feature_results <- NULL
  }
  list(
    trajectory_curve_result = trajectory_curve_result,
    projections = projections,
    feature_results = feature_results,
    feature_error = feature_error,
    error = NULL
  )
}

app_volcano_plot <- function(data, x_col, y_col, title) {
  plot_data <- dplyr::mutate(
    data,
    volcano_y = -log10(pmax(.data[[y_col]], .Machine$double.xmin))
  )
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data[[x_col]], y = volcano_y, color = significance, text = feature)
  ) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    ggplot2::scale_color_manual(values = c("Negative" = "blue", "Positive" = "red", "Pass p-value cutoff" = "grey", "Not significant" = "black", "Not Significant" = "black")) +
    ggplot2::labs(title = title, x = "Median difference", y = "-log10 adjusted p-value") +
    ggplot2::theme_minimal(base_size = 12)
  plotly::ggplotly(plot, tooltip = c("text", "x", "y")) |>
    plotly::config(responsive = TRUE)
}

app_heatmap_data <- function(data, group_col, ordered = TRUE) {
  matrix <- stats::xtabs(stats::reformulate(c("feature", group_col), response = "value"), data = data)
  matrix[!is.finite(matrix)] <- 0
  if (ordered) {
    row_order <- if (nrow(matrix) > 1) stats::hclust(stats::dist(matrix))$order else 1
    column_order <- if (ncol(matrix) > 1) stats::hclust(stats::dist(t(matrix)))$order else 1
    matrix <- matrix[row_order, column_order, drop = FALSE]
  }
  out <- as.data.frame(as.table(matrix))
  names(out) <- c("feature", group_col, "value")
  out
}

app_heatmap_plot <- function(data, group_col, group_label, value_label, ordered = TRUE, rotate_x = FALSE) {
  plot_data <- app_heatmap_data(data, group_col, ordered)
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data[[group_col]], y = feature, fill = value, text = paste(feature, .data[[group_col]], round(value, 3)))
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "green", mid = "white", high = "red", midpoint = 0) +
    ggplot2::labs(x = group_label, y = "Feature", fill = value_label) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = if (rotate_x) ggplot2::element_text(angle = 45, hjust = 1) else ggplot2::element_text(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  plotly::ggplotly(plot, tooltip = "text") |>
    plotly::config(responsive = TRUE)
}
