#' Plot reduced Cell Painting coordinates
#'
#' @param dimR A `cellpainting_dimred` object.
#' @param compound_name_map Optional named character vector used to recode compound labels.
#'
#' @return A `ggplot` object.
#' @export
plotDimred <- function(
  dimR,
  compound_name_map = NULL
) {
  if (!inherits(dimR, "cellpainting_dimred")) {
    stop("`dimR` must be a `cellpainting_dimred` object created by `dimReduce()`.")
  }

  plot_tbl <- tibble::as_tibble(dimR$reduced_data)
  validate_xy_cols(plot_tbl)
  validate_dimred_cols(plot_tbl, "Compound")

  plot_tbl$is_unknown <- is_unknown_value(plot_tbl$Compound)
  plot_tbl$compound_label <- recode_compound_labels(plot_tbl$Compound, compound_name_map)

  unknown_tbl <- plot_tbl[plot_tbl$is_unknown, , drop = FALSE]
  known_tbl <- plot_tbl[!plot_tbl$is_unknown, , drop = FALSE]

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = unknown_tbl,
      ggplot2::aes(x = x, y = y),
      color = "darkgray",
      alpha = 1,
      size = 1,
      shape = 19
    ) +
    ggplot2::geom_point(
      data = known_tbl,
      ggplot2::aes(x = x, y = y, color = compound_label),
      alpha = 1,
      size = 1,
      shape = 19
    ) +
    ggplot2::scale_color_viridis_d(option = "viridis") +
    ggplot2::labs(
      title = "UMAP Distribution: Known and Unknown Compounds",
      subtitle = "Unknown compound points shown in gray.",
      x = "DIM1",
      y = "DIM2",
      color = "Compound"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "right")
}

#' Plot DBSCAN clusters for unknown compounds
#'
#' @param predClust A `cellpainting_clusters` object.
#'
#' @return A `ggplot` object.
#' @export
plotUnknownClusters <- function(predClust) {
  if (!inherits(predClust, "cellpainting_clusters")) {
    stop("`predClust` must be a `cellpainting_clusters` object created by `compoundCluster()`.")
  }

  if (is.null(predClust$unknown_assignment) || nrow(predClust$unknown_assignment) == 0) {
    stop("`predClust` does not contain clustered unknown compounds.")
  }

  unknown_assignment <- tibble::as_tibble(predClust$unknown_assignment)
  validate_xy_cols(unknown_assignment)
  validate_dimred_cols(unknown_assignment, "cluster")

  unknown_centroids <- tibble::as_tibble(predClust$unknown_centroids)

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = unknown_assignment,
      ggplot2::aes(x = x, y = y, color = cluster),
      shape = 19,
      size = 0.5
    ) +
    ggplot2::geom_point(
      data = unknown_centroids,
      ggplot2::aes(x = mean_x, y = mean_y),
      color = "black",
      fill = "gold",
      size = 5,
      shape = 23
    ) +
    geom_text_label(
      data = unknown_centroids,
      ggplot2::aes(x = mean_x, y = mean_y, label = cluster),
      vjust = -1,
      color = "black",
      fontface = "bold",
      size = 6
    ) +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::labs(
      title = "DBSCAN Clustering of Unknown Compound Points",
      subtitle = "Yellow diamonds represent cluster centroids. Cluster 0 is noise.",
      x = "DIM1",
      y = "DIM2",
      color = "Cluster"
    ) +
    ggplot2::theme_bw(base_size = 12)
}

#' Plot DBSCAN clusters for one known compound
#'
#' @param predClust A `cellpainting_clusters` object.
#' @param compound Compound name to plot.
#'
#' @return A `ggplot` object.
#' @export
plotKnownClusters <- function(predClust, compound) {
  if (!inherits(predClust, "cellpainting_clusters")) {
    stop("`predClust` must be a `cellpainting_clusters` object created by `compoundCluster()`.")
  }

  compound <- validate_plot_compound(compound, names(predClust$known_assignments))
  sub_points <- tibble::as_tibble(predClust$known_assignments[[compound]])
  summary_tbl <- tibble::as_tibble(predClust$known_clusters[[compound]])
  if (nrow(summary_tbl) == 0) {
    stop("No known cluster summary is available for `compound = \"", compound, "\"`.")
  }

  validate_xy_cols(sub_points)
  validate_dimred_cols(sub_points, "Concentration")
  sub_points$concentration_factor <- as.factor(sub_points$Concentration)

  unique_concs <- sort(unique(sub_points$Concentration))
  discrete_colors <- viridisLite::viridis(length(unique_concs), option = "viridis")
  names(discrete_colors) <- as.character(unique_concs)

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = sub_points,
      ggplot2::aes(x = x, y = y, color = concentration_factor),
      alpha = 0.7,
      size = 1.4
    ) +
    ggplot2::geom_point(
      data = summary_tbl,
      ggplot2::aes(x = mean_x, y = mean_y),
      shape = 23,
      fill = "gold",
      color = "black",
      size = 4,
      stroke = 1.5
    ) +
    geom_text_label(
      data = summary_tbl,
      ggplot2::aes(x = mean_x, y = mean_y, label = format_conc(mean_concentration)),
      color = "black",
      fontface = "bold",
      size = 6,
      bg.color = "white",
      bg.r = 0.15,
      box.padding = 0.5,
      point.padding = 0.3
    ) +
    ggplot2::scale_color_manual(
      values = discrete_colors,
      labels = format_conc(names(discrete_colors)),
      name = "Concentration (M)"
    ) +
    ggplot2::labs(
      title = paste("DBSCAN clustering:", compound),
      subtitle = "Diamonds represent cluster centroids labeled by mean concentration.",
      x = "DIM1",
      y = "DIM2"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

#' Plot a concentration trajectory curve for one known compound
#'
#' @param predCurve A `cellpainting_curves` object.
#' @param compound Compound name to plot.
#'
#' @return A `ggplot` object.
#' @export
plotKnownCurves <- function(predCurve, compound) {
  if (!inherits(predCurve, "cellpainting_curves")) {
    stop("`predCurve` must be a `cellpainting_curves` object.")
  }

  compound <- validate_plot_compound(compound, names(predCurve$curves))
  plot_data <- prepare_known_curve_plot_data(predCurve, compound)
  sub_points <- plot_data$points
  summary_tbl <- plot_data$summary
  trajectory_curve_data <- plot_data$curve

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = sub_points,
      ggplot2::aes(x = x, y = y, color = Concentration),
      alpha = 0.15,
      size = 0.8
    ) +
    ggplot2::geom_segment(
      data = summary_tbl,
      ggplot2::aes(x = mean_x, y = mean_y, xend = proj_x, yend = proj_y),
      linetype = "dashed",
      color = "darkred",
      linewidth = 0.6
    ) +
    ggplot2::geom_path(
      data = trajectory_curve_data,
      ggplot2::aes(x = x, y = y, color = estimated_concentration),
      linewidth = 1.8,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      data = summary_tbl,
      ggplot2::aes(x = mean_x, y = mean_y),
      shape = 21,
      color = "black",
      fill = "gold",
      size = 4,
      stroke = 1
    ) +
    geom_text_label(
      data = summary_tbl,
      ggplot2::aes(x = proj_x, y = proj_y, label = format_conc(mean_concentration)),
      color = "black",
      fontface = "bold",
      size = 5.5,
      nudge_x = -0.8,
      nudge_y = 0.8,
      bg.color = "white",
      bg.r = 0.1
    ) +
    ggplot2::scale_color_viridis_c(
      option = "viridis",
      name = "Concentration (M)",
      labels = format_conc,
      guide = ggplot2::guide_colorbar()
    ) +
    ggplot2::labs(
      title = paste("Concentration Trajectory Curve:", compound),
      subtitle = "Labels show centroid concentrations mapped onto the curve.",
      x = "DIM1",
      y = "DIM2"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

#' Plot a known-compound curve as an interactive Plotly object
#'
#' @param predCurve A `cellpainting_curves` object.
#' @param compound Compound name to plot.
#' @param plot_limits Optional list with numeric `x` and `y` limits.
#'
#' @return A Plotly htmlwidget.
#' @export
plotKnownCurvesPlotly <- function(predCurve, compound, plot_limits = NULL) {
  if (!inherits(predCurve, "cellpainting_curves")) {
    stop("`predCurve` must be a `cellpainting_curves` object.")
  }

  compound <- validate_plot_compound(compound, names(predCurve$curves))
  plot_data <- prepare_known_curve_plot_data(predCurve, compound)
  sub_points <- plot_data$points
  summary_tbl <- plot_data$summary
  trajectory_curve_data <- plot_data$curve
  point_x <- sub_points$x
  point_y <- sub_points$y
  point_text <- paste0(
    "Concentration: ", format_conc(sub_points$Concentration),
    "<br>x: ", signif(point_x, 5), "<br>y: ", signif(point_y, 5)
  )
  centroid_text <- paste0(
    "Mean concentration: ", format_conc(summary_tbl$mean_concentration),
    "<br>x: ", signif(summary_tbl$mean_x, 5), "<br>y: ", signif(summary_tbl$mean_y, 5)
  )

  plot <- plotly::plot_ly()
  plot <- plotly::add_trace(
    plot,
    x = point_x,
    y = point_y,
    type = "scatter",
    mode = "markers",
    name = paste(compound, "wells"),
    text = point_text,
    hovertemplate = "%{text}<extra></extra>",
    marker = list(color = "red", opacity = 0.3, size = 5)
  )
  plot <- add_plotly_concentration_curve(plot, trajectory_curve_data)
  plot <- add_plotly_projection_segments(
    plot,
    summary_tbl$mean_x,
    summary_tbl$mean_y,
    summary_tbl$proj_x,
    summary_tbl$proj_y
  )
  plot <- plotly::add_trace(
    plot,
    x = summary_tbl$mean_x,
    y = summary_tbl$mean_y,
    type = "scatter",
    mode = "markers+text",
    name = "Known centroids",
    showlegend = FALSE,
    text = format_conc(summary_tbl$mean_concentration),
    textposition = "top center",
    textfont = list(color = "black", size = 12, family = "Arial Black"),
    hovertext = centroid_text,
    hovertemplate = "%{hovertext}<extra></extra>",
    marker = list(symbol = "diamond", color = "gold", size = 12, line = list(color = "black", width = 1.2))
  )

  plotly_trajectory_layout(
    plot,
    title = paste("Concentration Trajectory Curve:", compound),
    subtitle = "Labels show centroid concentrations mapped onto the curve.",
    plot_limits = plot_limits
  )
}

#' Plot unknown-cluster projections onto one known compound curve
#'
#' @param predCurve A `cellpainting_curves` object.
#' @param projections Projection results returned by `curveProject()`.
#' @param compound Compound name to plot.
#'
#' @return A `ggplot` object.
#' @export
plotUnknownProjections <- function(predCurve, projections, compound) {
  plot_data <- prepare_unknown_projection_plot_data(predCurve, projections, compound)
  proj_tbl <- plot_data$projections
  trajectory_curve_data <- plot_data$curve
  compound_raw_data <- plot_data$compound_points
  unknown_assignment <- plot_data$unknown_points

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = compound_raw_data,
      ggplot2::aes(x = x, y = y),
      color = "red",
      alpha = 0.3,
      size = 0.8
    ) +
    ggplot2::geom_point(
      data = unknown_assignment,
      ggplot2::aes(x = x, y = y),
      color = "darkgrey",
      alpha = 0.2,
      size = 0.5
    ) +
    ggplot2::geom_path(
      data = trajectory_curve_data,
      ggplot2::aes(x = x, y = y, color = estimated_concentration),
      linewidth = 1.6,
      alpha = 0.9,
      lineend = "round"
    ) +
    ggplot2::geom_segment(
      data = proj_tbl,
      ggplot2::aes(x = unknown_x, y = unknown_y, xend = projected_x, yend = projected_y),
      linetype = "dashed",
      color = "steelblue",
      linewidth = 0.5,
      alpha = 0.8
    ) +
    ggplot2::geom_point(
      data = proj_tbl,
      ggplot2::aes(x = unknown_x, y = unknown_y),
      shape = 21,
      fill = "gold",
      color = "black",
      size = 3.5,
      stroke = 1.2
    ) +
    geom_text_label(
      data = proj_tbl,
      ggplot2::aes(
        x = unknown_x,
        y = unknown_y,
        label = paste0("Cluster ", unknown_cluster, ", Est. Conc. ", format_conc(estimated_concentration))
      ),
      color = "black",
      fontface = "bold",
      size = 5.5,
      bg.color = "white",
      bg.r = 0.1,
      box.padding = 0.6,
      point.padding = 0.3
    ) +
    ggplot2::scale_color_viridis_c(option = "viridis", name = "Concentration (M)", labels = format_conc) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = paste("Concentration Trajectory Curve & Unknown Cluster Projections:", plot_data$compound),
      subtitle = paste("Gold points = Unknown centroids | Gray points = Unknown | Red points =", plot_data$compound),
      x = "DIM1",
      y = "DIM2"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

#' Plot unknown-cluster projections as a native Plotly widget
#'
#' @param predCurve A `cellpainting_curves` object.
#' @param projections Projection results returned by `curveProject()`.
#' @param compound Compound name to plot.
#' @param plot_limits Optional list with numeric `x` and `y` limits.
#'
#' @return A Plotly htmlwidget.
#' @export
plotUnknownProjectionsPlotly <- function(predCurve, projections, compound, plot_limits = NULL) {
  plot_data <- prepare_unknown_projection_plot_data(predCurve, projections, compound)
  proj_tbl <- plot_data$projections
  trajectory_curve_data <- plot_data$curve
  compound_points <- plot_data$compound_points
  unknown_points <- plot_data$unknown_points
  compound_x <- compound_points$x
  compound_y <- compound_points$y
  unknown_x <- unknown_points$x
  unknown_y <- unknown_points$y
  projection_text <- paste0(
    "Cluster ", proj_tbl$unknown_cluster,
    "<br>Estimated concentration: ", format_conc(proj_tbl$estimated_concentration),
    "<br>x: ", signif(proj_tbl$unknown_x, 5), "<br>y: ", signif(proj_tbl$unknown_y, 5)
  )
  plot <- plotly::plot_ly()
  plot <- plotly::add_trace(
    plot,
    x = compound_x,
    y = compound_y,
    type = "scatter",
    mode = "markers",
    name = paste(plot_data$compound, "wells"),
    showlegend = FALSE,
    text = paste0("x: ", signif(compound_x, 5), "<br>y: ", signif(compound_y, 5)),
    hovertemplate = "%{text}<extra></extra>",
    marker = list(color = "red", opacity = 0.3, size = 5)
  )
  plot <- plotly::add_trace(
    plot,
    x = unknown_x,
    y = unknown_y,
    type = "scatter",
    mode = "markers",
    name = "Unknown wells",
    showlegend = FALSE,
    text = paste0("x: ", signif(unknown_x, 5), "<br>y: ", signif(unknown_y, 5)),
    hovertemplate = "%{text}<extra></extra>",
    marker = list(color = "darkgrey", opacity = 0.2, size = 4)
  )

  plot <- add_plotly_concentration_curve(plot, trajectory_curve_data)
  plot <- add_plotly_projection_segments(
    plot,
    proj_tbl$unknown_x,
    proj_tbl$unknown_y,
    proj_tbl$projected_x,
    proj_tbl$projected_y
  )
  plot <- plotly::add_trace(
    plot,
    x = proj_tbl$unknown_x,
    y = proj_tbl$unknown_y,
    type = "scatter",
    mode = "markers+text",
    name = "Unknown centroids",
    showlegend = FALSE,
    text = paste0("Cluster ", proj_tbl$unknown_cluster, ", Est. Conc. ", format_conc(proj_tbl$estimated_concentration)),
    textposition = "top center",
    textfont = list(color = "black", size = 12, family = "Arial Black"),
    hovertext = projection_text,
    hovertemplate = "%{hovertext}<extra></extra>",
    marker = list(symbol = "diamond", color = "gold", size = 12, line = list(color = "black", width = 1.2))
  )

  plotly_trajectory_layout(
    plot,
    title = paste("Concentration Trajectory Curve & Unknown Cluster Projections:", plot_data$compound),
    subtitle = paste("Gold points = Unknown centroids | Gray points = Unknown | Red points =", plot_data$compound),
    plot_limits = plot_limits
  )
}

add_plotly_concentration_curve <- function(plot, trajectory_curve_data) {
  concentration_curve_segments <- trajectory_curve_data[-nrow(trajectory_curve_data), , drop = FALSE]
  concentration_curve_segments$xend <- trajectory_curve_data$x[-1]
  concentration_curve_segments$yend <- trajectory_curve_data$y[-1]
  concentration_curve_segments$segment_concentration <- (
    concentration_curve_segments$estimated_concentration + trajectory_curve_data$estimated_concentration[-1]
  ) / 2
  color_range <- range(concentration_curve_segments$segment_concentration, finite = TRUE)
  palette <- viridisLite::viridis(256, option = "viridis")
  color_index <- if (diff(color_range) == 0) {
    rep(1L, nrow(concentration_curve_segments))
  } else {
    pmax(1L, pmin(256L, floor((concentration_curve_segments$segment_concentration - color_range[1]) / diff(color_range) * 255) + 1L))
  }

  for (i in seq_len(nrow(concentration_curve_segments))) {
    plot <- plotly::add_trace(
      plot,
      x = c(concentration_curve_segments$x[i], concentration_curve_segments$xend[i]),
      y = c(concentration_curve_segments$y[i], concentration_curve_segments$yend[i]),
      type = "scatter",
      mode = "lines",
      showlegend = FALSE,
      text = paste0("Concentration: ", format_conc(concentration_curve_segments$segment_concentration[i])),
      hovertemplate = "%{text}<extra></extra>",
      line = list(color = palette[color_index[i]], width = 2.5)
    )
  }

  plotly::add_trace(
    plot,
    x = concentration_curve_segments$x,
    y = concentration_curve_segments$y,
    type = "scatter",
    mode = "markers",
    showlegend = FALSE,
    hoverinfo = "skip",
    marker = list(
      color = concentration_curve_segments$segment_concentration,
      cmin = color_range[1],
      cmax = color_range[2],
      colorscale = "Viridis",
      showscale = TRUE,
      size = 0,
      opacity = 0,
      colorbar = list(title = list(text = "Concentration (M)"), tickformat = ".2e")
    )
  )
}

add_plotly_projection_segments <- function(plot, x, y, xend, yend) {
  plotly::add_trace(
    plot,
    x = as.vector(rbind(x, xend, NA_real_)),
    y = as.vector(rbind(y, yend, NA_real_)),
    type = "scatter",
    mode = "lines",
    showlegend = FALSE,
    hoverinfo = "skip",
    line = list(color = "steelblue", dash = "dash", width = 1)
  )
}

plotly_trajectory_layout <- function(plot, title, subtitle, plot_limits) {
  layout_args <- list(
    title = list(text = paste0(title, "<br><sup>", subtitle, "</sup>")),
    xaxis = list(title = "DIM1"),
    yaxis = list(title = "DIM2", scaleanchor = "x", scaleratio = 1),
    legend = list(orientation = "v"),
    font = list(size = 12)
  )
  if (!is.null(plot_limits)) {
    layout_args$xaxis$range <- plot_limits$x
    layout_args$yaxis$range <- plot_limits$y
  }

  do.call(plotly::layout, c(list(plot), layout_args)) %>%
    plotly::config(responsive = TRUE)
}

prepare_unknown_projection_plot_data <- function(pred_curve, projections, compound) {
  if (!inherits(pred_curve, "cellpainting_curves")) {
    stop("`predCurve` must be a `cellpainting_curves` object.")
  }
  validate_projection_results(projections)

  compound <- validate_plot_compound(compound, names(pred_curve$curves))
  proj_tbl <- projections$projections[projections$projections$compound == compound, , drop = FALSE]
  if (nrow(proj_tbl) == 0) {
    stop("No unknown projection results are available for `compound = \"", compound, "\"`.")
  }

  trajectory_curve_data <- add_estimated_concentration(pred_curve$curves[[compound]])
  compound_raw_data <- tibble::as_tibble(pred_curve$point_assignments[[compound]])
  unknown_assignment <- tibble::as_tibble(projections$unknown_assignment)
  validate_xy_cols(compound_raw_data)
  validate_xy_cols(unknown_assignment)

  list(
    compound = compound,
    projections = proj_tbl,
    curve = trajectory_curve_data,
    compound_points = compound_raw_data,
    unknown_points = unknown_assignment
  )
}

format_conc <- function(x) {
  vapply(x, function(val) {
    if (is.na(val)) {
      return(NA_character_)
    }

    numeric_val <- suppressWarnings(as.numeric(val))
    if (is.na(numeric_val)) {
      return(as.character(val))
    }

    if (abs(numeric_val) > 0 && abs(numeric_val) < 0.01) {
      formatC(numeric_val, format = "e", digits = 2)
    } else {
      formatC(numeric_val, format = "f", digits = 2)
    }
  }, character(1))
}

geom_text_label <- function(...) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    return(ggrepel::geom_text_repel(...))
  }

  args <- list(...)
  args[c("bg.color", "bg.colour", "bg.r", "box.padding", "point.padding")] <- NULL
  do.call(ggplot2::geom_text, args)
}

recode_compound_labels <- function(x, compound_name_map = NULL) {
  labels <- as.character(x)
  if (is.null(compound_name_map)) {
    return(labels)
  }

  if (is.null(names(compound_name_map)) || any(names(compound_name_map) == "")) {
    stop("`compound_name_map` must be a named character vector.")
  }

  mapped <- compound_name_map[labels]
  labels[!is.na(mapped)] <- unname(mapped[!is.na(mapped)])
  labels
}

validate_plot_compound <- function(compound, available) {
  if (!is.character(compound) || length(compound) != 1 || is.na(compound) || identical(compound, "")) {
    stop("`compound` must be a single compound name.")
  }

  if (!compound %in% available) {
    stop(
      "`compound = \"", compound, "\"` was not found. Available compounds: ",
      paste(available, collapse = ", ")
    )
  }

  compound
}

add_estimated_concentration <- function(trajectory_curve) {
  trajectory_curve_data <- tibble::as_tibble(trajectory_curve$curve_data)
  concentration_references <- dplyr::distinct(
    tibble::as_tibble(trajectory_curve$ref_points),
    known_arc_length_ref,
    .keep_all = TRUE
  )
  trajectory_curve_data <- trajectory_curve_data[finite_trajectory_curve_points(trajectory_curve_data), , drop = FALSE]
  concentration_references <- concentration_references[
    finite_trajectory_reference_points(concentration_references),
    ,
    drop = FALSE
  ]
  if (nrow(trajectory_curve_data) < 2 || nrow(concentration_references) < 2) {
    stop("Curve data and reference points must contain at least 2 finite points.")
  }

  trajectory_curve_data$estimated_concentration <- stats::approx(
    x = concentration_references$known_arc_length_ref,
    y = concentration_references$concentration_ref,
    xout = trajectory_curve_data$arc_length,
    rule = 2
  )$y

  trajectory_curve_data
}

prepare_known_curve_plot_data <- function(pred_curve, compound) {
  sub_points <- tibble::as_tibble(pred_curve$point_assignments[[compound]])
  validate_xy_cols(sub_points)
  validate_dimred_cols(sub_points, "Concentration")
  trajectory_curve <- pred_curve$curves[[compound]]
  list(
    points = sub_points,
    summary = tibble::as_tibble(trajectory_curve$summary),
    curve = add_estimated_concentration(trajectory_curve)
  )
}

validate_projection_results <- function(projections) {
  required_names <- c("projections", "unknown_assignment", "unknown_centroids")
  missing_names <- setdiff(required_names, names(projections))
  if (length(missing_names) > 0) {
    stop("`projections` must be a projection list returned by `curveProject()`.")
  }

  required_cols <- c(
    "unknown_cluster", "compound", "estimated_concentration",
    "projected_x", "projected_y", "unknown_x", "unknown_y"
  )
  validate_dimred_cols(tibble::as_tibble(projections$projections), required_cols)

  invisible(projections)
}
