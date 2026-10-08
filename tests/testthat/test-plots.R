test_that("plot functions validate input classes", {
  expect_error(plotDimred(list()), "`dimR` must be")
  expect_error(plotUnknownClusters(list()), "`predClust` must be")
  expect_error(plotKnownClusters(list(), "A"), "`predClust` must be")
  expect_error(plotKnownCurves(list(), "A"), "`predCurve` must be")
  expect_error(plotKnownCurvesPlotly(list(), "A"), "`predCurve` must be")
  expect_error(plotUnknownProjections(list(), list(), "A"), "`predCurve` must be")
})

test_that("plot functions return expected plot objects", {
  dimred <- structure(
    list(
      reduced_data = tibble::tibble(
        WellId = paste0("W", 1:8),
        Compound = c("DrugA", "DrugA", "DrugA", "DrugA", NA, NA, NA, NA),
        Concentration = c(0.01, 0.01, 1, 1, NA, NA, NA, NA),
        x = c(0, 0.1, 1, 1.1, 2, 2.1, 3, 3.1),
        y = c(0, 0.1, 1, 1.1, 2, 2.1, 3, 3.1)
      ),
      params = list()
    ),
    class = "cellpainting_dimred"
  )

  clusters <- structure(
    list(
      known_assignments = list(
        DrugA = tibble::tibble(
          WellId = paste0("W", 1:4),
          Compound = "DrugA",
          Concentration = c(0.01, 0.01, 1, 1),
          x = c(0, 0.1, 1, 1.1),
          y = c(0, 0.1, 1, 1.1),
          cluster = as.factor(c(1, 1, 2, 2))
        )
      ),
      known_clusters = list(
        DrugA = tibble::tibble(
          cluster = as.factor(c(1, 2)),
          mean_x = c(0.05, 1.05),
          mean_y = c(0.05, 1.05),
          mean_concentration = c(0.01, 1),
          point_count = c(2, 2)
        )
      ),
      unknown_assignment = tibble::tibble(
        WellId = paste0("W", 5:8),
        Compound = NA_character_,
        x = c(2, 2.1, 3, 3.1),
        y = c(2, 2.1, 3, 3.1),
        cluster = as.factor(c(1, 1, 2, 2))
      ),
      unknown_centroids = tibble::tibble(
        cluster = as.factor(c(1, 2)),
        mean_x = c(2.05, 3.05),
        mean_y = c(2.05, 3.05),
        point_count = c(2, 2)
      ),
      params = list()
    ),
    class = "cellpainting_clusters"
  )

  curves <- structure(
    list(
      curves = list(
        DrugA = list(
          summary = tibble::tibble(
            cluster = as.factor(c(1, 2)),
            mean_x = c(0.05, 1.05),
            mean_y = c(0.05, 1.05),
            mean_concentration = c(0.01, 1),
            point_count = c(2, 2),
            proj_x = c(0.05, 1.05),
            proj_y = c(0.05, 1.05)
          ),
          ref_points = tibble::tibble(
            known_arc_length_ref = c(0, 1),
            concentration_ref = c(0.01, 1)
          ),
          curve_data = tibble::tibble(
            x = c(0, 0.5, 1),
            y = c(0, 0.5, 1),
            arc_length = c(0, 0.5, 1)
          ),
          compound = "DrugA"
        )
      ),
      point_assignments = list(
        DrugA = tibble::tibble(
          WellId = paste0("W", 1:4),
          Compound = "DrugA",
          Concentration = c(0.01, 0.01, 1, 1),
          x = c(0, 0.1, 1, 1.1),
          y = c(0, 0.1, 1, 1.1),
          lambda = c(0, 0.1, 0.9, 1)
        )
      ),
      params = list()
    ),
    class = "cellpainting_curves"
  )

  projections <- list(
    projections = tibble::tibble(
      unknown_cluster = as.factor(c(1, 2)),
      compound = "DrugA",
      distance = c(1, 2),
      arc_length = c(0.25, 0.75),
      estimated_concentration = c(0.25, 0.75),
      projected_x = c(0.25, 0.75),
      projected_y = c(0.25, 0.75),
      closest_x = c(0, 1),
      closest_y = c(0, 1),
      unknown_x = c(2.05, 3.05),
      unknown_y = c(2.05, 3.05)
    ),
    unknown_assignment = clusters$unknown_assignment,
    unknown_centroids = clusters$unknown_centroids
  )

  dimred_plot <- plotDimred(dimred)
  expect_s3_class(dimred_plot, "ggplot")
  expect_identical(dimred_plot$labels$x, "DIM1")
  expect_identical(dimred_plot$labels$y, "DIM2")
  expect_s3_class(plotUnknownClusters(clusters), "ggplot")
  expect_s3_class(plotKnownClusters(clusters, "DrugA"), "ggplot")
  expect_s3_class(plotKnownCurves(curves, "DrugA"), "ggplot")
  known_curve_plot <- plotKnownCurvesPlotly(curves, "DrugA")
  expect_s3_class(known_curve_plot, "plotly")
  expect_s3_class(plotUnknownProjections(curves, projections, "DrugA"), "ggplot")

  projection_plot <- plotUnknownProjectionsPlotly(curves, projections, "DrugA")
  expect_s3_class(projection_plot, "plotly")
  projection_traces <- plotly::plotly_build(projection_plot)$x$data
  expect_true(any(vapply(projection_traces, function(trace) identical(trace$name, "Unknown wells"), logical(1))))
  expect_true(any(vapply(projection_traces, function(trace) identical(trace$name, "Unknown centroids"), logical(1))))
  expect_false(any(vapply(
    projection_traces,
    function(trace) trace$name %in% c("Unknown wells", "DrugA wells") && isTRUE(trace$showlegend),
    logical(1)
  )))
  unknown_centroid_index <- which(vapply(
    projection_traces,
    function(trace) identical(trace$name, "Unknown centroids"),
    logical(1)
  ))[1]
  unknown_centroid_trace <- projection_traces[[unknown_centroid_index]]
  expect_false(isTRUE(unknown_centroid_trace$showlegend))
  expect_true(any(vapply(projection_traces, function(trace) isTRUE(trace$marker$showscale), logical(1))))

  known_curve_traces <- plotly::plotly_build(known_curve_plot)$x$data
  expect_true(any(vapply(known_curve_traces, function(trace) identical(trace$name, "DrugA wells"), logical(1))))
  expect_true(any(vapply(known_curve_traces, function(trace) identical(trace$name, "Known centroids"), logical(1))))
  known_centroid_index <- which(vapply(
    known_curve_traces,
    function(trace) identical(trace$name, "Known centroids"),
    logical(1)
  ))[1]
  known_centroid_trace <- known_curve_traces[[known_centroid_index]]
  expect_false(isTRUE(known_centroid_trace$showlegend))
  expect_true(any(vapply(known_curve_traces, function(trace) isTRUE(trace$marker$showscale), logical(1))))
})

test_that("compound-specific plots report available compounds", {
  clusters <- structure(
    list(
      known_assignments = list(DrugA = tibble::tibble()),
      known_clusters = list(DrugA = tibble::tibble()),
      params = list()
    ),
    class = "cellpainting_clusters"
  )

  curves <- structure(
    list(
      curves = list(DrugA = list()),
      point_assignments = list(DrugA = tibble::tibble()),
      params = list()
    ),
    class = "cellpainting_curves"
  )

  expect_error(plotKnownClusters(clusters, "DrugB"), "Available compounds: DrugA")
  expect_error(plotKnownCurves(curves, "DrugB"), "Available compounds: DrugA")
})
