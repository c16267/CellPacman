#' cellpacman: PAinted Cell and coMpound ANalysis
#'
#' Concentration-response analysis of Cell Painting morphological profiles.
#' The workflow is a chain of S3 objects, each produced by one exported
#' function:
#'
#' 1. [loadCellPainting()] reads and validates a well-level feature table and
#'    its plate metadata (`cellpainting_data`).
#' 2. [normalize()] z-scores every feature within plate
#'    (`cellpainting_normalized`).
#' 3. [dimReduce()] embeds the normalized profiles in two dimensions with
#'    UMAP (`cellpainting_dimred`).
#' 4. [compoundCluster()] runs DBSCAN separately for each known compound and
#'    jointly for the unknown wells (`cellpainting_clusters`).
#' 5. [curveEstimate()] fits a principal curve through the
#'    concentration-ordered cluster centroids of each known compound
#'    (`cellpainting_curves`).
#' 6. [curveProject()] projects unknown-cluster centroids onto every curve and
#'    interpolates an effective concentration.
#' 7. [selectFeatureCurve()], [selectFeatureCluster()],
#'    [compareCompoundFeatures()] and [compareClusterFeatures()] identify the
#'    features behind a trajectory or a cluster contrast.
#'
#' Plotting helpers ([plotDimred()], [plotKnownClusters()],
#' [plotUnknownClusters()], [plotKnownCurves()], [plotUnknownProjections()]
#' and their Plotly variants) visualize each stage. The bundled example screen
#' is available through [exampleDataPath()]. An interactive Shiny front end
#' that drives this package is distributed separately as the
#' `cellpacman-shiny` application; see `vignette("shiny-app")`.
#'
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  "Concentration", "adjusted_p_value", "bin", "cluster", "comparison", "compound",
  "compound_label", "concentration_factor", "distance", "estimated_concentration",
  "feature", "known_arc_length_ref", "lambda", "mean_concentration", "mean_x", "mean_y",
  "median_difference", "p_value", "proj_x", "proj_y", "projected_x", "projected_y",
  "significance", "unknown_cluster", "unknown_x", "unknown_y", "value", "x", "y"
))
