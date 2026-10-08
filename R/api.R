#' Load Cell Painting feature and metadata tables
#'
#' @param file.feature Path to a feature table or a data frame.
#' @param file.metadata Path to a metadata table or a data frame.
#'
#' @details `Compound` is a required metadata column. `NA` values mark unknown
#'   compounds; non-missing values identify known compounds.
#'
#' @return A `cellpainting_data` object containing feature and metadata tables.
#' @export
loadCellPainting <- function(file.feature, file.metadata) {
  feature_tbl <- validate_feature_table(read_cellpainting_table(file.feature))
  metadata_tbl <- validate_metadata_table(read_cellpainting_table(file.metadata))
  
  missing_metadata_ids <- setdiff(feature_tbl$WellId, metadata_tbl$WellId)
  if (length(missing_metadata_ids) > 0) {
    stop("`file.metadata` is missing one or more `WellId` values present in `file.feature`.")
  }
  
  structure(
    list(
      feature_data = feature_tbl,
      metadata = metadata_tbl
    ),
    class = "cellpainting_data"
  )
}

#' @rdname loadCellPainting
#' @export
load <- function(file.feature, file.metadata) {
  loadCellPainting(file.feature, file.metadata)
}

#' Normalize Cell Painting features
#'
#' @param data A `cellpainting_data` object.
#' @param method Normalization method. Currently only `"z.score"` is supported.
#' @param by Optional grouping columns used to compute normalization per group.
#'
#' @return A `cellpainting_normalized` object.
#' @export
normalize <- function(data, method = "z.score", by = "Plate") {
  if (!identical(method, "z.score")) {
    stop("Only `method = \"z.score\"` is currently supported.")
  }
  if (!inherits(data, "cellpainting_data")) {
    stop("`data` must be a `cellpainting_data` object created by `loadCellPainting()`.")
  }
  
  feature_cols <- get_feature_cols(data$feature_data)
  norm_tbl <- dplyr::left_join(data$metadata, data$feature_data, by = "WellId")
  
  if (length(by) == 0) {
    stop("`by` must specify at least 1 column.")
  }
  
  invalid_by <- setdiff(by, names(norm_tbl))
  if (length(invalid_by) > 0) {
    stop("Invalid `by` columns: ", paste(invalid_by, collapse = ", "))
  }
  
  norm_tbl <- norm_tbl %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(feature_cols), zscore_vector)) %>%
    dplyr::ungroup()
  
  structure(
    list(
      feature_data = tibble::as_tibble(norm_tbl[, c("WellId", feature_cols), drop = FALSE]),
      metadata = data$metadata,
      normalization = list(method = method, by = by)
    ),
    class = c("cellpainting_normalized", "cellpainting_data")
  )
}

#' Run dimensionality reduction
#'
#' @param data.norm A `cellpainting_normalized` object.
#' @param method Dimensionality reduction method. Currently only `"umap"` is supported.
#' @param random_seed Seed passed to the dimensionality reduction backend.
#' @param ... Additional arguments passed to the dimensionality reduction backend.
#'
#' @details Dimensionality reduction always computes two coordinates, stored as
#'   `x` and `y` in result$reduced_data.
#'
#' @return A `cellpainting_dimred` object.
#' @export
dimReduce <- function(data.norm, method = "umap", random_seed = 1, ...) {
  if (!identical(method, "umap")) {
    stop("Only `method = \"umap\"` is currently supported.")
  }
  if (!inherits(data.norm, "cellpainting_normalized")) {
    stop("`data.norm` must be a `cellpainting_normalized` object created by `normalize()`.")
  }
  feature_cols <- get_feature_cols(data.norm$feature_data)
  feature_matrix <- as.matrix(data.norm$feature_data[, feature_cols, drop = FALSE])
  if (nrow(feature_matrix) < 2 || ncol(feature_matrix) < 1) {
    stop("`data.norm` must contain at least 2 observations and 1 feature.")
  }
  if (any(!is.finite(feature_matrix))) {
    stop("Normalized feature data must contain only finite values before dimensionality reduction.")
  }
  
  set.seed(random_seed)
  dimred_res <- do.call(
    umap::umap,
    c(
      list(
        d = feature_matrix,
        n_components = 2,
        min_dist = 0.25
      ),
      list(...)
    )
  )
  
  embedding <- tibble::tibble(
    x = dimred_res$layout[, 1],
    y = dimred_res$layout[, 2]
  )
  
  embedding$WellId <- data.norm$feature_data$WellId
  dimred_tbl <- dplyr::left_join(data.norm$metadata, embedding, by = "WellId")
  
  structure(
    list(
      metadata = data.norm$metadata,
      feature_data = data.norm$feature_data,
      reduced_data = tibble::as_tibble(dimred_tbl),
      reduction = list(method = method, n_dim = 2L, seed = random_seed)
    ),
    class = "cellpainting_dimred"
  )
}

#' Cluster compounds in reduced space
#'
#' @param dimR A `cellpainting_dimred` object.
#' @param method Clustering method. Currently only `"dbscan"` is supported.
#' @param eps DBSCAN epsilon.
#' @param minPts DBSCAN `minPts`.
#'
#' @details `Compound` is required. Rows where `Compound` is `NA` are clustered
#'   together as unknown wells; non-missing values are clustered per compound.
#'
#' @return A `cellpainting_clusters` object.
#' @export
compoundCluster <- function(
    dimR,
    method = "dbscan",
    eps = 0.5,
    minPts = 5
) {
  if (!identical(method, "dbscan")) {
    stop("Only `method = \"dbscan\"` is currently supported.")
  }
  if (!inherits(dimR, "cellpainting_dimred")) {
    stop("`dimR` must be a `cellpainting_dimred` object created by `dimReduce()`.")
  }
  validate_positive_number(eps, "eps")
  minPts <- validate_positive_integer(minPts, "minPts")
  
  dim_tbl <- tibble::as_tibble(dimR$reduced_data)
  validate_dimred_cols(dim_tbl, c("Compound", "Concentration"))
  validate_xy_cols(dim_tbl)
  
  known_mask <- !is_unknown_value(dim_tbl$Compound)
  known_tbl <- dim_tbl[known_mask, , drop = FALSE]
  unknown_tbl <- dim_tbl[!known_mask, , drop = FALSE]
  
  known_results <- vector("list", 0)
  known_assignments <- vector("list", 0)
  
  for (cmp in unique(stats::na.omit(known_tbl$Compound))) {
    sub_data <- known_tbl[known_tbl$Compound == cmp, , drop = FALSE]
    if (nrow(sub_data) < minPts) {
      next
    }
    
    db_res <- dbscan::dbscan(sub_data[, c("x", "y"), drop = FALSE], eps = eps, minPts = minPts)
    sub_data$cluster <- as.factor(db_res$cluster)
    clean_data <- sub_data[sub_data$cluster != "0", , drop = FALSE]
    
    if (nrow(clean_data) == 0) {
      next
    }
    
    cluster_summary <- clean_data %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(
        mean_x = mean(x, na.rm = TRUE),
        mean_y = mean(y, na.rm = TRUE),
        mean_concentration = mean(Concentration, na.rm = TRUE),
        point_count = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::arrange(mean_concentration)
    
    known_results[[as.character(cmp)]] <- cluster_summary
    known_assignments[[as.character(cmp)]] <- tibble::as_tibble(sub_data)
  }
  
  unknown_assignment <- tibble::tibble()
  unknown_centroids <- tibble::tibble(
    cluster = factor(),
    mean_x = numeric(),
    mean_y = numeric(),
    point_count = integer()
  )
  if (nrow(unknown_tbl) >= minPts) {
    db_res <- dbscan::dbscan(unknown_tbl[, c("x", "y"), drop = FALSE], eps = eps, minPts = minPts)
    unknown_tbl$cluster <- as.factor(db_res$cluster)
    unknown_assignment <- tibble::as_tibble(unknown_tbl)
    unknown_centroids <- unknown_tbl %>%
      dplyr::filter(cluster != "0") %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(
        mean_x = mean(x, na.rm = TRUE),
        mean_y = mean(y, na.rm = TRUE),
        point_count = dplyr::n(),
        .groups = "drop"
      )
  }
  
  structure(
    list(
      reduced_data = dim_tbl,
      known_clusters = known_results,
      known_assignments = known_assignments,
      unknown_assignment = unknown_assignment,
      unknown_centroids = unknown_centroids,
      params = list(
        method = method,
        eps = eps,
        minPts = minPts
      )
    ),
    class = "cellpainting_clusters"
  )
}

#' Estimate per-compound concentration trajectory curves
#'
#' @param dimR A `cellpainting_dimred` object.
#' @param method Curve-fitting method. Currently only `"princurve"` is supported.
#' @param eps DBSCAN epsilon.
#' @param minPts DBSCAN `minPts`.
#' @param smoother Smoother passed to `princurve::principal_curve()`.
#'
#' @return A `cellpainting_curves` object.
#' @export
curveEstimate <- function(
    dimR,
    method = "princurve",
    eps = 0.5,
    minPts = 5,
    smoother = "lowess"
) {
  if (!identical(method, "princurve")) {
    stop("Only `method = \"princurve\"` is currently supported.")
  }
  validate_positive_number(eps, "eps")
  minPts <- validate_positive_integer(minPts, "minPts")
  
  clusters <- compoundCluster(
    dimR = dimR,
    method = "dbscan",
    eps = eps,
    minPts = minPts
  )
  
  concentration_curves <- vector("list", 0)
  point_assignments <- vector("list", 0)
  
  for (cmp in names(clusters$known_clusters)) {
    cluster_summary <- clusters$known_clusters[[cmp]]
    sub_data <- clusters$known_assignments[[cmp]] %>%
      dplyr::filter(cluster != "0")
    
    if (nrow(cluster_summary) < 2 || nrow(sub_data) < minPts) {
      next
    }
    validate_xy_cols(sub_data)
    
    fit <- princurve::principal_curve(
      as.matrix(sub_data[, c("x", "y"), drop = FALSE]),
      start = as.matrix(cluster_summary[, c("mean_x", "mean_y"), drop = FALSE]),
      smoother = smoother
    )
    
    fitted_curve_points <- fit$s[fit$ord, , drop = FALSE]
    curve_arc_length <- fit$lambda[fit$ord]
    finite_point_mask <- stats::complete.cases(fitted_curve_points) &
      is.finite(fitted_curve_points[, 1]) &
      is.finite(fitted_curve_points[, 2]) &
      is.finite(curve_arc_length)
    if (sum(finite_point_mask) < 2) {
      next
    }
    fitted_curve_points <- fitted_curve_points[finite_point_mask, , drop = FALSE]
    curve_arc_length <- curve_arc_length[finite_point_mask]
    
    centroid_matrix <- as.matrix(cluster_summary[, c("mean_x", "mean_y"), drop = FALSE])
    if (nrow(centroid_matrix) < 2 || any(!is.finite(centroid_matrix))) {
      next
    }
    centroid_proj <- princurve::project_to_curve(centroid_matrix, fitted_curve_points, stretch = 0)
    
    synced_lambdas <- vapply(
      seq_len(nrow(centroid_matrix)),
      function(i) {
        idx <- which.min(rowSums((fitted_curve_points - matrix(centroid_matrix[i, ], nrow(fitted_curve_points), 2, byrow = TRUE))^2))
        curve_arc_length[idx]
      },
      numeric(1)
    )
    
    cluster_summary$proj_x <- centroid_proj$s[, 1]
    cluster_summary$proj_y <- centroid_proj$s[, 2]
    
    sub_data$lambda <- fit$lambda
    median_trajectory_position <- stats::median(sub_data$lambda, na.rm = TRUE)
    sub_data$trajectory_group <- ifelse(sub_data$lambda > median_trajectory_position, "High", "Low")
    
    concentration_curves[[cmp]] <- list(
      summary = tibble::as_tibble(cluster_summary),
      ref_points = tibble::tibble(
        known_arc_length_ref = synced_lambdas,
        concentration_ref = cluster_summary$mean_concentration
      ),
      curve_data = tibble::tibble(
        x = fitted_curve_points[, 1],
        y = fitted_curve_points[, 2],
        arc_length = curve_arc_length
      ),
      fit = fit,
      compound = cmp
    )
    
    point_assignments[[cmp]] <- tibble::as_tibble(sub_data)
  }
  
  structure(
    list(
      curves = concentration_curves,
      point_assignments = point_assignments,
      cluster_results = clusters,
      params = list(
        method = method,
        eps = eps,
        minPts = minPts,
        smoother = smoother
      )
    ),
    class = "cellpainting_curves"
  )
}

#' Project unknown compound clusters onto known concentration trajectory curves
#'
#' @param predCurve A `cellpainting_curves` object.
#' @param newData Optional `cellpainting_dimred` object.
#' @param eps DBSCAN epsilon for unknown clustering.
#' @param minPts DBSCAN `minPts` for unknown clustering.
#'
#' @return A list containing all projections and per-cluster best matches.
#' @export
curveProject <- function(
    predCurve,
    newData = NULL,
    eps = 0.5,
    minPts = 5
) {
  if (!inherits(predCurve, "cellpainting_curves")) {
    stop("`predCurve` must be a `cellpainting_curves` object.")
  }
  validate_positive_number(eps, "eps")
  minPts <- validate_positive_integer(minPts, "minPts")
  trajectory_curve_result <- predCurve
  
  clusters <- if (is.null(newData)) {
    trajectory_curve_result$cluster_results
  } else {
    compoundCluster(
      dimR = newData,
      method = "dbscan",
      eps = eps,
      minPts = minPts
    )
  }
  
  na_centroids <- clusters$unknown_centroids
  if (is.null(na_centroids) || nrow(na_centroids) == 0) {
    stop("No unknown compound clusters were available for projection.")
  }
  validate_xy_cols(tibble::tibble(x = na_centroids$mean_x, y = na_centroids$mean_y))
  
  centroid_matrix <- as.matrix(na_centroids[, c("mean_x", "mean_y"), drop = FALSE])
  projection_results <- lapply(names(trajectory_curve_result$curves), function(cmp) {
    compound_curve <- trajectory_curve_result$curves[[cmp]]
    trajectory_curve_data <- tibble::as_tibble(compound_curve$curve_data)
    trajectory_curve_data <- trajectory_curve_data[finite_trajectory_curve_points(trajectory_curve_data), , drop = FALSE]
    if (nrow(trajectory_curve_data) < 2) {
      return(NULL)
    }
    trajectory_curve_points <- as.matrix(trajectory_curve_data[, c("x", "y"), drop = FALSE])

    concentration_references <- dplyr::distinct(
      compound_curve$ref_points,
      known_arc_length_ref,
      .keep_all = TRUE
    )
    concentration_references <- concentration_references[
      finite_trajectory_reference_points(concentration_references),
      ,
      drop = FALSE
    ]
    if (nrow(concentration_references) < 2) {
      return(NULL)
    }

    projected_centroids <- princurve::project_to_curve(centroid_matrix, trajectory_curve_points, stretch = 0)
    squared_distances <- outer(rowSums(centroid_matrix^2), rowSums(trajectory_curve_points^2), "+") -
      2 * tcrossprod(centroid_matrix, trajectory_curve_points)
    closest_indices <- max.col(-squared_distances, ties.method = "first")
    arc_lengths <- trajectory_curve_data$arc_length[closest_indices]
    estimated_concentrations <- stats::approx(
      x = concentration_references$known_arc_length_ref,
      y = concentration_references$concentration_ref,
      xout = arc_lengths,
      rule = 2
    )$y

    tibble::tibble(
      unknown_cluster = na_centroids$cluster,
      compound = cmp,
      distance = sqrt(rowSums((centroid_matrix - projected_centroids$s)^2)),
      arc_length = arc_lengths,
      estimated_concentration = estimated_concentrations,
      projected_x = projected_centroids$s[, 1],
      projected_y = projected_centroids$s[, 2],
      closest_x = trajectory_curve_points[closest_indices, 1],
      closest_y = trajectory_curve_points[closest_indices, 2],
      unknown_x = centroid_matrix[, 1],
      unknown_y = centroid_matrix[, 2]
    )
  })
  
  all_matches <- dplyr::bind_rows(projection_results)
  if (nrow(all_matches) == 0) {
    stop("No fitted concentration trajectory curves had enough finite reference points for projection.")
  }
  all_matches <- all_matches %>%
    dplyr::mutate(
      .centroid_order = match(as.character(.data$unknown_cluster), as.character(na_centroids$cluster)),
      .compound_order = match(.data$compound, names(trajectory_curve_result$curves))
    ) %>%
    dplyr::arrange(.data$.centroid_order, .data$.compound_order) %>%
    dplyr::select(-dplyr::all_of(c(".centroid_order", ".compound_order")))
  best_matches <- all_matches %>%
    dplyr::group_by(unknown_cluster) %>%
    dplyr::slice_min(distance, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  list(
    projections = all_matches,
    best_matches = best_matches,
    unknown_assignment = clusters$unknown_assignment,
    unknown_centroids = na_centroids
  )
}

#' Select cluster-discriminating features among unknown clusters
#'
#' @param predClust A `cellpainting_clusters` object.
#' @param data.norm A `cellpainting_normalized` object.
#' @param method Statistical test. Currently only `"wilcox"` is supported.
#' @param p.adjust Multiple-testing correction method.
#' @param cutoff Adjusted p-value cutoff.
#' @param medianDiffCutoff Median-difference cutoff.
#'
#' @return A tibble of pairwise cluster feature statistics.
#' @export
selectFeatureCluster <- function(
    predClust,
    data.norm,
    method = "wilcox",
    p.adjust = "bonferroni",
    cutoff = 0.05,
    medianDiffCutoff = 0.5
) {
  if (!identical(method, "wilcox")) {
    stop("Only `method = \"wilcox\"` is currently supported.")
  }
  
  if (!inherits(predClust, "cellpainting_clusters")) {
    stop("`predClust` must be a `cellpainting_clusters` object.")
  }
  cluster_obj <- predClust
  if (!inherits(data.norm, "cellpainting_normalized")) {
    stop("`data.norm` must be a `cellpainting_normalized` object created by `normalize()`.")
  }
  
  unknown_tbl <- cluster_obj$unknown_assignment
  if (is.null(unknown_tbl) || nrow(unknown_tbl) == 0) {
    stop("`predClust` does not contain clustered unknown compounds.")
  }
  
  merged <- dplyr::left_join(
    dplyr::left_join(data.norm$metadata, data.norm$feature_data, by = "WellId"),
    unknown_tbl[, c("WellId", "cluster"), drop = FALSE],
    by = "WellId"
  ) %>%
    dplyr::filter(!is.na(cluster), cluster != "0")
  
  clusters <- sort(unique(merged$cluster))
  feature_cols <- get_feature_cols(data.norm$feature_data)
  if (length(clusters) < 2 || length(feature_cols) == 0) {
    return(empty_feature_cluster_results())
  }
  
  comparisons <- vector("list", 0)
  
  for (i in seq_len(length(clusters) - 1)) {
    for (j in (i + 1):length(clusters)) {
      c1 <- clusters[i]
      c2 <- clusters[j]
      pair_data <- merged[merged$cluster %in% c(c1, c2), , drop = FALSE]
      
      stats_list <- lapply(feature_cols, function(feat) {
        x1 <- pair_data[[feat]][pair_data$cluster == c1]
        x2 <- pair_data[[feat]][pair_data$cluster == c2]
        
        if (!is.numeric(x1) || !is.numeric(x2) || length(stats::na.omit(x1)) < 2 || length(stats::na.omit(x2)) < 2) {
          return(NULL)
        }
        
        p_val <- stats::wilcox.test(x1, x2, exact = FALSE)$p.value
        median_difference <- stats::median(x2, na.rm = TRUE) - stats::median(x1, na.rm = TRUE)
        
        tibble::tibble(
          comparison = paste0("cluster_", c1, "_vs_", c2),
          feature = feat,
          p_value = p_val,
          median_difference = median_difference
        )
      })
      
      comp_df <- dplyr::bind_rows(stats_list)
      if (nrow(comp_df) == 0) {
        next
      }
      
      comp_df$adjusted_p_value <- stats::p.adjust(comp_df$p_value, method = p.adjust)
      comp_df$significance <- dplyr::case_when(
        comp_df$adjusted_p_value < cutoff & comp_df$median_difference > medianDiffCutoff ~ "Up",
        comp_df$adjusted_p_value < cutoff & comp_df$median_difference < -medianDiffCutoff ~ "Down",
        comp_df$adjusted_p_value < cutoff ~ "Pass p-value cutoff",
        TRUE ~ "Not Significant"
      )
      
      comparisons[[paste0(c1, "_", c2)]] <- comp_df
    }
  }
  
  out <- dplyr::bind_rows(comparisons)
  if (nrow(out) == 0) {
    return(empty_feature_cluster_results())
  }
  
  out %>%
    dplyr::arrange(comparison, adjusted_p_value)
}

#' Select curve-associated features for known compounds
#'
#' @param predCurve A `cellpainting_curves` object.
#' @param data.norm A `cellpainting_normalized` object.
#' @param method Statistical test. Currently only `"wilcox"` is supported.
#' @param p.adjust Multiple-testing correction method.
#' @param cutoff Adjusted p-value cutoff.
#' @param medianDiffCutoff Median-difference cutoff.
#' @param n_bins Number of bins used for heatmap-ready summaries.
#'
#' @return A named list of per-compound feature-selection results.
#' @export
selectFeatureCurve <- function(
    predCurve,
    data.norm,
    method = "wilcox",
    p.adjust = "bonferroni",
    cutoff = 0.05,
    medianDiffCutoff = 0.5,
    n_bins = 20
) {
  if (!identical(method, "wilcox")) {
    stop("Only `method = \"wilcox\"` is currently supported.")
  }
  
  if (!inherits(predCurve, "cellpainting_curves")) {
    stop("`predCurve` must be a `cellpainting_curves` object.")
  }
  trajectory_curve_result <- predCurve
  if (!inherits(data.norm, "cellpainting_normalized")) {
    stop("`data.norm` must be a `cellpainting_normalized` object created by `normalize()`.")
  }
  n_bins <- validate_positive_integer(n_bins, "n_bins")
  feature_cols <- get_feature_cols(data.norm$feature_data)
  results <- vector("list", 0)
  
  for (cmp in names(trajectory_curve_result$point_assignments)) {
    trajectory_points <- trajectory_curve_result$point_assignments[[cmp]]
    if (nrow(trajectory_points) == 0) {
      next
    }
    
    sub_data <- dplyr::left_join(
      trajectory_points[, c("WellId", "lambda", "trajectory_group"), drop = FALSE],
      dplyr::left_join(data.norm$metadata, data.norm$feature_data, by = "WellId"),
      by = "WellId"
    )
    
    stats_list <- lapply(feature_cols, function(feat) {
      group_low <- sub_data[[feat]][sub_data$trajectory_group == "Low"]
      group_high <- sub_data[[feat]][sub_data$trajectory_group == "High"]
      
      if (
        !is.numeric(group_low) ||
        !is.numeric(group_high) ||
        length(stats::na.omit(group_low)) < 2 ||
        length(stats::na.omit(group_high)) < 2 ||
        length(unique(stats::na.omit(group_low))) < 2 ||
        length(unique(stats::na.omit(group_high))) < 2
      ) {
        return(NULL)
      }
      
      test_res <- stats::wilcox.test(group_low, group_high, exact = FALSE)
      med_diff <- stats::median(group_high, na.rm = TRUE) - stats::median(group_low, na.rm = TRUE)
      
      tibble::tibble(
        feature = feat,
        p_value = test_res$p.value,
        median_difference = med_diff
      )
    })
    
    stats_df <- dplyr::bind_rows(stats_list)
    if (nrow(stats_df) == 0) {
      next
    }
    
    stats_df$adjusted_p_value <- stats::p.adjust(stats_df$p_value, method = p.adjust)
    stats_df$significance <- dplyr::case_when(
      stats_df$adjusted_p_value < cutoff & stats_df$median_difference > medianDiffCutoff ~ "Positive",
      stats_df$adjusted_p_value < cutoff & stats_df$median_difference < -medianDiffCutoff ~ "Negative",
      stats_df$adjusted_p_value < cutoff ~ "Pass p-value cutoff",
      TRUE ~ "Not Significant"
    )
    
    top_10 <- stats_df %>%
      dplyr::arrange(adjusted_p_value) %>%
      dplyr::slice_head(n = 10)
    
    sub_data_binned <- sub_data %>%
      dplyr::arrange(lambda) %>%
      dplyr::mutate(bin = dplyr::ntile(lambda, n_bins))
    
    top_feats <- top_10$feature
    heatmap_ready <- sub_data_binned %>%
      dplyr::group_by(bin) %>%
      dplyr::summarise(
        dplyr::across(dplyr::all_of(top_feats), \(x) mean(x, na.rm = TRUE)),
        mean_concentration = mean(Concentration, na.rm = TRUE),
        .groups = "drop"
      )
    
    results[[cmp]] <- list(
      stats = stats_df,
      top_10 = top_10,
      heatmap_data = heatmap_ready
    )
  }
  
  results
}
