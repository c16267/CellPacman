# Pairwise feature comparisons and per-group feature summaries. These
# functions back the "Compare compounds" and "Compare clusters" panels of the
# companion Shiny application, but they are ordinary analysis functions that
# operate on the objects returned by normalize() and compoundCluster().

#' Compare features between two compound groups
#'
#' Tests every normalized feature for a location shift between the wells of
#' two compounds with a two-sided Wilcoxon rank-sum test (normal
#' approximation), reports the difference in medians
#' (`second` minus `first`) as the effect size, and adjusts the p-values
#' across features.
#'
#' @param data.norm A `cellpainting_normalized` object created by
#'   [normalize()].
#' @param first,second Compound names as they appear in the `Compound`
#'   metadata column. `NA` (or the string `"__unknown__"`) denotes the pooled
#'   unannotated wells, i.e. all wells whose `Compound` is `NA`.
#' @param p.adjust Multiple-testing correction passed to [stats::p.adjust()].
#' @param cutoff Adjusted p-value cutoff used to call a feature significant.
#' @param medianDiffCutoff Absolute median-difference cutoff used to label a
#'   significant feature `"Positive"` or `"Negative"`.
#'
#' @return A tibble with one row per feature, ordered by adjusted p-value:
#'   `feature`, `median_difference`, `p_value`, `adjusted_p_value`, and
#'   `significance` (`"Positive"`, `"Negative"`, `"Pass p-value cutoff"`, or
#'   `"Not significant"`). Features with fewer than two finite values in
#'   either group are omitted.
#' @seealso [compareClusterFeatures()] for the analogous comparison between
#'   two DBSCAN clusters, [summarizeCompoundFeatures()] for the per-compound
#'   medians used to display the significant features.
#' @examples
#' data <- loadCellPainting(exampleDataPath("features"), exampleDataPath("metadata"))
#' data.norm <- normalize(data, by = "Plate")
#' head(compareCompoundFeatures(data.norm, "DMSO", "Cabo"))
#' @export
compareCompoundFeatures <- function(
    data.norm,
    first,
    second,
    p.adjust = "BH",
    cutoff = 0.05,
    medianDiffCutoff = 0.5
) {
  if (!inherits(data.norm, "cellpainting_normalized")) {
    stop("`data.norm` must be a `cellpainting_normalized` object created by `normalize()`.")
  }
  first <- validate_compound_group(first, "first")
  second <- validate_compound_group(second, "second")
  if (identical(first, second)) {
    stop("Select two different compounds to compare.")
  }

  comparison_data <- dplyr::left_join(
    data.norm$metadata,
    data.norm$feature_data,
    by = "WellId"
  )
  comparison_data$comparison_group <- ifelse(
    is_unknown_value(comparison_data$Compound),
    "__unknown__",
    as.character(comparison_data$Compound)
  )
  comparison_data <- comparison_data[
    comparison_data$comparison_group %in% c(first, second),
    ,
    drop = FALSE
  ]

  group_counts <- table(comparison_data$comparison_group)
  if (!all(c(first, second) %in% names(group_counts)) || any(group_counts[c(first, second)] < 2)) {
    stop("Each selected compound must contain at least two observations.")
  }

  pairwise_feature_tests(
    comparison_data,
    get_feature_cols(data.norm$feature_data),
    group_col = "comparison_group",
    first = first,
    second = second,
    p.adjust = p.adjust,
    cutoff = cutoff,
    medianDiffCutoff = medianDiffCutoff
  )
}

#' Compare features between two DBSCAN clusters
#'
#' Tests every normalized feature for a location shift between the wells of
#' two clusters of the same compound (or of the unknown wells), using the same
#' Wilcoxon rank-sum test, median-difference effect size, and multiplicity
#' correction as [compareCompoundFeatures()].
#'
#' @inheritParams compareCompoundFeatures
#' @param assignments A cluster assignment table with `WellId` and `cluster`
#'   columns, i.e. one element of `predClust$known_assignments` or
#'   `predClust$unknown_assignment` from [compoundCluster()].
#' @param first,second Cluster labels to compare, as character strings
#'   (`"1"`, `"2"`, ...). Cluster `"0"` is DBSCAN noise.
#'
#' @return A tibble with the same columns as [compareCompoundFeatures()].
#' @examples
#' \donttest{
#' data <- loadCellPainting(exampleDataPath("features"), exampleDataPath("metadata"))
#' data.norm <- normalize(data, by = "Plate")
#' dimR <- dimReduce(data.norm, random_seed = 1)
#' predClust <- compoundCluster(dimR, eps = 0.5, minPts = 5)
#' unknown <- predClust$unknown_assignment
#' head(compareClusterFeatures(data.norm, unknown, "1", "2"))
#' }
#' @export
compareClusterFeatures <- function(
    data.norm,
    assignments,
    first,
    second,
    p.adjust = "BH",
    cutoff = 0.05,
    medianDiffCutoff = 0.5
) {
  if (!inherits(data.norm, "cellpainting_normalized")) {
    stop("`data.norm` must be a `cellpainting_normalized` object created by `normalize()`.")
  }
  validate_dimred_cols(assignments, c("WellId", "cluster"))
  first <- as.character(first)
  second <- as.character(second)
  if (identical(first, second)) {
    stop("Select two different clusters to compare.")
  }

  comparison_data <- dplyr::left_join(
    data.norm$feature_data,
    assignments[, c("WellId", "cluster"), drop = FALSE],
    by = "WellId"
  )
  comparison_data$cluster <- as.character(comparison_data$cluster)
  comparison_data <- comparison_data[
    !is.na(comparison_data$cluster) & comparison_data$cluster %in% c(first, second),
    ,
    drop = FALSE
  ]

  pairwise_feature_tests(
    comparison_data,
    get_feature_cols(data.norm$feature_data),
    group_col = "cluster",
    first = first,
    second = second,
    p.adjust = p.adjust,
    cutoff = cutoff,
    medianDiffCutoff = medianDiffCutoff
  )
}

#' Summarize selected features by compound or by cluster
#'
#' Computes the median normalized value of each selected feature within each
#' group, in the long format expected by heatmap functions. Unannotated wells
#' form the `"Unknown"` compound group; DBSCAN noise (cluster `"0"`) is
#' excluded from cluster summaries.
#'
#' @inheritParams compareClusterFeatures
#' @param features Character vector of feature names to summarize, typically
#'   the significant features returned by [compareCompoundFeatures()] or
#'   [compareClusterFeatures()].
#'
#' @return A tibble with columns `value`, `feature`, and `compound` (or
#'   `cluster`), one row per feature-group combination.
#' @examples
#' data <- loadCellPainting(exampleDataPath("features"), exampleDataPath("metadata"))
#' data.norm <- normalize(data, by = "Plate")
#' top <- head(compareCompoundFeatures(data.norm, "DMSO", "Cabo")$feature, 5)
#' summarizeCompoundFeatures(data.norm, top)
#' @export
summarizeCompoundFeatures <- function(data.norm, features) {
  if (!inherits(data.norm, "cellpainting_normalized")) {
    stop("`data.norm` must be a `cellpainting_normalized` object created by `normalize()`.")
  }
  if (length(features) == 0) {
    return(tibble::tibble(value = numeric(), feature = character(), compound = character()))
  }

  heatmap_data <- dplyr::left_join(data.norm$metadata, data.norm$feature_data, by = "WellId")
  heatmap_data$compound <- ifelse(
    is_unknown_value(heatmap_data$Compound),
    "Unknown",
    as.character(heatmap_data$Compound)
  )
  summary_data <- heatmap_data %>%
    dplyr::group_by(compound) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(features), ~ stats::median(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  summary_data %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(features),
      names_to = "feature",
      values_to = "value",
      cols_vary = "slowest"
    ) %>%
    dplyr::select(dplyr::all_of(c("value", "feature", "compound")))
}

#' @rdname summarizeCompoundFeatures
#' @export
summarizeClusterFeatures <- function(data.norm, assignments, features) {
  if (!inherits(data.norm, "cellpainting_normalized")) {
    stop("`data.norm` must be a `cellpainting_normalized` object created by `normalize()`.")
  }
  validate_dimred_cols(assignments, c("WellId", "cluster"))
  if (length(features) == 0) {
    return(tibble::tibble(value = numeric(), feature = character(), cluster = character()))
  }

  summary_data <- dplyr::left_join(
    data.norm$feature_data,
    assignments[, c("WellId", "cluster"), drop = FALSE],
    by = "WellId"
  ) %>%
    dplyr::filter(!is.na(cluster), cluster != "0") %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(features), ~ stats::median(.x, na.rm = TRUE)), .groups = "drop")

  summary_data %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(features),
      names_to = "feature",
      values_to = "value",
      cols_vary = "slowest"
    ) %>%
    dplyr::mutate(cluster = as.character(cluster)) %>%
    dplyr::select(dplyr::all_of(c("value", "feature", "cluster")))
}

# Shared engine: per-feature Wilcoxon rank-sum test between two groups of a
# wide table, median-difference effect size, and significance labels.
pairwise_feature_tests <- function(data, feature_cols, group_col, first, second, p.adjust, cutoff, medianDiffCutoff) {
  results <- lapply(feature_cols, function(feature) {
    first_values <- data[[feature]][data[[group_col]] == first]
    second_values <- data[[feature]][data[[group_col]] == second]
    first_values <- first_values[is.finite(first_values)]
    second_values <- second_values[is.finite(second_values)]

    if (length(first_values) < 2 || length(second_values) < 2) {
      return(NULL)
    }

    tibble::tibble(
      feature = feature,
      median_difference = stats::median(second_values) - stats::median(first_values),
      p_value = stats::wilcox.test(first_values, second_values, exact = FALSE)$p.value
    )
  })
  results <- dplyr::bind_rows(results)
  if (nrow(results) == 0) {
    return(tibble::tibble(
      feature = character(), median_difference = numeric(), p_value = numeric(),
      adjusted_p_value = numeric(), significance = character()
    ))
  }

  results %>%
    dplyr::mutate(
      adjusted_p_value = stats::p.adjust(p_value, method = p.adjust),
      significance = dplyr::case_when(
        adjusted_p_value < cutoff & median_difference > medianDiffCutoff ~ "Positive",
        adjusted_p_value < cutoff & median_difference < -medianDiffCutoff ~ "Negative",
        adjusted_p_value < cutoff ~ "Pass p-value cutoff",
        TRUE ~ "Not significant"
      )
    ) %>%
    dplyr::arrange(adjusted_p_value)
}

validate_compound_group <- function(x, name) {
  if (length(x) != 1) {
    stop("`", name, "` must be a single compound name.")
  }
  if (is.na(x)) {
    return("__unknown__")
  }
  if (!is.character(x) || x == "") {
    stop("`", name, "` must be a single compound name.")
  }
  x
}
