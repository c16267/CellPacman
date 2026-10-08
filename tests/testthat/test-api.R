test_that("read_cellpainting_table reads data frames and validates paths", {
  tbl <- feature_fixture()

  expect_s3_class(read_cellpainting_table(tbl), "tbl_df")
  expect_error(read_cellpainting_table("missing.tsv"), "existing file")
})

test_that("suppressWarningMessage suppresses only the requested warning", {
  expect_silent(suppressWarningMessage(warning("ignore"), "ignore"))
  expect_warning(suppressWarningMessage(warning("keep"), "ignore"), "keep")
})

test_that("loadCellPainting validates feature and metadata inputs", {
  expect_s3_class(loadCellPainting(feature_fixture(), metadata_fixture()), "cellpainting_data")

  bad_feature <- feature_fixture()
  bad_feature$WellId[2] <- bad_feature$WellId[1]
  expect_error(loadCellPainting(bad_feature, metadata_fixture()), "duplicate `WellId`")

  bad_metadata <- dplyr::select(metadata_fixture(), -Plate)
  expect_error(loadCellPainting(feature_fixture(), bad_metadata), "missing required columns: Plate")

  missing_compound <- dplyr::select(metadata_fixture(), -Compound)
  expect_error(loadCellPainting(feature_fixture(), missing_compound), "missing required columns: Compound")

  missing_metadata_id <- dplyr::filter(metadata_fixture(), WellId != "W1")
  expect_error(loadCellPainting(feature_fixture(), missing_metadata_id), "missing one or more `WellId`")
})

test_that("normalize handles zero-variance features by group", {
  data <- loadCellPainting(feature_fixture(), metadata_fixture())
  norm <- normalize(data, by = "Plate")

  expect_s3_class(norm, "cellpainting_normalized")
  expect_true(all(is.finite(norm$feature_data$FeatureA)))
  expect_equal(norm$feature_data$FeatureA[1:3], c(0, 0, 0))
})

test_that("NA Compound values alone identify unknown wells", {
  expect_identical(
    is_unknown_value(c(NA_character_, "", "NA", "DrugA")),
    c(TRUE, FALSE, FALSE, FALSE)
  )
})

test_that("compoundCluster validates parameters and empty cluster cases", {
  dimred <- structure(
    list(
      reduced_data = tibble::tibble(
        WellId = paste0("W", 1:4),
        Compound = c("DrugA", "DrugA", NA, NA),
        Concentration = c(0.1, 1, NA, NA),
        x = c(0, 0.1, 10, 10.1),
        y = c(0, 0.1, 10, 10.1)
      )
    ),
    class = "cellpainting_dimred"
  )

  expect_error(compoundCluster(dimred, eps = 0), "`eps`")
  expect_error(compoundCluster(dimred, minPts = 1.5), "`minPts`")
  expect_error(compoundCluster(structure(
    list(reduced_data = dplyr::select(dimred$reduced_data, -Compound)),
    class = "cellpainting_dimred"
  )), "Compound")
  expect_error(compoundCluster(structure(
    list(reduced_data = dplyr::select(dimred$reduced_data, -y)),
    class = "cellpainting_dimred"
  )), "Missing required columns: y")
  expect_error(
    compoundCluster(structure(list(reduced_data = dplyr::select(dimred$reduced_data, -Concentration)), class = "cellpainting_dimred")),
    "Concentration"
  )

  clusters <- compoundCluster(dimred, eps = 0.001, minPts = 3)
  expect_s3_class(clusters, "cellpainting_clusters")
  expect_length(clusters$known_clusters, 0)
  expect_equal(nrow(clusters$unknown_assignment), 0)
  expect_equal(nrow(clusters$unknown_centroids), 0)
})

test_that("selectFeatureCluster returns empty results for fewer than two unknown clusters", {
  data <- loadCellPainting(feature_fixture(), metadata_fixture())
  norm <- normalize(data, by = "Plate")
  clusters <- structure(
    list(
      unknown_assignment = tibble::tibble(
        WellId = paste0("W", 1:3),
        cluster = factor(c(1, 1, 1))
      )
    ),
    class = "cellpainting_clusters"
  )

  out <- selectFeatureCluster(clusters, norm, medianDiffCutoff = 3.9)
  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 0)
  expect_named(out, c("comparison", "feature", "p_value", "median_difference", "adjusted_p_value", "significance"))
})

test_that("selectFeatureCluster median differences are robust to outliers", {
  norm <- structure(
    list(
      metadata = tibble::tibble(
        WellId = paste0("W", 1:6),
        Plate = "P1",
        Compound = NA_character_,
        Concentration = NA_real_
      ),
      feature_data = tibble::tibble(
        WellId = paste0("W", 1:6),
        FeatureA = c(-3, -2, -1, 1, 2, 100)
      )
    ),
    class = c("cellpainting_normalized", "cellpainting_data")
  )
  clusters <- structure(
    list(
      unknown_assignment = tibble::tibble(
        WellId = paste0("W", 1:6),
        cluster = factor(c(1, 1, 1, 2, 2, 2))
      )
    ),
    class = "cellpainting_clusters"
  )

  out <- selectFeatureCluster(clusters, norm, medianDiffCutoff = 3.9)
  expect_equal(nrow(out), 1)
  expect_true(is.finite(out$median_difference))
  expect_equal(out$median_difference, 4)
})

test_that("compareCompoundFeatures compares known and unknown groups", {
  data <- loadCellPainting(feature_fixture(), metadata_fixture())
  norm <- normalize(data, by = "Plate")

  out <- compareCompoundFeatures(norm, "DrugA", "__unknown__")

  expect_named(out, c("feature", "median_difference", "p_value", "adjusted_p_value", "significance"))
  expect_equal(nrow(out), 2)
  expect_true(all(is.finite(out$median_difference)))
  expect_true(all(out$adjusted_p_value >= out$p_value))
  expect_error(compareCompoundFeatures(norm, "DrugA", "DrugA"), "different compounds")
})

test_that("summarizeCompoundFeatures summarizes all compound groups", {
  data <- loadCellPainting(feature_fixture(), metadata_fixture())
  norm <- normalize(data, by = "Plate")

  out <- summarizeCompoundFeatures(norm, "FeatureA")

  expect_named(out, c("value", "feature", "compound"))
  expect_setequal(out$compound, c("DrugA", "DrugB", "Unknown"))
  expect_equal(nrow(out), 3)
})

test_that("cluster feature comparisons use selected cluster assignments", {
  norm <- structure(
    list(
      metadata = tibble::tibble(WellId = paste0("W", 1:6)),
      feature_data = tibble::tibble(
        WellId = paste0("W", 1:6),
        FeatureA = c(-3, -2, -1, 1, 2, 3)
      )
    ),
    class = c("cellpainting_normalized", "cellpainting_data")
  )
  assignments <- tibble::tibble(
    WellId = paste0("W", 1:6),
    cluster = factor(c(1, 1, 1, 2, 2, 2))
  )

  out <- compareClusterFeatures(norm, assignments, "1", "2")
  heatmap <- summarizeClusterFeatures(norm, assignments, "FeatureA")

  expect_named(out, c("feature", "median_difference", "p_value", "adjusted_p_value", "significance"))
  expect_equal(out$median_difference, 4)
  expect_setequal(heatmap$cluster, c("1", "2"))
  expect_error(compareClusterFeatures(norm, assignments, "1", "1"), "different clusters")
})

test_that("selectFeatureCurve returns statistics and trajectory heatmap data", {
  norm <- structure(
    list(
      metadata = tibble::tibble(
        WellId = paste0("W", 1:6),
        Concentration = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
      ),
      feature_data = tibble::tibble(
        WellId = paste0("W", 1:6),
        FeatureA = c(1, 2, 3, 4, 5, 6),
        FeatureB = c(6, 5, 4, 3, 2, 1)
      )
    ),
    class = c("cellpainting_normalized", "cellpainting_data")
  )
  curves <- structure(
    list(
      point_assignments = list(
        DrugA = tibble::tibble(
          WellId = paste0("W", 1:6),
          lambda = 1:6,
          trajectory_group = c("Low", "Low", "Low", "High", "High", "High")
        )
      )
    ),
    class = "cellpainting_curves"
  )

  out <- selectFeatureCurve(curves, norm, n_bins = 3, medianDiffCutoff = 2.9)

  expect_named(out, "DrugA")
  expect_named(out$DrugA, c("stats", "top_10", "heatmap_data"))
  expect_setequal(out$DrugA$stats$feature, c("FeatureA", "FeatureB"))
  expect_equal(nrow(out$DrugA$heatmap_data), 3)
})

test_that("curveProject reports curves without enough finite reference points", {
  pred_curve <- structure(
    list(
      curves = list(
        DrugA = list(
          ref_points = tibble::tibble(known_arc_length_ref = 1, concentration_ref = 0.1),
          curve_data = tibble::tibble(x = c(0, 1), y = c(0, 1), arc_length = c(0, 1))
        )
      ),
      cluster_results = list(
        unknown_assignment = tibble::tibble(WellId = "W1", cluster = factor(1), x = 0, y = 0),
        unknown_centroids = tibble::tibble(cluster = factor(1), mean_x = 0.5, mean_y = 0.5, point_count = 2)
      )
    ),
    class = "cellpainting_curves"
  )

  expect_error(curveProject(pred_curve), "enough finite reference points")
})

test_that("curveProject projects all unknown centroids onto each curve", {
  curve_data <- tibble::tibble(
    x = seq(0, 10),
    y = 0,
    arc_length = seq(0, 10)
  )
  pred_curve <- structure(
    list(
      curves = list(
        DrugA = list(
          ref_points = tibble::tibble(
            known_arc_length_ref = c(0, 10),
            concentration_ref = c(0.1, 1)
          ),
          curve_data = curve_data
        )
      ),
      cluster_results = list(
        unknown_assignment = tibble::tibble(
          WellId = c("W1", "W2"),
          cluster = factor(c(1, 2)),
          x = c(2, 8),
          y = 1
        ),
        unknown_centroids = tibble::tibble(
          cluster = factor(c(1, 2)),
          mean_x = c(2, 8),
          mean_y = c(1, 1),
          point_count = c(1L, 1L)
        )
      )
    ),
    class = "cellpainting_curves"
  )

  out <- curveProject(pred_curve)

  expect_equal(nrow(out$projections), 2)
  expect_equal(as.character(out$projections$unknown_cluster), c("1", "2"))
  expect_equal(out$projections$arc_length, c(2, 8))
  expect_equal(out$projections$estimated_concentration, c(0.28, 0.82))
  expect_equal(out$projections$distance, c(1, 1), tolerance = 1e-6)
})

test_that("readCellPainting dispatches on the file extension", {
  csv_path <- tempfile(fileext = ".csv")
  tsv_path <- tempfile(fileext = ".tsv")
  txt_path <- tempfile(fileext = ".bin")
  on.exit(unlink(c(csv_path, tsv_path, txt_path)), add = TRUE)
  readr::write_csv(feature_fixture(), csv_path)
  readr::write_tsv(feature_fixture(), tsv_path)
  readr::write_tsv(feature_fixture(), txt_path)

  expect_equal(readCellPainting(csv_path), feature_fixture())
  expect_equal(readCellPainting(tsv_path), feature_fixture())
  expect_equal(readCellPainting(txt_path, ext = "txt"), feature_fixture())
  expect_error(readCellPainting(txt_path), "Unsupported file type")
  expect_error(readCellPainting("missing.csv"), "existing file")
})

test_that("example data paths resolve to the bundled tables", {
  expect_true(file.exists(exampleDataPath("features")))
  expect_true(file.exists(exampleDataPath("metadata")))
  expect_error(exampleDataPath("other"))

  metadata <- readCellPainting(exampleDataPath("metadata"))
  expect_true(all(c("WellId", "Plate", "Compound", "Concentration") %in% names(metadata)))
  expect_equal(nrow(metadata), 1920L)
  expect_setequal(stats::na.omit(unique(metadata$Compound)), c("Axit", "Cabo", "Dacti", "DMSO"))
})

test_that("compareCompoundFeatures accepts NA for the pooled unknown wells", {
  data <- loadCellPainting(feature_fixture(), metadata_fixture())
  norm <- normalize(data, by = "Plate")

  via_na <- compareCompoundFeatures(norm, "DrugA", NA)
  via_sentinel <- compareCompoundFeatures(norm, "DrugA", "__unknown__")
  expect_equal(via_na, via_sentinel)
  expect_error(compareCompoundFeatures(norm, NA, NA), "different compounds")
  expect_error(compareCompoundFeatures(norm, "DrugA", "Missing"), "at least two observations")
})
