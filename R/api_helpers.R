suppressWarningMessage <- function(expr, message) {
  withCallingHandlers(
    expr,
    warning = function(condition) {
      if (identical(conditionMessage(condition), message)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

#' Read a Cell Painting feature or metadata table from disk
#'
#' Reads one delimited or Excel table into a tibble. The format is chosen
#' from `ext`, which defaults to the file extension: comma-separated (`csv`),
#' tab-separated (`tsv`, `txt`), or Excel (`xlsx`). Column names are kept
#' exactly as stored, so feature names containing spaces (as exported by
#' Harmony and similar software) are preserved.
#'
#' [loadCellPainting()] calls this function internally when it is given file
#' paths; `readCellPainting()` is exported so that applications can parse an
#' uploaded file before validation, for example to preview it.
#'
#' @param file Path to the table.
#' @param ext File type used to pick the parser. Defaults to the extension of
#'   `file`; pass it explicitly when the file name does not carry the right
#'   extension (e.g. a temporary upload).
#'
#' @return A tibble with the table contents.
#' @examples
#' metadata <- readCellPainting(exampleDataPath("metadata"))
#' metadata[, c("WellId", "Plate", "Compound", "Concentration")]
#' @export
readCellPainting <- function(file, ext = NULL) {
  if (!is.character(file) || length(file) != 1 || is.na(file) || !file.exists(file)) {
    stop("`file` must be the path to an existing file.")
  }
  if (is.null(ext)) {
    ext <- tools::file_ext(file)
  }
  ext <- tolower(ext)
  out <- switch(
    ext,
    csv = readr::read_csv(file, show_col_types = FALSE),
    tsv = readr::read_tsv(file, show_col_types = FALSE),
    txt = readr::read_tsv(file, show_col_types = FALSE),
    xlsx = readxl::read_excel(file),
    stop("Unsupported file type `", ext, "`. Use csv, tsv, txt, or xlsx.")
  )
  tibble::as_tibble(out)
}

read_cellpainting_table <- function(x) {
  if (inherits(x, "data.frame")) {
    return(tibble::as_tibble(x))
  }

  if (!is.character(x) || length(x) != 1 || !file.exists(x)) {
    stop("Input must be a data frame or a path to an existing file.")
  }

  readCellPainting(x)
}

validate_feature_table <- function(tbl) {
  if (ncol(tbl) < 2) {
    stop("`file.feature` must contain `WellId` plus at least one numeric feature column.")
  }
  
  if (!identical(names(tbl)[1], "WellId")) {
    stop("`file.feature` must contain `WellId` as its first column.")
  }
  
  feature_cols <- names(tbl)[-1]
  non_numeric <- feature_cols[!vapply(tbl[feature_cols], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop("All columns in `file.feature` after `WellId` must be numeric feature columns.")
  }
  
  if (anyDuplicated(tbl$WellId) > 0) {
    stop("`file.feature` contains duplicate `WellId` values.")
  }
  
  tibble::as_tibble(tbl)
}

validate_metadata_table <- function(tbl, required_cols = c("WellId", "Plate", "Compound", "Concentration")) {
  if (!"WellId" %in% names(tbl)) {
    stop("`file.metadata` must contain a `WellId` column.")
  }
  
  missing_cols <- setdiff(required_cols, names(tbl))
  if (length(missing_cols) > 0) {
    stop("`file.metadata` is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (anyDuplicated(tbl$WellId) > 0) {
    stop("`file.metadata` contains duplicate `WellId` values.")
  }
  
  if ("Concentration" %in% names(tbl) && !is.numeric(tbl$Concentration)) {
    stop("`file.metadata` column `Concentration` must be numeric.")
  }
  
  tibble::as_tibble(tbl)
}

get_feature_cols <- function(feature_tbl) {
  feature_cols <- setdiff(names(feature_tbl), "WellId")
  if (length(feature_cols) == 0) {
    stop("No feature columns were found.")
  }
  
  non_numeric <- feature_cols[!vapply(feature_tbl[feature_cols], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop("Feature columns must be numeric.")
  }
  
  feature_cols
}

zscore_vector <- function(col) {
  if (!is.numeric(col)) {
    return(col)
  }
  
  sd_val <- stats::sd(col, na.rm = TRUE)
  if (is.na(sd_val) || sd_val == 0) {
    return(rep(0, length(col)))
  }
  
  (col - mean(col, na.rm = TRUE)) / sd_val
}

is_unknown_value <- function(x) {
  is.na(x)
}

validate_dimred_cols <- function(tbl, cols) {
  missing_cols <- setdiff(cols, names(tbl))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
}

validate_positive_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || !is.finite(x) || x <= 0) {
    stop("`", name, "` must be a single finite number greater than 0.")
  }
  
  invisible(x)
}

validate_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || !is.finite(x) || x < 1 || x != as.integer(x)) {
    stop("`", name, "` must be a single positive integer.")
  }
  
  invisible(as.integer(x))
}

validate_xy_cols <- function(tbl) {
  xy_cols <- c("x", "y")
  validate_dimred_cols(tbl, xy_cols)
  
  non_numeric <- xy_cols[!vapply(tbl[xy_cols], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop("`x` and `y` must be numeric columns.")
  }
  
  non_finite <- xy_cols[vapply(tbl[xy_cols], function(x) any(!is.finite(x)), logical(1))]
  if (length(non_finite) > 0) {
    stop("`x` and `y` must contain only finite values.")
  }
  
  invisible(xy_cols)
}

empty_feature_cluster_results <- function() {
  tibble::tibble(
    comparison = character(),
    feature = character(),
    p_value = numeric(),
    median_difference = numeric(),
    adjusted_p_value = numeric(),
    significance = character()
  )
}

finite_trajectory_curve_points <- function(trajectory_curve_data) {
  required <- c("x", "y", "arc_length")
  validate_dimred_cols(trajectory_curve_data, required)
  stats::complete.cases(trajectory_curve_data[, required, drop = FALSE]) &
    is.finite(trajectory_curve_data$x) &
    is.finite(trajectory_curve_data$y) &
    is.finite(trajectory_curve_data$arc_length)
}

finite_trajectory_reference_points <- function(ref_points) {
  required <- c("known_arc_length_ref", "concentration_ref")
  validate_dimred_cols(ref_points, required)
  stats::complete.cases(ref_points[, required, drop = FALSE]) &
    is.finite(ref_points$known_arc_length_ref) &
    is.finite(ref_points$concentration_ref)
}
