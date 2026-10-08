#' Paths to the bundled example Cell Painting tables
#'
#' `cellpacman` ships a small Cell Painting screen in `inst/extdata`: 1,920
#' wells from five 384-well plates, with 664 well-level morphological features
#' (`example_features.tsv`) and the matching plate annotation
#' (`example_metadata.tsv`). Each plate holds three annotated compounds
#' (`Axit` = axitinib, `Cabo` = cabozantinib, `Dacti` = dactinomycin) at
#' eight concentrations in duplicate, 16 wells of the `DMSO` vehicle control,
#' and 320 unannotated wells whose `Compound` is `NA` and which therefore
#' enter the workflow as unknown compounds.
#'
#' @param table Which table to locate: `"features"` or `"metadata"`.
#'
#' @return A single file path to a tab-separated table that can be passed to
#'   [loadCellPainting()].
#' @examples
#' exampleDataPath("features")
#' exampleDataPath("metadata")
#' @export
exampleDataPath <- function(table = c("features", "metadata")) {
  table <- match.arg(table)
  example_data_path(switch(
    table,
    features = "example_features.tsv",
    metadata = "example_metadata.tsv"
  ))
}

# Resolve a file name inside the installed `extdata` directory. `system.file()`
# also resolves correctly under `devtools::load_all()`.
example_data_path <- function(file) {
  path <- system.file("extdata", file, package = "cellpacman", mustWork = FALSE)
  if (!nzchar(path)) {
    stop("Example data file `", file, "` was not found in the installed package.")
  }
  path
}
