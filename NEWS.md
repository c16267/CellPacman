# cellpacman 0.1.0

## Package structure

* The package is now a standalone analysis library. The Shiny interface was
  moved out of the package into the sibling `cellpacman-shiny` application,
  which calls the exported API; `shiny`, `bslib`, `reactable`,
  `shinycssloaders`, `future`, `promises`, and `digest` are no longer
  package dependencies, and `runCellPACMAN()` has been replaced by
  `shiny::runApp("cellpacman-shiny")`.
* The example screen moved from the non-standard `test_data/` directory to
  `inst/extdata/` so that it is installed with the package; the duplicate copy
  under `vignettes/` was removed. `exampleDataPath()` returns the installed
  paths.

## New exported functions

* `readCellPainting()` reads a single `csv`/`tsv`/`txt`/`xlsx` table (the
  parser used by `loadCellPainting()`).
* `compareCompoundFeatures()` and `compareClusterFeatures()` test every
  feature between two compounds (where `NA` denotes the pooled unknown wells)
  or between two DBSCAN clusters, with configurable `p.adjust`, `cutoff`, and
  `medianDiffCutoff`.
* `summarizeCompoundFeatures()` and `summarizeClusterFeatures()` return
  per-group feature medians in long format for heatmaps.
* `exampleDataPath()` locates the bundled example tables.

## Documentation

* `vignette("getting-started")` was rewritten around the bundled example data
  and the package's own plotting functions, and now states the model behind
  each step (within-plate standardization, UMAP, per-compound DBSCAN,
  principal-curve trajectories with arc-length calibration, projection and
  interpolation, rank-based feature selection).
* New `vignette("shiny-app")` documents the companion application.

# cellpacman 0.0.1

* Initial version with the analysis API and an embedded Shiny interface.
