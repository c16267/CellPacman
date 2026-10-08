# cellpacman-shiny

Shiny front end for the [`cellpacman`](..) R package (PACMAN, PAinted Cell
and coMpound ANalysis). The application lives in this directory of the
package repository and is excluded from the package build (`.Rbuildignore`). The application is a standard Shiny app
directory that **depends on the installed `cellpacman` package**: it performs
no analysis of its own, and every stage calls the package's exported API.

```
cellpacman-shiny/
├── app.R                 entry point: reads the configuration, returns cellpacman_app()
├── R/                    sourced automatically by shiny::runApp()
│   ├── launch.R          app_config(), cellpacman_app(), app_www_dir()
│   ├── ui.R              app_ui(): bslib navbar page with the four tabs
│   ├── ui_helpers.R      plot/upload output containers
│   ├── server.R          app_server(): renders tables and plots from stage results
│   ├── server_helpers.R  upload parsing, analysis wrappers run in workers, volcano/heatmap plots
│   ├── workflow.R        app_workflow_server(): stage orchestration, caching, invalidation
│   └── tasks.R           app_task(): ExtendedTask wrapper around future_promise()
├── www/
│   ├── cellpacman.css    layout rules (tab grids, plot cards, disabled nav links)
│   └── cellpacman.js     `setNavDisabled` message handler used to gate the tabs
├── tests/
│   ├── testthat.R        runner
│   └── testthat/         server, workflow, task and browser tests (+ apps/responsive fixture)
├── Dockerfile            container image serving the app (built from the project root)
└── README.md
```

## Requirements

* R >= 4.1 with the `cellpacman` package installed
  (`remotes::install_github("c16267/cellpacman")` or
  `install.packages("cellpacman_<version>.tar.gz", repos = NULL, type = "source")`).
* Application packages: `shiny` (>= 1.8.1), `bslib` (>= 0.10.0), `htmltools`,
  `reactable`, `shinycssloaders` (>= 1.1.0), `plotly`, `ggplot2`, `dplyr`,
  `tidyr`, `future`, `promises`, `digest`; optional: `ggpubr` (plot theme),
  `ggrepel` (labels), `shinytest2` + `later` (browser test).

```r
install.packages(c("shiny", "bslib", "htmltools", "reactable", "shinycssloaders",
                   "plotly", "ggplot2", "dplyr", "tidyr", "future", "promises", "digest"))
```

## Running

```r
shiny::runApp("cellpacman-shiny")                 # from the repository root
shiny::runApp("cellpacman-shiny", port = 3838)
```

Configuration is read from environment variables so that the same directory
works unchanged on every host:

| Variable | Default | Effect |
| --- | --- | --- |
| `CELLPACMAN_EXAMPLE_DATA` | `false` | `true` pre-loads the example tables bundled with the package. |
| `CELLPACMAN_WORKERS` | `1` | Number of background R processes for parsing and analysis. |

```r
Sys.setenv(CELLPACMAN_EXAMPLE_DATA = "true", CELLPACMAN_WORKERS = "2")
shiny::runApp("cellpacman-shiny")
```

## How it works

* **Stages.** Data → Dimension Reduction → Clustering → Trajectory. Each tab
  is enabled by the server only after the previous stage has completed
  (`www/cellpacman.js` toggles the navbar links).
* **Background computation.** Uploads are parsed and validated in a
  background R process; dimension reduction, clustering and trajectory
  analysis also run in the background (`future::multisession` +
  `promises::future_promise()` inside `shiny::ExtendedTask`), so the
  interface stays responsive. Analysis buttons switch to a busy state in the
  browser before any R work starts.
* **Caching.** Each session keeps the most recent successful result per stage.
  Repeating a request with identical data and settings reuses it; replacing an
  input table invalidates every downstream stage; results of superseded jobs
  are discarded.
* **Workers load the installed package.** After editing `cellpacman`,
  reinstall it and restart the app so the workers use the new code. The app's
  own `R/` files are re-read on every `shiny::runApp()`.

Each interface element maps to one package function:

| Tab / panel | `cellpacman` function |
| --- | --- |
| Data: upload + preview | `readCellPainting()`, `loadCellPainting()` |
| Dimension Reduction: UMAP | `normalize()`, `dimReduce()`, `plotDimred()` |
| Dimension Reduction: Compare Compounds | `compareCompoundFeatures()`, `summarizeCompoundFeatures()` |
| Clustering: DBSCAN plot | `compoundCluster()`, `plotKnownClusters()` |
| Clustering: Compare Clusters | `compareClusterFeatures()`, `summarizeClusterFeatures()` |
| Trajectory: curve / projections | `curveEstimate()`, `curveProject()`, `plotKnownCurvesPlotly()`, `plotUnknownProjectionsPlotly()` |
| Trajectory: feature panels | `selectFeatureCurve()` |

## Deploying

The serving R process must be able to spawn child processes (workers) and
must have `cellpacman` installed.

* **Shiny Server / Posit Connect** — copy this directory to the server and
  install `cellpacman` plus the application packages in the server's library.
* **shinyapps.io / `rsconnect::deployApp()`** — rsconnect must be able to
  rebuild `cellpacman` on the server; install it with
  `remotes::install_github("c16267/cellpacman")` rather than from a local
  tarball before deploying.
* **Docker** — the image installs the package from the source tarball in
  the repository root and serves this app on port 3838. Build it from the
  **repository root** so that both are in the build context:

```sh
docker build -f cellpacman-shiny/Dockerfile -t cellpacman-shiny .
docker run --rm -p 3838:3838 -e CELLPACMAN_EXAMPLE_DATA=true cellpacman-shiny
```

## Tests

```r
# from this directory
testthat::test_dir("tests/testthat")   # server, workflow and task tests
shinytest2::test_app()                 # the same, plus the headless-Chrome browser test
```

`tests/testthat/helper-app.R` loads `R/` into the test environment exactly as
`shiny::runApp()` does; `tests/testthat/apps/responsive/` is a copy of the app
with artificially slow workers used by the browser test to verify that uploads
and analyses acknowledge work before it finishes.
