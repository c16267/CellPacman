# Application construction and runtime configuration.
#
# cellpacman_app() assembles the shiny.appobj from the UI, server, and workflow
# modules defined in the other R/ files. app.R calls it with settings read by
# app_config(); tests call it directly.

# Read runtime settings from environment variables (see app.R).
app_config <- function(env = Sys.getenv) {
  use_example_data <- tolower(env("CELLPACMAN_EXAMPLE_DATA", "false")) %in% c("true", "1", "yes")
  workers <- suppressWarnings(as.integer(env("CELLPACMAN_WORKERS", "1")))
  if (is.na(workers) || workers < 1L) {
    workers <- 1L
  }
  list(use_example_data = use_example_data, workers = workers)
}

#' Build the PACMAN Shiny application
#'
#' @param use_example_data Whether to pre-load the example feature and
#'   metadata tables bundled with the cellpacman package
#'   (`cellpacman::exampleDataPath()`). Uploaded files replace either table.
#' @param workers Number of background R processes. File parsing, validation,
#'   and every analysis stage run in these workers through a
#'   `future::multisession` plan so that the interface stays responsive. The
#'   plan is installed when the app starts and the previous plan is restored
#'   when it stops.
#' @param ... Additional arguments passed to `shiny::shinyApp()`, e.g.
#'   `options = list(port = 3838, host = "0.0.0.0")`.
#'
#' @return A `shiny.appobj`; print it or pass it to `shiny::runApp()`.
cellpacman_app <- function(use_example_data = FALSE, workers = 1L, ...) {
  if (!is.logical(use_example_data) || length(use_example_data) != 1 || is.na(use_example_data)) {
    stop("`use_example_data` must be TRUE or FALSE.")
  }
  if (!is.numeric(workers) || length(workers) != 1 || is.na(workers) || workers < 1 || workers != as.integer(workers)) {
    stop("`workers` must be a single positive integer.")
  }
  workers <- as.integer(workers)

  previous_plan <- NULL
  shiny::shinyApp(
    ui = app_ui(use_example_data = use_example_data),
    server = function(input, output, session) {
      app_server(input, output, session, use_example_data = use_example_data)
    },
    onStart = function() {
      # Uploaded feature tables are large (tens of MB for a 384-well screen).
      options(shiny.maxRequestSize = 100 * 1024^2)
      previous_plan <<- future::plan()
      # I(1L) forces a separate process even on a single-core allocation.
      future::plan(future::multisession, workers = I(workers))
      shiny::onStop(function() future::plan(previous_plan))
    },
    ...
  )
}

# Directory holding the static assets (www/cellpacman.css, www/cellpacman.js).
# shiny::runApp() sets the working directory to the app directory, so the
# default resolves there; a different host (e.g. the test app under
# tests/testthat/apps) sets options(cellpacman.app.www = ...).
app_www_dir <- function() {
  dir <- getOption("cellpacman.app.www", "www")
  if (!dir.exists(dir)) {
    stop(
      "Static assets directory `", dir, "` was not found. Run the app from its ",
      "own directory or set `options(cellpacman.app.www = \"<path>/www\")`."
    )
  }
  normalizePath(dir)
}
