# PACMAN (PAinted Cell and coMpound ANalysis) -- Shiny application entry point.
#
# This directory is a standard Shiny app that drives the `cellpacman` R
# package: every analysis step calls the package's exported API, and the files
# in R/ (sourced automatically by shiny::runApp()) only orchestrate the
# workflow and render results.
#
#   shiny::runApp("cellpacman-shiny")                    # from the project root
#   shiny::runApp("cellpacman-shiny", port = 3838)
#
# Configuration is read from environment variables so that the same file works
# unchanged on a laptop, Shiny Server, Posit Connect, or Docker:
#
#   CELLPACMAN_EXAMPLE_DATA  "true" to pre-load the example tables bundled with
#                            the cellpacman package (default "false")
#   CELLPACMAN_WORKERS       number of background R processes used for parsing
#                            and analysis (default 1)
#
# The *installed* cellpacman package must be available to the R process that
# serves the app, because the background workers load the installed namespace.

library(shiny)
library(cellpacman)

config <- app_config()

cellpacman_app(
  use_example_data = config$use_example_data,
  workers = config$workers
)
