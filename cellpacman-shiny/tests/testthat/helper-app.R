# Load the application's R/ directory into the test environment, exactly as
# shiny::runApp() does before sourcing app.R. testthat runs with the working
# directory set to tests/testthat, so the app root is two levels up.
library(shiny)
library(cellpacman)

app_dir <- normalizePath(file.path("..", ".."), mustWork = TRUE)
options(cellpacman.app.www = file.path(app_dir, "www"))
shiny::loadSupport(app_dir, renv = environment(), globalrenv = NULL)
