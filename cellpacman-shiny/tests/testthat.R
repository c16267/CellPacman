# Run the application tests from the app directory:
#
#   testthat::test_dir("tests/testthat")          # server, workflow, task tests
#   shinytest2::test_app()                        # the same, plus the browser test
#
# tests/testthat/helper-app.R loads the app's R/ files into the test
# environment; the browser test additionally needs shinytest2 and Chrome.
library(testthat)

testthat::test_dir(
  file.path("tests", "testthat"),
  reporter = testthat::default_reporter(),
  stop_on_failure = TRUE
)
