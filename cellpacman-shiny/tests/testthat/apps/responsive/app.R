# Exercise the real UI/workflow with predictable worker latency. Only this test
# app replaces the worker factory; the scientific functions are unchanged.
#
# shiny::runApp() sets the working directory to this directory, so the main
# application lives four levels up.
library(shiny)
library(cellpacman)

main_app <- normalizePath(file.path("..", "..", "..", ".."), mustWork = TRUE)
options(cellpacman.app.www = file.path(main_app, "www"))
shiny::loadSupport(main_app, renv = environment(), globalrenv = NULL)

original_worker <- app_task_worker
slow_analysis <- function(fun) {
  force(fun)
  function(...) {
    Sys.sleep(2)
    fun(...)
  }
}
# app_task() looks `app_task_worker` up in this environment at call time.
app_task_worker <- function(fun) original_worker(slow_analysis(fun))

app <- cellpacman_app()
server <- app$serverFuncSource()
shiny::shinyApp(
  ui = shiny::tagList(
    app_ui(),
    shiny::div(style = "display: none;", shiny::textOutput("ping_response"))
  ),
  server = function(input, output, session) {
    server(input, output, session)
    output$ping_response <- shiny::renderText(input$ping_request)
    shiny::outputOptions(output, "ping_response", suspendWhenHidden = FALSE)
  },
  onStart = app$onStart
)
