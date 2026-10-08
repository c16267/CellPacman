wait_for_app <- function(predicate, session, timeout = 20) {
  deadline <- Sys.time() + timeout
  repeat {
    later::run_now(0.02)
    session$flushReact()
    if (isTRUE(shiny::isolate(predicate()))) return(invisible(NULL))
    if (Sys.time() > deadline) stop("Timed out waiting for app task")
  }
}
