# Keep the worker closure separate from session state: only ordinary arguments
# and the analysis function should be serialized to the background process.
app_task_worker <- function(fun) {
  force(fun)
  function(id, args, key) {
    rng_kind <- RNGkind()
    promises::future_promise({
      # Preserve the RNG algorithm used by the synchronous API. UMAP sets a
      # fixed seed, whose result also depends on RNGkind(), not just the seed.
      do.call(RNGkind, as.list(rng_kind))
      tryCatch(
        list(id = id, key = key, value = do.call(fun, args), error = NULL),
        error = function(e) list(id = id, key = key, value = NULL, error = conditionMessage(e))
      )
    }, seed = TRUE)
  }
}

# One bounded, session-local cache entry per stage. Replacing a running request
# retains only the latest pending request, rather than queueing obsolete uploads.
app_task <- function(fun, session = shiny::getDefaultReactiveDomain()) {
  task <- shiny::ExtendedTask$new(app_task_worker(fun))
  state <- shiny::reactiveVal(list(status = "empty", value = NULL, error = NULL))
  revision <- 0L
  cache <- NULL
  pending <- NULL
  running <- FALSE
  closed <- FALSE

  start <- function(request) {
    running <<- TRUE
    task$invoke(request$id, request$args, request$key)
  }

  shiny::observe({
    status <- task$status()
    if (!status %in% c("success", "error")) return()
    result <- tryCatch(task$result(), error = function(e) {
      list(id = NA_integer_, error = conditionMessage(e))
    })
    running <<- FALSE
    if (closed) return()
    if (!is.null(pending)) {
      request <- pending
      pending <<- NULL
      start(request)
      return()
    }
    if (identical(result$id, revision)) {
      if (is.null(result$error)) {
        cache <<- list(key = result$key, value = result$value)
        state(list(status = "ready", value = result$value, error = NULL))
      } else {
        state(list(status = "error", value = NULL, error = result$error))
      }
    } else if (identical(status, "error") && identical(shiny::isolate(state())$status, "running")) {
      state(list(status = "error", value = NULL, error = result$error))
    }
  })

  session$onSessionEnded(function() {
    closed <<- TRUE
    pending <<- NULL
    cache <<- NULL
  })

  list(
    invoke = function(args, key) {
      revision <<- revision + 1L
      pending <<- NULL
      if (!is.null(cache) && identical(cache$key, key)) {
        state(list(status = "ready", value = cache$value, error = NULL))
      } else {
        state(list(status = "running", value = NULL, error = NULL))
        request <- list(id = revision, args = args, key = key)
        if (running) pending <<- request else start(request)
      }
      invisible(NULL)
    },
    invalidate = function() {
      revision <<- revision + 1L
      pending <<- NULL
      cache <<- NULL
      state(list(status = "empty", value = NULL, error = NULL))
    },
    status = shiny::reactive(state()$status),
    result = shiny::reactive({
      current <- state()
      shiny::req(current$status != "running", cancelOutput = "progress")
      shiny::validate(shiny::need(current$status != "error", current$error))
      shiny::req(current$status == "ready")
      current$value
    })
  )
}
