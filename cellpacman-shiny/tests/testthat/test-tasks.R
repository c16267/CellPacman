test_that("background tasks leave the session responsive and cache completed work", {
  old_plan <- future::plan(future::multisession, workers = I(1L))
  on.exit(future::plan(old_plan), add = TRUE)
  gate <- tempfile()
  on.exit(unlink(gate), add = TRUE)

  shiny::testServer(function(input, output, session) {
    task <- app_task(function(gate) {
      deadline <- Sys.time() + 15
      while (!file.exists(gate) && Sys.time() < deadline) Sys.sleep(0.02)
      list(pid = Sys.getpid(), finished = Sys.time())
    })
    output$ping <- shiny::renderText(input$ping)
  }, {
    task$invoke(list(gate = gate), key = "same-input")
    session$flushReact()
    expect_identical(task$status(), "running")
    session$setInputs(ping = "responsive")
    expect_identical(output$ping, "responsive")
    expect_identical(task$status(), "running")
    file.create(gate)
    wait_for_app(function() task$status() == "ready", session)
    first <- task$result()
    expect_false(identical(first$pid, Sys.getpid()))

    task$invoke(list(gate = gate), key = "same-input")
    expect_identical(task$status(), "ready")
    expect_identical(task$result(), first)
  })
})

test_that("obsolete results are discarded and pending requests are coalesced", {
  old_plan <- future::plan(future::multisession, workers = I(1L))
  on.exit(future::plan(old_plan), add = TRUE)
  gate <- tempfile()
  on.exit(unlink(gate), add = TRUE)

  shiny::testServer(function(input, output, session) {
    task <- app_task(function(value, gate) {
      deadline <- Sys.time() + 15
      while (!file.exists(gate) && Sys.time() < deadline) Sys.sleep(0.02)
      value
    })
    seen <- shiny::reactiveVal(list())
    shiny::observeEvent(task$result(), { seen(c(seen(), list(task$result()))) })
  }, {
    task$invoke(list(value = "old", gate = gate), key = 1)
    session$flushReact()
    task$invalidate()
    expect_identical(task$status(), "empty")
    task$invoke(list(value = "superseded", gate = gate), key = 2)
    task$invoke(list(value = "latest", gate = gate), key = 3)
    file.create(gate)
    wait_for_app(function() task$status() == "ready", session)
    expect_identical(task$result(), "latest")
    expect_identical(seen(), list("latest"))
  })
})

test_that("task errors are visible and a later request can recover", {
  shiny::testServer(function(input, output, session) {
    task <- app_task(function(fail) {
      if (fail) stop("Analysis failed")
      "recovered"
    })
  }, {
    task$invoke(list(fail = TRUE), key = 1)
    wait_for_app(function() task$status() == "error", session)
    expect_error(task$result(), "Analysis failed")
    task$invoke(list(fail = FALSE), key = 2)
    wait_for_app(function() task$status() == "ready", session)
    expect_identical(task$result(), "recovered")
  })
})
