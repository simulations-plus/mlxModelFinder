test_that("custom metric_func runs and influences selection", {
  skip_if_no_monolix()

  # Trivial metric: always return 0 (lowest). Verifies signature + execution.
  my_metric <- function() 0

  out <- findModel(
    project = demo_proj, library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "clearance"),
    settings = tiny_settings(list(max_models = 2)),
    metric_func = my_metric
  )

  expect_true("metric" %in% names(out$best_model))
  expect_equal(out$best_model$metric, 0)
})

test_that("custom initial_func is called with project path and returns a path", {
  skip_if_no_monolix()

  called <- FALSE
  my_init <- function(p) {
    called <<- TRUE
    # do nothing; return same project path
    p
  }

  out <- findModel(
    project = demo_proj, library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "clearance"),
    settings = tiny_settings(list(max_models = 1)),
    initial_func = my_init
  )

  expect_true(called)
  expect_s3_class(out$results, "data.frame")
})
