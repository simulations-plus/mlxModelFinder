test_that("reproducible best_model with fixed seed", {
  skip_if_no_monolix()

  settings <- tiny_settings(list(max_models = 2))  # allow >1 model but still tiny

  run1 <- findModel(
    project = demo_proj,
    library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "clearance"),
    filters = tiny_pk_filters(),
    settings = settings
  )
  run2 <- findModel(
    project = demo_proj,
    library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "clearance"),
    filters = tiny_pk_filters(),
    settings = settings
  )

  # Compare best models by structural choices and metric
  expect_equal(names(run1$best_model), names(run2$best_model))
  if (nrow(run1$best_model) == 1 && nrow(run2$best_model) == 1) {
    expect_equal(run1$best_model[ , setdiff(names(run1$best_model), "metric"), drop = FALSE],
                 run2$best_model[ , setdiff(names(run2$best_model), "metric"), drop = FALSE])
    expect_equal(run1$best_model$metric, run2$best_model$metric, tolerance = 1e-8)
  }
})
