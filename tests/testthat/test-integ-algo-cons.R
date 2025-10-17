test_that("exhaustive and decision_tree agree on tiny space", {
  skip_if_no_monolix()

  req <- list(administration = "oral", parametrization = "clearance")
  filters <- tiny_pk_filters()

  ex <- findModel(
    project = demo_proj, library = "pk",
    requiredFilters = req, filters = filters,
    algorithm = "exhaustive_search",
    settings  = tiny_settings(list(max_models = 100))  # enough to cover tiny space
  )
  dt <- findModel(
    project = demo_proj, library = "pk",
    requiredFilters = req, filters = filters,
    algorithm = "decision_tree",
    settings  = tiny_settings(list(max_models = 10))
  )

  # Compare best structural choices (ignore metric floating jitter)
  common <- intersect(names(ex$best_model), names(dt$best_model))
  common <- setdiff(common, "metric")
  expect_equal(ex$best_model[ , common, drop = FALSE],
               dt$best_model[ , common, drop = FALSE],
               ignore_attr = TRUE)
})
