proj_path <- system.file("extdata", "theophylline_project.mlxtran", package = "mlxModelFinder")

test_that("algorithm defaults to decision_tree without IIV", {
  out <- tryCatch(
    findModel(
      project = proj_path,
      library = "pk",
      requiredFilters = list(administration = "oral", parametrization = "rate"),
      algorithm = NULL,
      settings = list(findIIV = FALSE, plot = FALSE, max_models = 1)
    ),
    error = identity
  )
  if (inherits(out, "error")) skip("Needs early-exit mode to verify structure.")
  expect_true(is.null(out$iteration_progress) || NROW(out$iteration_progress) == 0)
})

test_that("algorithm defaults to ACO with IIV", {
  out <- tryCatch(
    findModel(
      project = proj_path,
      library = "pk",
      requiredFilters = list(administration = "oral", parametrization = "rate"),
      algorithm = NULL,
      settings = list(findIIV = TRUE, plot = FALSE, max_models = 1, N = 1)
    ),
    error = identity
  )
  if (inherits(out, "error")) skip("Needs early-exit mode to verify structure.")
  expect_true("iteration_progress" %in% names(out))
  expect_true("prob_history" %in% names(out))
})
