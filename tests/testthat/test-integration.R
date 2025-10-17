proj_path <- system.file("extdata", "theophylline_project.mlxtran", package = "mlxModelFinder")

test_that("integration: runs with real Monolix project", {
  skip_on_cran()
  skip_if_not_installed("lixoftConnectors")

  if (!file.exists(proj_path)) skip("Demo project not found.")

  out <- findModel(
    project = proj_path,
    library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "clearance"),
    settings = list(findIIV = FALSE, findError = FALSE, plot = FALSE, seed = 123456)
  )

  expect_type(out, "list")
  expect_s3_class(out$results, "data.frame")
  expect_true(nrow(out$results) >= 1)
  expect_s3_class(out$best_model, "data.frame")
  expect_true(nrow(out$best_model) == 1)
})
