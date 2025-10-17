# tests/testthat/test-findModel-input-validation.R
local_uses_testthat3 <- function() {
  if (utils::packageVersion("testthat") < "3.1.0") {
    testthat::skip("Requires testthat 3e for snapshotting")
  }
}

proj_path <- system.file("extdata", "theophylline_project.mlxtran", package = "mlxModelFinder")

test_that("errors if required args are missing", {
  expect_snapshot(error = TRUE, findModel())
  expect_snapshot(error = TRUE, findModel(project = NULL))
  expect_snapshot(error = TRUE, findModel(project = proj_path))
  expect_snapshot(error = TRUE, findModel(project = proj_path, library = "pk"))
})

test_that("errors on unsupported library", {
  expect_snapshot(error = TRUE, findModel(
    project = proj_path,
    library = "pd",
    requiredFilters = list(administration = "oral", parametrization = "clearance")
  ))
})

test_that("requiredFilters must include administration + parametrization with valid values", {
  expect_snapshot(error = TRUE, findModel(
    project = proj_path, library = "pk",
    requiredFilters = list()
  ))
  expect_snapshot(error = TRUE, findModel(
    project = proj_path, library = "pk",
    requiredFilters = list(administration = "oral")
  ))
  expect_snapshot(error = TRUE, findModel(
    project = proj_path, library = "pk",
    requiredFilters = list(parametrization = "rate")
  ))
  expect_snapshot(error = TRUE, findModel(
    project = proj_path, library = "pk",
    requiredFilters = list(administration = "iv_push", parametrization = "clearance")
  ))
  expect_snapshot(error = TRUE, findModel(
    project = proj_path, library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "foo")
  ))
})

test_that("optional filters validate known choices and bioavailability rule", {
  expect_snapshot(error = TRUE, findModel(
    project = proj_path, library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "rate"),
    filters = list(absorption = c("firstOrder", "bogus"))
  ))

  expect_snapshot(error = TRUE, findModel(
    project = proj_path, library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "rate"),
    filters = list(bioavailability = "true")
  ))

  expect_no_error(findModel(
    project = proj_path, library = "pk",
    requiredFilters = list(administration = "oralBolus", parametrization = "rate"),
    filters = list(bioavailability = "true"),
    algorithm = "decision_tree",
    settings = list(plot = FALSE, max_models = 1)
  ))
})

test_that("algorithm argument validates choices", {
  expect_snapshot(error = TRUE, findModel(
    project = proj_path, library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "rate"),
    algorithm = "genetic"
  ))

  expect_no_error(findModel(
    project = proj_path, library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "rate"),
    algorithm = "decision_tree",
    settings = list(plot = FALSE, max_models = 1)
  ))
})
