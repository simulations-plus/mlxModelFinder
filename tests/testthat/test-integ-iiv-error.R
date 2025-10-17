test_that("settings$findIIV toggles algorithm and returns ACO fields", {
  skip_if_no_monolix()

  # No IIV -> decision_tree by default, no ACO fields
  out1 <- findModel(
    project = demo_proj, library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "clearance"),
    settings = tiny_settings(list(findIIV = FALSE))
  )
  expect_true(is.null(out1$iteration_progress) || NROW(out1$iteration_progress) == 0)
  expect_true(is.null(out1$prob_history) || length(out1$prob_history) == 0)

  # With IIV -> ACO by default, ACO fields present (even if small)
  out2 <- findModel(
    project = demo_proj, library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "clearance"),
    settings = tiny_settings(list(findIIV = TRUE, max_models = 2))
  )
  expect_true("iteration_progress" %in% names(out2))
  expect_true("prob_history" %in% names(out2))
})

test_that("findError toggles residual error exploration without breaking runs", {
  skip_if_no_monolix()

  base <- list(administration = "oral", parametrization = "clearance")

  off <- findModel(
    project = demo_proj, library = "pk",
    requiredFilters = base,
    settings = tiny_settings(list(findError = FALSE))
  )
  on <- findModel(
    project = demo_proj, library = "pk",
    requiredFilters = base,
    settings = tiny_settings(list(findError = TRUE))
  )

  expect_s3_class(off$results, "data.frame")
  expect_s3_class(on$results, "data.frame")
})
