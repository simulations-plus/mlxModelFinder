test_that("smoke: default oral/clearance runs", {
  skip_if_no_monolix()

  out <- findModel(
    project = demo_proj,
    library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "clearance"),
    settings = tiny_settings()
  )

  expect_type(out, "list")
  expect_s3_class(out$results, "data.frame")
  expect_true(nrow(out$results) >= 1)
  expect_s3_class(out$best_model, "data.frame")
  expect_true(nrow(out$best_model) %in% c(0L, 1L))
})
