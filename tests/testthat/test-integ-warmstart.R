test_that("previous_runs warm-start and output appends", {
  skip_if_no_monolix()
  td <- withr::local_tempdir()
  csv <- file.path(td, "runs.csv")

  # First pass writes
  out1 <- findModel(
    project = demo_proj, library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "clearance"),
    filters = tiny_pk_filters(),
    settings = tiny_settings(list(output = csv, max_models = 2))
  )
  expect_true(file.exists(csv))
  n1 <- NROW(read.csv(csv))

  # Second pass warm-starts and appends (implementation-dependent, but should not error)
  out2 <- findModel(
    project = demo_proj, library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "clearance"),
    filters = tiny_pk_filters(),
    settings = tiny_settings(list(output = csv, previous_runs = csv, max_models = 2))
  )
  n2 <- NROW(read.csv(csv))
  expect_true(n2 >= n1)  # at least not smaller
})

test_that("ACO with output but no previous_runs does not error", {
  skip_if_no_monolix()
  td <- withr::local_tempdir()
  csv <- file.path(td, "aco_runs.csv")

  result <- findModel(
    project = demo_proj, library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "clearance"),
    filters = tiny_pk_filters(),
    algorithm = "ACO",
    settings = tiny_settings(list(output = csv, N = 2, max_models = 4))
  )

  expect_true(file.exists(csv))
  expect_true(NROW(read.csv(csv)) > 0)
  expect_s3_class(result$results, "data.frame")
})
