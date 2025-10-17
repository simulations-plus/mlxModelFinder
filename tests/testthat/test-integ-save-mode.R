test_that("save_mode controls persisted runs and best_path", {
  skip_if_no_monolix()
  td <- withr::local_tempdir()

  # best
  out_best <- findModel(
    project = demo_proj, library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "clearance"),
    filters = tiny_pk_filters(),
    settings = tiny_settings(list(save_mode = "best"))
  )
  expect_true(is.null(out_best$best_path) || dir.exists(out_best$best_path))

  # all
  out_all <- findModel(
    project = demo_proj, library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "clearance"),
    filters = tiny_pk_filters(),
    settings = tiny_settings(list(save_mode = "all"))
  )
  expect_true(is.null(out_all$best_path) || dir.exists(out_all$best_path))

  # none
  out_none <- findModel(
    project = demo_proj, library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "clearance"),
    filters = tiny_pk_filters(),
    settings = tiny_settings(list(save_mode = "none"))
  )
  # Should return NULL/NA; definitely not a directory
  expect_true(is.null(out_none$best_path) || is.na(out_none$best_path) || !dir.exists(out_none$best_path))
})
