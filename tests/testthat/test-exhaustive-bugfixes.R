# tests/testthat/test-exhaustive-bugfixes.R
# Tests for exhaustive search bugfixes:
# 1. fixed_iiv applied when iiv = FALSE
# 2. save_mode = "all" saves individual run files

# demo_proj is defined in helper-monolix.R

# --- fixed_iiv with iiv = FALSE (integration tests) ---

test_that("fixed_iiv=FALSE turns off IIV when iiv=FALSE in exhaustive search", {
  skip_if_no_monolix()

  custom_func <- create_simple_custom_func()
  param_mapping <- list(
    common_params = c("V", "Cl", "ka"),
    filter_params = list(
      absorption = list(firstOrder = character(0))
    )
  )

  # Use save_mode = "best" so the best model is re-run on the main session
  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = custom_func,
    param_mapping = param_mapping,
    fixed_iiv = list(ka = FALSE),
    settings = list(
      iiv = FALSE,
      error = FALSE,
      plot = FALSE,
      save_mode = "best"
    )
  )

  expect_equal(nrow(result$results), 1)

  # After findModel, the main session has the best model loaded with correct IIV
  var_state <- lixoftConnectors::getIndividualParameterModel()$variability$id
  expect_false(var_state[["ka"]])
  # V and Cl should still have IIV (default is TRUE)
  expect_true(var_state[["V"]])
  expect_true(var_state[["Cl"]])
})

test_that("fixed_iiv=TRUE forces IIV on when iiv=FALSE in exhaustive search", {
  skip_if_no_monolix()

  custom_func <- create_simple_custom_func()
  param_mapping <- list(
    common_params = c("V", "Cl", "ka"),
    filter_params = list(
      absorption = list(firstOrder = character(0))
    )
  )

  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = custom_func,
    param_mapping = param_mapping,
    fixed_iiv = list(V = TRUE, Cl = TRUE, ka = TRUE),
    settings = list(
      iiv = FALSE,
      error = FALSE,
      plot = FALSE,
      save_mode = "best"
    )
  )

  expect_equal(nrow(result$results), 1)

  # All params should have IIV forced on
  var_state <- lixoftConnectors::getIndividualParameterModel()$variability$id
  expect_true(var_state[["V"]])
  expect_true(var_state[["Cl"]])
  expect_true(var_state[["ka"]])
})

test_that("fixed_iiv with iiv=FALSE works for pk library exhaustive search", {
  skip_if_no_monolix()

  result <- findModel(
    project = demo_proj,
    library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "clearance"),
    filters = list(
      absorption = "firstOrder",
      distribution = "1compartment",
      elimination = "linear",
      delay = "noDelay"
    ),
    algorithm = "exhaustive_search",
    fixed_iiv = list(ka = FALSE),
    settings = list(
      iiv = FALSE,
      error = FALSE,
      plot = FALSE,
      save_mode = "best"
    )
  )

  expect_equal(nrow(result$results), 1)

  # ka should have IIV turned off
  var_state <- lixoftConnectors::getIndividualParameterModel()$variability$id
  expect_false(var_state[["ka"]])
})

# --- save_mode = "all" tests ---

test_that("save_mode=all creates individual run files in exhaustive search", {
  skip_if_no_monolix()

  custom_func <- create_simple_custom_func()
  param_mapping <- list(
    common_params = c("V", "Cl", "ka"),
    filter_params = list(
      absorption = list(firstOrder = character(0))
    )
  )

  # Get result_folder BEFORE findModel changes the project state
  lixoftConnectors::loadProject(demo_proj)
  result_folder <- lixoftConnectors::getProjectSettings()$directory
  autobuild_dir <- file.path(result_folder, "autoBuild")

  # Clean up any previous autoBuild contents
  if (dir.exists(autobuild_dir)) {
    old_runs <- list.files(autobuild_dir, pattern = "^run_.*\\.mlxtran$", full.names = TRUE)
    file.remove(old_runs)
  }

  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = custom_func,
    param_mapping = param_mapping,
    settings = list(
      iiv = FALSE,
      error = FALSE,
      plot = FALSE,
      save_mode = "all"
    )
  )

  expect_equal(nrow(result$results), 1)

  # Check that run files were created
  run_files <- list.files(autobuild_dir, pattern = "^run_.*\\.mlxtran$")
  expect_true(length(run_files) >= 1)

  # best_run should also exist (post-processing saves it for "all" mode)
  expect_true(file.exists(file.path(autobuild_dir, "best_run.mlxtran")))
})

test_that("save_mode=all creates correct number of run files with IIV search", {
  skip_if_no_monolix()

  custom_func <- create_simple_custom_func()
  param_mapping <- list(
    common_params = c("V", "Cl", "ka"),
    filter_params = list(
      absorption = list(firstOrder = character(0))
    )
  )

  # Get result_folder BEFORE findModel changes the project state
  lixoftConnectors::loadProject(demo_proj)
  result_folder <- lixoftConnectors::getProjectSettings()$directory
  autobuild_dir <- file.path(result_folder, "autoBuild")

  # Clean up any previous autoBuild contents
  if (dir.exists(autobuild_dir)) {
    old_runs <- list.files(autobuild_dir, pattern = "^run_.*\\.mlxtran$", full.names = TRUE)
    file.remove(old_runs)
  }

  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = custom_func,
    param_mapping = param_mapping,
    fixed_iiv = list(V = TRUE, Cl = TRUE),
    settings = list(
      iiv = TRUE,
      error = FALSE,
      plot = FALSE,
      save_mode = "all"
    )
  )

  # V and Cl fixed, ka free -> 2^1 = 2 combos
  expect_equal(nrow(result$results), 2)

  # Should have 2 individual run files
  run_files <- list.files(autobuild_dir, pattern = "^run_.*\\.mlxtran$")
  expect_equal(length(run_files), 2)
})
