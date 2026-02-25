# tests/testthat/test-obs-mapping.R
# Tests for the obs_mapping feature

# --- Validation tests (no Monolix needed) ---

test_that("validate_obs_mapping rejects non-list input", {
  expect_error(
    validate_obs_mapping("mapping"),
    "non-empty list"
  )
  expect_error(
    validate_obs_mapping(42),
    "non-empty list"
  )
})

test_that("validate_obs_mapping rejects empty list", {
  expect_error(
    validate_obs_mapping(list()),
    "non-empty list"
  )
})

test_that("validate_obs_mapping rejects non-list entries", {
  expect_error(
    validate_obs_mapping(list("not a list")),
    "must be a list"
  )
})

test_that("validate_obs_mapping rejects missing obsId", {
  expect_error(
    validate_obs_mapping(list(list(modelOutput = "Cc"))),
    "obsId must be a single character"
  )
})

test_that("validate_obs_mapping rejects missing modelOutput", {
  expect_error(
    validate_obs_mapping(list(list(obsId = "1"))),
    "modelOutput must be a single character"
  )
})

test_that("validate_obs_mapping rejects non-character obsId", {
  expect_error(
    validate_obs_mapping(list(list(obsId = 1, modelOutput = "Cc"))),
    "obsId must be a single character"
  )
})

test_that("validate_obs_mapping rejects non-character modelOutput", {
  expect_error(
    validate_obs_mapping(list(list(obsId = "1", modelOutput = 42))),
    "modelOutput must be a single character"
  )
})

test_that("validate_obs_mapping rejects non-character observationName", {
  expect_error(
    validate_obs_mapping(list(list(obsId = "1", modelOutput = "Cc", observationName = 42))),
    "observationName must be a single character"
  )
})

test_that("validate_obs_mapping accepts valid single mapping", {
  expect_true(validate_obs_mapping(
    list(list(obsId = "1", modelOutput = "Cc"))
  ))
})

test_that("validate_obs_mapping accepts valid mapping with observationName", {
  expect_true(validate_obs_mapping(
    list(list(obsId = "1", modelOutput = "Cc", observationName = "concentration"))
  ))
})

test_that("validate_obs_mapping accepts valid multiple mappings", {
  expect_true(validate_obs_mapping(list(
    list(obsId = "1", modelOutput = "Cc", observationName = "concentration"),
    list(obsId = "2", modelOutput = "Level")
  )))
})

test_that("validate_obs_mapping validates each entry independently", {
  # First entry valid, second missing modelOutput
  expect_error(
    validate_obs_mapping(list(
      list(obsId = "1", modelOutput = "Cc"),
      list(obsId = "2")
    )),
    "obs_mapping\\[\\[2\\]].*modelOutput"
  )
})

# --- validate_obsidtouse skips for custom + obs_mapping ---

test_that("validate_obsidtouse skips when library=custom and obs_mapping provided", {
  # Should return TRUE without needing Monolix (no getData call)
  expect_true(validate_obsidtouse(
    settings = list(),
    library = "custom",
    obs_mapping = list(list(obsId = "1", modelOutput = "Cc"))
  ))
})

# --- Integration tests (require Monolix) ---

test_that("obs_mapping works with custom model exhaustive search", {
  skip_if_no_monolix()

  custom_func <- create_simple_custom_func()

  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = custom_func,
    obs_mapping = list(list(obsId = "concentration", modelOutput = "Cc")),
    settings = list(
      iiv = FALSE,
      error = FALSE,
      plot = FALSE,
      save_mode = "none"
    )
  )

  expect_type(result, "list")
  expect_true(nrow(result$results) > 0)

  # Verify the mapping was applied
  mapping <- lixoftConnectors::getMapping()
  expect_equal(mapping$mapping[[1]]$modelOutput, "Cc")
})

test_that("obs_mapping works with IIV and custom model exhaustive search", {
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
    obs_mapping = list(list(obsId = "concentration", modelOutput = "Cc")),
    settings = list(
      iiv = TRUE,
      error = FALSE,
      plot = FALSE,
      save_mode = "none"
    )
  )

  expect_type(result, "list")
  # 2^3 = 8 IIV combinations
  expect_equal(nrow(result$results), 8)
})

test_that("findModel works without obs_mapping (NULL default, backward compat)", {
  skip_if_no_monolix()

  custom_func <- create_simple_custom_func()

  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = custom_func,
    settings = list(
      iiv = FALSE,
      error = FALSE,
      plot = FALSE,
      save_mode = "none"
    )
  )

  expect_type(result, "list")
  expect_true(nrow(result$results) > 0)
})
