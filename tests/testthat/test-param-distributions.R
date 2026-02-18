# tests/testthat/test-param-distributions.R
# Tests for the param_distributions feature

# demo_proj is defined in helper-monolix.R

# --- Validation tests (no Monolix needed) ---

test_that("validate_param_distributions rejects non-list input", {
  expect_error(
    validate_param_distributions("logNormal"),
    "must be a named list"
  )
})

test_that("validate_param_distributions rejects unnamed list", {
  expect_error(
    validate_param_distributions(list("logNormal", "normal")),
    "must be a named list"
  )
})

test_that("validate_param_distributions rejects invalid distribution values", {
  expect_error(
    validate_param_distributions(list(V = "lognormal")),
    "Invalid distribution"
  )
  expect_error(
    validate_param_distributions(list(V = "logNormal", Cl = "foo")),
    "Invalid distribution"
  )
})

test_that("validate_param_distributions accepts all valid values", {
  expect_true(validate_param_distributions(list(V = "normal")))
  expect_true(validate_param_distributions(list(V = "logNormal")))
  expect_true(validate_param_distributions(list(V = "logitNormal")))
  expect_true(validate_param_distributions(
    list(V = "logNormal", Cl = "normal", ka = "logitNormal")
  ))
})

# --- Integration tests (require Monolix) ---

test_that("param_distributions works with exhaustive search", {
  skip_if_no_monolix()

  custom_func <- create_simple_custom_func()

  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = custom_func,
    param_distributions = list(V = "logNormal", Cl = "normal", ka = "logNormal"),
    settings = list(
      iiv = FALSE,
      error = FALSE,
      plot = FALSE,
      save_mode = "none"
    )
  )

  expect_type(result, "list")
  expect_true(nrow(result$results) > 0)

  # Verify distributions were applied
  dists <- lixoftConnectors::getIndividualParameterModel()$distribution
  expect_equal(unname(dists["Cl"]), "normal")
  expect_equal(unname(dists["V"]), "logNormal")
})

test_that("param_distributions works with IIV and exhaustive search", {
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
    param_distributions = list(V = "normal", Cl = "normal"),
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

  # Verify distributions were applied on the last model run
  dists <- lixoftConnectors::getIndividualParameterModel()$distribution
  expect_equal(unname(dists["V"]), "normal")
  expect_equal(unname(dists["Cl"]), "normal")
})

test_that("param_distributions ignores params not in current model", {
  skip_if_no_monolix()

  custom_func <- create_simple_custom_func()

  # Include a parameter (Tk0) that doesn't exist in the model (V, Cl, ka)
  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = custom_func,
    param_distributions = list(V = "normal", Tk0 = "logitNormal"),
    settings = list(
      iiv = FALSE,
      error = FALSE,
      plot = FALSE,
      save_mode = "none"
    )
  )

  expect_type(result, "list")
  expect_true(nrow(result$results) > 0)

  # V should have been set, Tk0 silently ignored
  dists <- lixoftConnectors::getIndividualParameterModel()$distribution
  expect_equal(unname(dists["V"]), "normal")
})

test_that("findModel works without param_distributions (NULL default)", {
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

test_that("param_distributions works with pk library and decision tree", {
  skip_if_no_monolix()

  result <- findModel(
    project = demo_proj,
    library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "clearance"),
    filters = list(
      absorption = "firstOrder",
      distribution = "1compartment",
      elimination = "linear"
    ),
    algorithm = "decision_tree",
    param_distributions = list(V = "normal"),
    settings = list(
      iiv = FALSE,
      error = FALSE,
      plot = FALSE,
      save_mode = "none"
    )
  )

  expect_type(result, "list")
  expect_true(nrow(result$results) > 0)

  dists <- lixoftConnectors::getIndividualParameterModel()$distribution
  expect_equal(unname(dists["V"]), "normal")
})
