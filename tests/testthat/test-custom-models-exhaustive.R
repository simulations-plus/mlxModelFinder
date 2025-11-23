# tests/testthat/test-custom-models-exhaustive.R

# demo_proj is defined in helper-monolix.R

test_that("exhaustive search works with custom models (no IIV, no error)", {
  skip_if_no_monolix()
  
  custom_func <- create_simple_custom_func()
  
  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = c("firstOrder", "zeroOrder")),
    algorithm = "exhaustive_search",
    model_creation_func = custom_func,
    settings = list(
      iiv = FALSE,
      error = FALSE,
      plot = FALSE,
      max_models = 2,
      save_mode = "none"
    )
  )
  
  expect_type(result, "list")
  expect_true("results" %in% names(result))
  expect_true("best_model" %in% names(result))
  expect_s3_class(result$results, "data.frame")
  expect_true(nrow(result$results) > 0)
  expect_true("metric" %in% names(result$results))
})

test_that("exhaustive search works with custom models and IIV", {
  skip_if_no_monolix()
  
  custom_func <- create_simple_custom_func()
  param_mapping <- create_simple_param_mapping()
  
  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = custom_func,
    param_mapping = param_mapping,
    settings = list(
      iiv = TRUE,
      error = FALSE,
      plot = FALSE,
      max_models = 4,  # Will test multiple IIV combinations
      save_mode = "none"
    )
  )
  
  expect_type(result, "list")
  expect_true("results" %in% names(result))
  expect_s3_class(result$results, "data.frame")
  
  # Check that IIV parameters are in the results
  param_cols <- intersect(c("V", "Cl", "ka"), names(result$results))
  expect_true(length(param_cols) > 0)
})

test_that("exhaustive search works with custom models and error models", {
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
      error = TRUE,
      plot = FALSE,
      max_models = 4,  # Will test multiple error models
      save_mode = "none"
    )
  )
  
  expect_type(result, "list")
  expect_true("results" %in% names(result))
  expect_s3_class(result$results, "data.frame")
  expect_true("error" %in% names(result$results))
})

test_that("exhaustive search with custom models and both IIV and error", {
  skip_if_no_monolix()
  
  custom_func <- create_simple_custom_func()
  param_mapping <- create_simple_param_mapping()
  
  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = custom_func,
    param_mapping = param_mapping,
    settings = list(
      iiv = TRUE,
      error = TRUE,
      plot = FALSE,
      max_models = 8,
      save_mode = "none"
    )
  )
  
  expect_type(result, "list")
  expect_true("results" %in% names(result))
  expect_s3_class(result$results, "data.frame")
  
  # Should have both IIV params and error column
  expect_true("error" %in% names(result$results))
  param_cols <- intersect(c("V", "Cl", "ka"), names(result$results))
  expect_true(length(param_cols) > 0)
})

test_that("exhaustive search with filter-based custom models", {
  skip_if_no_monolix()
  
  custom_func <- create_filter_based_custom_func()
  
  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = c("firstOrder", "zeroOrder")),
    algorithm = "exhaustive_search",
    model_creation_func = custom_func,
    settings = list(
      iiv = FALSE,
      error = FALSE,
      plot = FALSE,
      max_models = 2,
      save_mode = "none"
    )
  )
  
  expect_type(result, "list")
  expect_true("results" %in% names(result))
  expect_s3_class(result$results, "data.frame")
  expect_true("absorption" %in% names(result$results))
  
  # Should have tested both absorption types
  expect_true(any(result$results$absorption == "firstOrder"))
  expect_true(any(result$results$absorption == "zeroOrder"))
})

test_that("exhaustive search respects max_models limit with custom models", {
  skip_if_no_monolix()
  
  custom_func <- create_simple_custom_func()
  param_mapping <- create_simple_param_mapping()
  
  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = c("firstOrder", "zeroOrder")),
    algorithm = "exhaustive_search",
    model_creation_func = custom_func,
    param_mapping = param_mapping,
    settings = list(
      iiv = TRUE,
      error = FALSE,
      plot = FALSE,
      max_models = 3,  # Limit to 3 models
      save_mode = "none"
    )
  )
  
  expect_type(result, "list")
  expect_true(nrow(result$results) <= 3)
})

test_that("exhaustive search with custom models returns best model", {
  skip_if_no_monolix()
  
  custom_func <- create_simple_custom_func()
  
  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = c("firstOrder", "zeroOrder")),
    algorithm = "exhaustive_search",
    model_creation_func = custom_func,
    settings = list(
      iiv = FALSE,
      error = FALSE,
      plot = FALSE,
      max_models = 2,
      save_mode = "none"
    )
  )
  
  expect_true("best_model" %in% names(result))
  expect_s3_class(result$best_model, "data.frame")
  expect_equal(nrow(result$best_model), 1)
  
  # Best model should have the minimum metric
  expect_equal(
    result$best_model$metric,
    min(result$results$metric)
  )
})

test_that("exhaustive search with custom models handles single filter value", {
  skip_if_no_monolix()
  
  custom_func <- create_simple_custom_func()
  
  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),  # Single value
    algorithm = "exhaustive_search",
    model_creation_func = custom_func,
    settings = list(
      iiv = FALSE,
      error = FALSE,
      plot = FALSE,
      max_models = 1,
      save_mode = "none"
    )
  )
  
  expect_type(result, "list")
  expect_true(nrow(result$results) >= 1)
})
