# tests/testthat/test-custom-models-aco.R

# demo_proj is defined in helper-monolix.R

test_that("ACO works with custom models (no IIV, no error)", {
  skip_if_no_monolix()
  
  custom_func <- create_simple_custom_func()
  
  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = c("firstOrder", "zeroOrder")),
    algorithm = "ACO",
    model_creation_func = custom_func,
    settings = list(
      iiv = FALSE,
      error = FALSE,
      plot = FALSE,
      N = 2,  # Small number of ants
      max_models = 4,
      save_mode = "none",
      xgboost = FALSE  # Disable for simpler testing
    )
  )
  
  expect_type(result, "list")
  expect_true("results" %in% names(result))
  expect_true("best_model" %in% names(result))
  expect_s3_class(result$results, "data.frame")
  expect_true(nrow(result$results) > 0)
  expect_true("metric" %in% names(result$results))
})

test_that("ACO works with custom models and IIV", {
  skip_if_no_monolix()
  
  custom_func <- create_simple_custom_func()
  param_mapping <- create_simple_param_mapping()
  
  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = c("firstOrder", "zeroOrder")),
    algorithm = "ACO",
    model_creation_func = custom_func,
    param_mapping = param_mapping,
    settings = list(
      iiv = TRUE,
      error = FALSE,
      plot = FALSE,
      N = 2,
      max_models = 6,
      save_mode = "none",
      xgboost = FALSE,
      initial_iiv_forced_iter = 1  # Reduce forced iterations for testing
    )
  )
  
  expect_type(result, "list")
  expect_true("results" %in% names(result))
  expect_s3_class(result$results, "data.frame")
  
  # Check that IIV parameters are in the results
  param_cols <- intersect(c("V", "Cl", "ka", "Tk0"), names(result$results))
  expect_true(length(param_cols) > 0)
})

test_that("ACO works with custom models and error models", {
  skip_if_no_monolix()
  
  custom_func <- create_simple_custom_func()
  
  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "ACO",
    model_creation_func = custom_func,
    settings = list(
      iiv = FALSE,
      error = TRUE,
      plot = FALSE,
      N = 2,
      max_models = 4,
      save_mode = "none",
      xgboost = FALSE
    )
  )
  
  expect_type(result, "list")
  expect_true("results" %in% names(result))
  expect_s3_class(result$results, "data.frame")
  expect_true("error" %in% names(result$results))
})

test_that("ACO with custom models and both IIV and error", {
  skip_if_no_monolix()
  
  custom_func <- create_simple_custom_func()
  param_mapping <- create_simple_param_mapping()
  
  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "ACO",
    model_creation_func = custom_func,
    param_mapping = param_mapping,
    settings = list(
      iiv = TRUE,
      error = TRUE,
      plot = FALSE,
      N = 2,
      max_models = 6,
      save_mode = "none",
      xgboost = FALSE,
      initial_iiv_forced_iter = 1
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

test_that("ACO with filter-based custom models", {
  skip_if_no_monolix()
  
  custom_func <- create_filter_based_custom_func()
  
  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = c("firstOrder", "zeroOrder")),
    algorithm = "ACO",
    model_creation_func = custom_func,
    settings = list(
      iiv = FALSE,
      error = FALSE,
      plot = FALSE,
      N = 2,
      max_models = 4,
      save_mode = "none",
      xgboost = FALSE
    )
  )
  
  expect_type(result, "list")
  expect_true("results" %in% names(result))
  expect_s3_class(result$results, "data.frame")
  expect_true("absorption" %in% names(result$results))
})

test_that("ACO respects max_models limit with custom models", {
  skip_if_no_monolix()
  
  custom_func <- create_simple_custom_func()
  
  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = c("firstOrder", "zeroOrder")),
    algorithm = "ACO",
    model_creation_func = custom_func,
    settings = list(
      iiv = FALSE,
      error = FALSE,
      plot = FALSE,
      N = 2,
      max_models = 3,  # Limit to 3 models
      save_mode = "none",
      xgboost = FALSE
    )
  )
  
  expect_type(result, "list")
  expect_true(nrow(result$results) <= 3)
})

test_that("ACO with custom models returns best model", {
  skip_if_no_monolix()
  
  custom_func <- create_simple_custom_func()
  
  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = c("firstOrder", "zeroOrder")),
    algorithm = "ACO",
    model_creation_func = custom_func,
    settings = list(
      iiv = FALSE,
      error = FALSE,
      plot = FALSE,
      N = 2,
      max_models = 4,
      save_mode = "none",
      xgboost = FALSE
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

test_that("ACO with custom models returns iteration progress", {
  skip_if_no_monolix()
  
  custom_func <- create_simple_custom_func()
  
  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = c("firstOrder", "zeroOrder")),
    algorithm = "ACO",
    model_creation_func = custom_func,
    settings = list(
      iiv = FALSE,
      error = FALSE,
      plot = FALSE,
      N = 2,
      max_models = 4,
      save_mode = "none",
      xgboost = FALSE
    )
  )
  
  # ACO should return iteration_progress
  expect_true("iteration_progress" %in% names(result))
})

test_that("ACO with custom models and param_mapping filters IIV correctly", {
  skip_if_no_monolix()
  
  custom_func <- create_filter_based_custom_func()
  param_mapping <- create_simple_param_mapping()
  
  result <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = c("firstOrder", "zeroOrder")),
    algorithm = "ACO",
    model_creation_func = custom_func,
    param_mapping = param_mapping,
    settings = list(
      iiv = TRUE,
      error = FALSE,
      plot = FALSE,
      N = 2,
      max_models = 6,
      save_mode = "none",
      xgboost = FALSE,
      initial_iiv_forced_iter = 1
    )
  )
  
  expect_type(result, "list")
  expect_s3_class(result$results, "data.frame")
  
  # Check that results have appropriate IIV columns
  # For firstOrder: should have V, Cl, ka
  # For zeroOrder: should have V, Cl, Tk0
  fo_rows <- result$results[result$results$absorption == "firstOrder", ]
  zo_rows <- result$results[result$results$absorption == "zeroOrder", ]
  
  if (nrow(fo_rows) > 0) {
    # firstOrder models should have ka column
    expect_true("ka" %in% names(fo_rows))
  }
  
  if (nrow(zo_rows) > 0) {
    # zeroOrder models should have Tk0 column
    expect_true("Tk0" %in% names(zo_rows))
  }
})

test_that("ACO with custom models handles elitism setting", {
  skip_if_no_monolix()
  
  custom_func <- create_simple_custom_func()
  
  # Test with elitism enabled
  result_elitism <- findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "ACO",
    model_creation_func = custom_func,
    settings = list(
      iiv = FALSE,
      error = FALSE,
      plot = FALSE,
      N = 2,
      max_models = 4,
      save_mode = "none",
      xgboost = FALSE,
      elitism = TRUE
    )
  )
  
  expect_type(result_elitism, "list")
  expect_true(nrow(result_elitism$results) > 0)
})
