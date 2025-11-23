# tests/testthat/test-custom-models-validation.R

# demo_proj is defined in helper-monolix.R

# Helper: create a simple valid custom model function
valid_custom_func <- function(library, filters) {
  # Return path to a simple 1-compartment model
  system.file("extdata", "theophylline_project.mlxtran", package = "mlxModelFinder")
}

# Helper: create a valid param_mapping
valid_param_mapping <- function() {
  list(
    common_params = c("V", "Cl"),
    filter_params = list(
      absorption = list(
        firstOrder = c("ka"),
        zeroOrder = c("Tk0")
      )
    )
  )
}

test_that("library='custom' requires model_creation_func", {
  skip_if_no_monolix()
  
  expect_snapshot(error = TRUE, findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search"
  ))
})

test_that("library='custom' rejects decision_tree algorithm", {
  skip_if_no_monolix()
  
  expect_snapshot(error = TRUE, findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "decision_tree",
    model_creation_func = valid_custom_func
  ))
})

test_that("library='custom' with iiv=TRUE requires param_mapping", {
  skip_if_no_monolix()
  
  expect_snapshot(error = TRUE, findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = valid_custom_func,
    settings = list(iiv = TRUE, plot = FALSE, max_models = 1)
  ))
})

test_that("param_mapping must have correct structure", {
  skip_if_no_monolix()
  
  # Missing filter_params
  expect_snapshot(error = TRUE, findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = valid_custom_func,
    param_mapping = list(common_params = c("V", "Cl")),
    settings = list(iiv = TRUE, plot = FALSE, max_models = 1)
  ))
  
  # Missing common_params
  expect_snapshot(error = TRUE, findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = valid_custom_func,
    param_mapping = list(filter_params = list()),
    settings = list(iiv = TRUE, plot = FALSE, max_models = 1)
  ))
  
  # common_params not a character vector
  expect_snapshot(error = TRUE, findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = valid_custom_func,
    param_mapping = list(
      common_params = list("V", "Cl"),  # Should be character vector
      filter_params = list()
    ),
    settings = list(iiv = TRUE, plot = FALSE, max_models = 1)
  ))
  
  # filter_params not a list
  expect_snapshot(error = TRUE, findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = valid_custom_func,
    param_mapping = list(
      common_params = c("V", "Cl"),
      filter_params = c("ka", "Tk0")  # Should be a list
    ),
    settings = list(iiv = TRUE, plot = FALSE, max_models = 1)
  ))
})

test_that("model_creation_func must have correct signature", {
  skip_if_no_monolix()
  
  # Function with wrong number of arguments
  bad_func_1arg <- function(library) {
    "path/to/model.txt"
  }
  
  expect_snapshot(error = TRUE, findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = bad_func_1arg
  ))
  
  # Function with too many arguments
  bad_func_3args <- function(library, filters, extra) {
    "path/to/model.txt"
  }
  
  expect_snapshot(error = TRUE, findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = bad_func_3args
  ))
})

test_that("model_creation_func must return valid path", {
  skip_if_no_monolix()
  
  # Function returns NULL
  func_returns_null <- function(library, filters) {
    NULL
  }
  
  expect_snapshot(error = TRUE, findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = func_returns_null,
    settings = tiny_settings()
  ))
  
  # Function returns non-existent path
  func_returns_bad_path <- function(library, filters) {
    "path/that/does/not/exist.txt"
  }
  
  expect_snapshot(error = TRUE, findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = func_returns_bad_path,
    settings = tiny_settings()
  ))
  
  # Function returns multiple paths
  func_returns_multiple <- function(library, filters) {
    c("path1.txt", "path2.txt")
  }
  
  expect_snapshot(error = TRUE, findModel(
    project = demo_proj,
    library = "custom",
    filters = list(absorption = "firstOrder"),
    algorithm = "exhaustive_search",
    model_creation_func = func_returns_multiple,
    settings = tiny_settings()
  ))
})

test_that("model_creation_func is ignored when library != 'custom'", {
  skip_if_no_monolix()
  
  # Should work and just ignore the custom function
  expect_warning(
    findModel(
      project = demo_proj,
      library = "pk",
      requiredFilters = list(administration = "oral", parametrization = "clearance"),
      filters = tiny_pk_filters(),
      model_creation_func = valid_custom_func,
      settings = tiny_settings()
    ),
    "model_creation_func.*ignored"
  )
})

test_that("custom library accepts exhaustive_search and ACO algorithms", {
  skip_if_no_monolix()
  
  # These should not error during validation (may fail later for other reasons)
  expect_no_error(
    tryCatch(
      findModel(
        project = demo_proj,
        library = "custom",
        filters = list(absorption = "firstOrder"),
        algorithm = "exhaustive_search",
        model_creation_func = valid_custom_func,
        settings = tiny_settings()
      ),
      error = function(e) {
        # Only re-throw if it's a validation error about algorithm
        if (grepl("algorithm.*not supported", e$message, ignore.case = TRUE)) {
          stop(e)
        }
      }
    )
  )
  
  expect_no_error(
    tryCatch(
      findModel(
        project = demo_proj,
        library = "custom",
        filters = list(absorption = "firstOrder"),
        algorithm = "ACO",
        model_creation_func = valid_custom_func,
        settings = tiny_settings()
      ),
      error = function(e) {
        # Only re-throw if it's a validation error about algorithm
        if (grepl("algorithm.*not supported", e$message, ignore.case = TRUE)) {
          stop(e)
        }
      }
    )
  )
})
