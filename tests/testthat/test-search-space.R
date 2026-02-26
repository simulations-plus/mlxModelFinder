# tests/testthat/test-search-space.R

test_that("Search space calculation for custom models without IIV", {
  # Simple case: 2 absorption types
  filters <- list(absorption = c("firstOrder", "zeroOrder"))
  space <- mlxModelFinder:::calculateSearchSpace("custom", filters, error = FALSE, iiv = FALSE)
  expect_equal(space, 2)
  
  # Multi-dimensional: 2 x 2
  filters <- list(
    absorption = c("firstOrder", "zeroOrder"),
    distribution = c("1cpt", "2cpt")
  )
  space <- mlxModelFinder:::calculateSearchSpace("custom", filters, error = FALSE, iiv = FALSE)
  expect_equal(space, 4)
  
  # Three dimensions: 2 x 2 x 3
  filters <- list(
    absorption = c("firstOrder", "zeroOrder"),
    distribution = c("1cpt", "2cpt"),
    elimination = c("linear", "MM", "combined")
  )
  space <- mlxModelFinder:::calculateSearchSpace("custom", filters, error = FALSE, iiv = FALSE)
  expect_equal(space, 12)
})

test_that("Search space calculation with error models", {
  filters <- list(absorption = c("firstOrder", "zeroOrder"))
  
  # Without error
  space_no_error <- mlxModelFinder:::calculateSearchSpace("custom", filters, error = FALSE, iiv = FALSE)
  
  # With error (4 error models)
  space_with_error <- mlxModelFinder:::calculateSearchSpace("custom", filters, error = TRUE, iiv = FALSE)
  
  expect_equal(space_with_error, space_no_error * 4)
  expect_equal(space_with_error, 8)  # 2 structural * 4 error

  # With error and 2 observations (4^2 = 16 error combinations)
  space_2obs <- mlxModelFinder:::calculateSearchSpace("custom", filters, error = TRUE, iiv = FALSE, n_obs = 2)
  expect_equal(space_2obs, space_no_error * 16)
  expect_equal(space_2obs, 32)  # 2 structural * 16 error

  # With error and 3 observations (4^3 = 64 error combinations)
  space_3obs <- mlxModelFinder:::calculateSearchSpace("custom", filters, error = TRUE, iiv = FALSE, n_obs = 3)
  expect_equal(space_3obs, space_no_error * 64)
  expect_equal(space_3obs, 128)  # 2 structural * 64 error
})

test_that("Search space calculation with IIV and param_mapping", {
  # Simple param_mapping: all models have same 3 parameters
  param_mapping <- list(
    common_params = c("V", "Cl"),
    filter_params = list(
      absorption = list(
        firstOrder = c("ka"),
        zeroOrder = c("Tk0")
      )
    )
  )
  
  filters <- list(absorption = c("firstOrder", "zeroOrder"))
  
  # Each model has 3 params (V, Cl, ka/Tk0)
  # Each param can have IIV on/off: 2^3 = 8 combinations
  # 2 structural models * 8 IIV combinations = 16 total
  space <- mlxModelFinder:::calculateSearchSpace("custom", filters, error = FALSE, iiv = TRUE, param_mapping)
  expect_equal(space, 16)
})

test_that("Search space calculation with varying parameter counts", {
  # Different models have different number of parameters
  param_mapping <- list(
    common_params = c("V", "Cl"),
    filter_params = list(
      absorption = list(
        firstOrder = c("ka"),
        zeroOrder = c("Tk0")
      ),
      distribution = list(
        "1cpt" = character(0),      # No additional params
        "2cpt" = c("Q", "V2")        # 2 additional params
      )
    )
  )
  
  filters <- list(
    absorption = c("firstOrder", "zeroOrder"),
    distribution = c("1cpt", "2cpt")
  )
  
  # Calculate expected space:
  # firstOrder + 1cpt: V, Cl, ka (3 params) -> 2^3 = 8
  # firstOrder + 2cpt: V, Cl, ka, Q, V2 (5 params) -> 2^5 = 32
  # zeroOrder + 1cpt: V, Cl, Tk0 (3 params) -> 2^3 = 8
  # zeroOrder + 2cpt: V, Cl, Tk0, Q, V2 (5 params) -> 2^5 = 32
  # Total: 8 + 32 + 8 + 32 = 80
  
  space <- mlxModelFinder:::calculateSearchSpace("custom", filters, error = FALSE, iiv = TRUE, param_mapping)
  expect_equal(space, 80)
})

test_that("Search space calculation with IIV and error models", {
  param_mapping <- list(
    common_params = c("V", "Cl"),
    filter_params = list(
      absorption = list(
        firstOrder = c("ka"),
        zeroOrder = c("Tk0")
      )
    )
  )
  
  filters <- list(absorption = c("firstOrder", "zeroOrder"))
  
  # Base space with IIV: 16 (2 structural * 8 IIV each)
  # With 4 error models: 16 * 4 = 64
  space <- mlxModelFinder:::calculateSearchSpace("custom", filters, error = TRUE, iiv = TRUE, param_mapping)
  expect_equal(space, 64)
})

test_that("Search space calculation without param_mapping uses conservative estimate", {
  filters <- list(absorption = c("firstOrder", "zeroOrder"))
  
  # Without param_mapping, should use conservative estimate
  # 2 structural * 8 (assumed 2^3) = 16
  space <- mlxModelFinder:::calculateSearchSpace("custom", filters, error = FALSE, iiv = TRUE, param_mapping = NULL)
  expect_equal(space, 16)
})

test_that("Search space calculation with only common params", {
  # All models have only common parameters
  param_mapping <- list(
    common_params = c("V", "Cl", "ka"),
    filter_params = list()
  )
  
  filters <- list(absorption = c("firstOrder", "zeroOrder"))
  
  # Both models have same 3 params: 2 structural * 2^3 = 16
  space <- mlxModelFinder:::calculateSearchSpace("custom", filters, error = FALSE, iiv = TRUE, param_mapping)
  expect_equal(space, 16)
})

test_that("Search space calculation with no common params", {
  # Only filter-specific parameters
  param_mapping <- list(
    common_params = character(0),
    filter_params = list(
      absorption = list(
        firstOrder = c("V", "Cl", "ka"),
        zeroOrder = c("V", "Cl", "Tk0")
      )
    )
  )
  
  filters <- list(absorption = c("firstOrder", "zeroOrder"))
  
  # Each model has 3 params: 2 structural * 2^3 = 16
  space <- mlxModelFinder:::calculateSearchSpace("custom", filters, error = FALSE, iiv = TRUE, param_mapping)
  expect_equal(space, 16)
})

test_that("Search space calculation with single filter value", {
  param_mapping <- list(
    common_params = c("V", "Cl"),
    filter_params = list(
      absorption = list(
        firstOrder = c("ka")
      )
    )
  )
  
  filters <- list(absorption = "firstOrder")
  
  # 1 structural model with 3 params: 2^3 = 8
  space <- mlxModelFinder:::calculateSearchSpace("custom", filters, error = FALSE, iiv = TRUE, param_mapping)
  expect_equal(space, 8)
})

test_that("Search space calculation with empty filters", {
  filters <- list()
  
  # No filters means 1 model
  space <- mlxModelFinder:::calculateSearchSpace("custom", filters, error = FALSE, iiv = FALSE)
  expect_equal(space, 1)
})

test_that("Search space calculation matches structural space when no IIV", {
  filters <- list(
    absorption = c("firstOrder", "zeroOrder"),
    distribution = c("1cpt", "2cpt")
  )
  
  param_mapping <- list(
    common_params = c("V", "Cl"),
    filter_params = list(
      absorption = list(
        firstOrder = c("ka"),
        zeroOrder = c("Tk0")
      )
    )
  )
  
  # Without IIV, param_mapping shouldn't matter
  space_no_mapping <- mlxModelFinder:::calculateSearchSpace("custom", filters, error = FALSE, iiv = FALSE, NULL)
  space_with_mapping <- mlxModelFinder:::calculateSearchSpace("custom", filters, error = FALSE, iiv = FALSE, param_mapping)
  
  expect_equal(space_no_mapping, 4)
  expect_equal(space_with_mapping, 4)
})

test_that("Search space calculation with complex multi-dimensional param_mapping", {
  param_mapping <- list(
    common_params = c("V"),
    filter_params = list(
      absorption = list(
        firstOrder = c("ka", "Cl"),
        zeroOrder = c("Tk0", "Cl"),
        sigmoid = c("ka", "Tk0s", "Cl")
      ),
      distribution = list(
        "1cpt" = character(0),
        "2cpt" = c("Q", "V2"),
        "3cpt" = c("Q", "V2", "Q2", "V3")
      )
    )
  )
  
  filters <- list(
    absorption = c("firstOrder", "zeroOrder"),
    distribution = c("1cpt", "2cpt")
  )
  
  # Calculate expected:
  # firstOrder + 1cpt: V, ka, Cl (3 params) -> 2^3 = 8
  # firstOrder + 2cpt: V, ka, Cl, Q, V2 (5 params) -> 2^5 = 32
  # zeroOrder + 1cpt: V, Tk0, Cl (3 params) -> 2^3 = 8
  # zeroOrder + 2cpt: V, Tk0, Cl, Q, V2 (5 params) -> 2^5 = 32
  # Total: 80
  
  space <- mlxModelFinder:::calculateSearchSpace("custom", filters, error = FALSE, iiv = TRUE, param_mapping)
  expect_equal(space, 80)
})

test_that("Search space calculation handles character(0) in filter_params", {
  param_mapping <- list(
    common_params = c("V", "Cl"),
    filter_params = list(
      absorption = list(
        firstOrder = c("ka"),
        zeroOrder = character(0)  # No additional params for zeroOrder
      )
    )
  )
  
  filters <- list(absorption = c("firstOrder", "zeroOrder"))
  
  # firstOrder: V, Cl, ka (3 params) -> 2^3 = 8
  # zeroOrder: V, Cl (2 params) -> 2^2 = 4
  # Total: 12
  
  space <- mlxModelFinder:::calculateSearchSpace("custom", filters, error = FALSE, iiv = TRUE, param_mapping)
  expect_equal(space, 12)
})
