# tests/testthat/test-param-mapping.R

test_that("getParamListFromMapping extracts common params", {
  param_mapping <- list(
    common_params = c("V", "Cl"),
    filter_params = list()
  )
  
  filters <- list(absorption = "firstOrder")
  
  result <- mlxModelFinder:::getParamListFromMapping(filters, param_mapping)
  
  expect_equal(sort(result), sort(c("V", "Cl")))
})

test_that("getParamListFromMapping extracts filter-specific params", {
  param_mapping <- list(
    common_params = c("V", "Cl"),
    filter_params = list(
      absorption = list(
        firstOrder = c("ka"),
        zeroOrder = c("Tk0")
      )
    )
  )
  
  filters <- list(absorption = "firstOrder")
  
  result <- mlxModelFinder:::getParamListFromMapping(filters, param_mapping)
  
  expect_equal(sort(result), sort(c("V", "Cl", "ka")))
})

test_that("getParamListFromMapping handles multiple filter categories", {
  param_mapping <- list(
    common_params = c("V", "Cl"),
    filter_params = list(
      absorption = list(
        firstOrder = c("ka"),
        zeroOrder = c("Tk0")
      ),
      distribution = list(
        "1compartment" = character(0),
        "2compartments" = c("Q", "V2")
      )
    )
  )
  
  filters <- list(
    absorption = "firstOrder",
    distribution = "2compartments"
  )
  
  result <- mlxModelFinder:::getParamListFromMapping(filters, param_mapping)
  
  expect_equal(sort(result), sort(c("V", "Cl", "ka", "Q", "V2")))
})

test_that("getParamListFromMapping returns unique parameters", {
  param_mapping <- list(
    common_params = c("V", "Cl", "V"),  # Duplicate V
    filter_params = list(
      absorption = list(
        firstOrder = c("ka", "V")  # V also in common_params
      )
    )
  )
  
  filters <- list(absorption = "firstOrder")
  
  result <- mlxModelFinder:::getParamListFromMapping(filters, param_mapping)
  
  # Should have unique values only
  expect_equal(length(result), length(unique(result)))
  expect_true("V" %in% result)
  expect_true("Cl" %in% result)
  expect_true("ka" %in% result)
})

test_that("getParamListFromMapping handles missing filter values gracefully", {
  param_mapping <- list(
    common_params = c("V", "Cl"),
    filter_params = list(
      absorption = list(
        firstOrder = c("ka")
      )
    )
  )
  
  # Filter value not in mapping
  filters <- list(absorption = "sigmoid")
  
  result <- mlxModelFinder:::getParamListFromMapping(filters, param_mapping)
  
  # Should still return common params
  expect_equal(sort(result), sort(c("V", "Cl")))
})

test_that("getParamListFromMapping handles missing filter categories gracefully", {
  param_mapping <- list(
    common_params = c("V", "Cl"),
    filter_params = list(
      absorption = list(
        firstOrder = c("ka")
      )
    )
  )
  
  # Filter category not in mapping
  filters <- list(
    absorption = "firstOrder",
    elimination = "linear"  # Not in param_mapping
  )
  
  result <- mlxModelFinder:::getParamListFromMapping(filters, param_mapping)
  
  # Should return common params + absorption params
  expect_equal(sort(result), sort(c("V", "Cl", "ka")))
})

test_that("getParamListFromMapping returns empty vector when param_mapping is NULL", {
  result <- mlxModelFinder:::getParamListFromMapping(list(absorption = "firstOrder"), NULL)
  
  expect_equal(result, character(0))
})

test_that("getParamListFromMapping handles empty common_params", {
  param_mapping <- list(
    common_params = NULL,
    filter_params = list(
      absorption = list(
        firstOrder = c("ka")
      )
    )
  )
  
  filters <- list(absorption = "firstOrder")
  
  result <- mlxModelFinder:::getParamListFromMapping(filters, param_mapping)
  
  expect_equal(result, "ka")
})

test_that("getParamListFromMapping handles empty filter_params", {
  param_mapping <- list(
    common_params = c("V", "Cl"),
    filter_params = NULL
  )
  
  filters <- list(absorption = "firstOrder")
  
  result <- mlxModelFinder:::getParamListFromMapping(filters, param_mapping)
  
  expect_equal(sort(result), sort(c("V", "Cl")))
})

test_that("getParamListFromMapping handles complex nested structure", {
  param_mapping <- list(
    common_params = c("V", "Cl"),
    filter_params = list(
      absorption = list(
        firstOrder = c("ka"),
        zeroOrder = c("Tk0"),
        transitCompartments = c("Ktr", "Mtt")
      ),
      distribution = list(
        "1compartment" = character(0),
        "2compartments" = c("Q", "V2"),
        "3compartments" = c("Q", "V2", "Q2", "V3")
      ),
      elimination = list(
        linear = character(0),
        MichaelisMenten = c("Vm", "Km")
      )
    )
  )
  
  filters <- list(
    absorption = "transitCompartments",
    distribution = "3compartments",
    elimination = "MichaelisMenten"
  )
  
  result <- mlxModelFinder:::getParamListFromMapping(filters, param_mapping)
  
  expected <- c("V", "Cl", "Ktr", "Mtt", "Q", "V2", "Q2", "V3", "Vm", "Km")
  expect_equal(sort(result), sort(unique(expected)))
})

test_that("getParamListFromMapping handles character(0) filter values", {
  param_mapping <- list(
    common_params = c("V", "Cl"),
    filter_params = list(
      distribution = list(
        "1compartment" = character(0),  # No additional params
        "2compartments" = c("Q", "V2")
      )
    )
  )
  
  filters <- list(distribution = "1compartment")
  
  result <- mlxModelFinder:::getParamListFromMapping(filters, param_mapping)
  
  # Should only have common params
  expect_equal(sort(result), sort(c("V", "Cl")))
})
