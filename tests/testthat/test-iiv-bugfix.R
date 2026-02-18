# tests/testthat/test-iiv-bugfix.R
# Tests for the IIV bugfix in exhaustive_search with custom models
# Bug: as.data.frame(vector) instead of as.data.frame(as.list(vector))
# caused IIV combinations to not be generated for custom models.

# demo_proj is defined in helper-monolix.R

test_that("exhaustive search generates IIV combinations for custom models", {
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
      save_mode = "none"
    )
  )

  # With 3 IIV parameters (V, Cl, ka), should get 2^3 = 8 combinations
  expect_equal(nrow(result$results), 8)

  # IIV columns should have both TRUE and FALSE values
  iiv_cols <- intersect(c("V", "Cl", "ka"), names(result$results))
  expect_length(iiv_cols, 3)

  for (col in iiv_cols) {
    unique_vals <- unique(na.omit(result$results[[col]]))
    expect_true(TRUE %in% unique_vals, info = paste(col, "should have TRUE"))
    expect_true(FALSE %in% unique_vals, info = paste(col, "should have FALSE"))
  }
})

test_that("exhaustive search generates IIV combinations for filter-based custom models", {
  skip_if_no_monolix()

  custom_func <- create_filter_based_custom_func()
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
      save_mode = "none"
    )
  )

  # firstOrder: V, Cl, ka -> 8 combos
  # zeroOrder: V, Cl, Tk0 -> 8 combos
  # Total: 16
  expect_equal(nrow(result$results), 16)

  # Both absorption types should be present
  expect_true(any(result$results$absorption == "firstOrder"))
  expect_true(any(result$results$absorption == "zeroOrder"))

  # IIV columns should vary
  expect_true(length(unique(na.omit(result$results$V))) == 2)
  expect_true(length(unique(na.omit(result$results$Cl))) == 2)
})

test_that("IIV combos correct count with single filter custom model", {
  skip_if_no_monolix()

  custom_func <- create_simple_custom_func()
  param_mapping <- list(
    common_params = c("V", "Cl", "ka"),
    filter_params = list(
      absorption = list(
        firstOrder = character(0)
      )
    )
  )

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
      save_mode = "none"
    )
  )

  # 1 structural model * 2^3 IIV combos = 8
  expect_equal(nrow(result$results), 8)
})
