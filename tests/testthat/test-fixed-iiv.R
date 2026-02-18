# tests/testthat/test-fixed-iiv.R
# Tests for the fixed_iiv feature

# demo_proj is defined in helper-monolix.R

# --- Validation tests (no Monolix needed) ---

test_that("validate_fixed_iiv rejects non-list input", {
  expect_error(
    validate_fixed_iiv(TRUE),
    "must be a named list"
  )
  expect_error(
    validate_fixed_iiv(c(Cl = TRUE, V = FALSE)),
    "must be a named list"
  )
})

test_that("validate_fixed_iiv rejects unnamed list", {
  expect_error(
    validate_fixed_iiv(list(TRUE, FALSE)),
    "must be a named list"
  )
})

test_that("validate_fixed_iiv rejects non-logical values", {
  expect_error(
    validate_fixed_iiv(list(Cl = "yes")),
    "must be TRUE or FALSE"
  )
  expect_error(
    validate_fixed_iiv(list(Cl = TRUE, V = 1)),
    "must be TRUE or FALSE"
  )
})

test_that("validate_fixed_iiv accepts valid input", {
  expect_true(validate_fixed_iiv(list(Cl = TRUE)))
  expect_true(validate_fixed_iiv(list(Cl = TRUE, V = FALSE)))
  expect_true(validate_fixed_iiv(list(ka = FALSE)))
})

# --- apply_fixed_iiv unit tests (no Monolix needed) ---

test_that("apply_fixed_iiv with NULL fixed_iiv generates all combos", {
  default_iiv <- data.frame(V = TRUE, Cl = TRUE, ka = TRUE)
  result <- apply_fixed_iiv(default_iiv, NULL)
  expect_equal(nrow(result), 8)  # 2^3
  expect_true(all(c("V", "Cl", "ka") %in% names(result)))
})

test_that("apply_fixed_iiv locks TRUE params", {
  default_iiv <- data.frame(V = TRUE, Cl = TRUE, ka = TRUE)
  result <- apply_fixed_iiv(default_iiv, list(V = TRUE))
  # 2 free params (Cl, ka) -> 2^2 = 4 combos, V always TRUE
  expect_equal(nrow(result), 4)
  expect_true(all(result$V == TRUE))
  expect_equal(length(unique(result$Cl)), 2)
  expect_equal(length(unique(result$ka)), 2)
})

test_that("apply_fixed_iiv locks FALSE params", {
  default_iiv <- data.frame(V = TRUE, Cl = TRUE, ka = TRUE)
  result <- apply_fixed_iiv(default_iiv, list(ka = FALSE))
  # 2 free params (V, Cl) -> 2^2 = 4 combos, ka always FALSE
  expect_equal(nrow(result), 4)
  expect_true(all(result$ka == FALSE))
  expect_equal(length(unique(result$V)), 2)
  expect_equal(length(unique(result$Cl)), 2)
})

test_that("apply_fixed_iiv locks multiple params", {
  default_iiv <- data.frame(V = TRUE, Cl = TRUE, ka = TRUE)
  result <- apply_fixed_iiv(default_iiv, list(V = TRUE, ka = FALSE))
  # 1 free param (Cl) -> 2^1 = 2 combos
  expect_equal(nrow(result), 2)
  expect_true(all(result$V == TRUE))
  expect_true(all(result$ka == FALSE))
  expect_equal(length(unique(result$Cl)), 2)
})

test_that("apply_fixed_iiv with all params fixed gives 1 combo", {
  default_iiv <- data.frame(V = TRUE, Cl = TRUE)
  result <- apply_fixed_iiv(default_iiv, list(V = TRUE, Cl = FALSE))
  expect_equal(nrow(result), 1)
  expect_true(result$V)
  expect_false(result$Cl)
})

test_that("apply_fixed_iiv ignores params not in default_iiv", {
  default_iiv <- data.frame(V = TRUE, Cl = TRUE)
  result <- apply_fixed_iiv(default_iiv, list(V = TRUE, Tk0 = FALSE))
  # Only V is locked, Cl is free -> 2^1 = 2
  expect_equal(nrow(result), 2)
  expect_true(all(result$V == TRUE))
  expect_false("Tk0" %in% names(result))
})

# --- Integration tests (require Monolix) ---

test_that("fixed_iiv locks params in exhaustive search with custom models", {
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
    fixed_iiv = list(V = TRUE, ka = FALSE),
    settings = list(
      iiv = TRUE,
      error = FALSE,
      plot = FALSE,
      save_mode = "none"
    )
  )

  # 1 free param (Cl) -> 2^1 = 2 combos
  expect_equal(nrow(result$results), 2)
  # V always TRUE, ka always FALSE
  expect_true(all(result$results$V == TRUE))
  expect_true(all(result$results$ka == FALSE))
  # Cl varies
  expect_equal(length(unique(result$results$Cl)), 2)
})

test_that("fixed_iiv works with pk library exhaustive search", {
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
    algorithm = "exhaustive_search",
    fixed_iiv = list(Cl = TRUE, V = TRUE),
    settings = list(
      iiv = TRUE,
      error = FALSE,
      plot = FALSE,
      save_mode = "none"
    )
  )

  expect_type(result, "list")
  expect_true(nrow(result$results) > 0)
  # Cl and V should always be TRUE
  expect_true(all(result$results$Cl == TRUE))
  expect_true(all(result$results$V == TRUE))
})

test_that("fixed_iiv works with decision tree", {
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
    fixed_iiv = list(Cl = TRUE),
    settings = list(
      iiv = TRUE,
      error = FALSE,
      plot = FALSE,
      save_mode = "none"
    )
  )

  expect_type(result, "list")
  expect_true(nrow(result$results) > 0)
})

test_that("findModel works with fixed_iiv = NULL (default)", {
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
    settings = list(
      iiv = TRUE,
      error = FALSE,
      plot = FALSE,
      save_mode = "none"
    )
  )

  # Without fixed_iiv, all 2^3 = 8 combos
  expect_equal(nrow(result$results), 8)
})
