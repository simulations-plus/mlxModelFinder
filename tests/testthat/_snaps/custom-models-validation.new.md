# library='custom' requires model_creation_func

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search")
    Message
      i Search space size: 1 models
    Condition
      Error:
      ! Error: 'model_creation_func' must be provided as a function when library='custom'.

# library='custom' rejects decision_tree algorithm

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "decision_tree", model_creation_func = valid_custom_func)
    Message
      i Search space size: 1 models
    Condition
      Error:
      ! Error: library='custom' is only supported with algorithm='exhaustive_search' or 'ACO'.

# library='custom' with iiv=TRUE requires param_mapping

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = valid_custom_func,
      settings = list(iiv = TRUE, plot = FALSE, max_models = 1))
    Message
      i Search space size: 8 models
    Condition
      Error:
      ! Error: 'param_mapping' must be provided as a list when library='custom' and settings$iiv=TRUE.

# param_mapping must have correct structure

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = valid_custom_func,
      param_mapping = list(common_params = c("V", "Cl")), settings = list(iiv = TRUE,
        plot = FALSE, max_models = 1))
    Message
      i Search space size: 4 models
    Condition
      Error:
      ! Error: 'param_mapping' must contain 'filter_params' and 'common_params' elements.

---

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = valid_custom_func,
      param_mapping = list(filter_params = list()), settings = list(iiv = TRUE, plot = FALSE,
        max_models = 1))
    Message
      i Search space size: 1 models
    Condition
      Error:
      ! Error: 'param_mapping' must contain 'filter_params' and 'common_params' elements.

---

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = valid_custom_func,
      param_mapping = list(common_params = list("V", "Cl"), filter_params = list()),
      settings = list(iiv = TRUE, plot = FALSE, max_models = 1))
    Message
      i Search space size: 4 models
    Condition
      Error:
      ! Error: 'param_mapping$common_params' must be a character vector.

---

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = valid_custom_func,
      param_mapping = list(common_params = c("V", "Cl"), filter_params = c("ka",
        "Tk0")), settings = list(iiv = TRUE, plot = FALSE, max_models = 1))
    Message
      i Search space size: 4 models
    Condition
      Error:
      ! Error: 'param_mapping$filter_params' must be a list.

# model_creation_func must have correct signature

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = bad_func_1arg)
    Message
      i Search space size: 1 models
    Condition
      Error:
      ! Error: 'model_creation_func' must accept exactly 2 arguments: library and filters.

---

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = bad_func_3args)
    Message
      i Search space size: 1 models
    Condition
      Error:
      ! Error: 'model_creation_func' must accept exactly 2 arguments: library and filters.

# model_creation_func must return valid path

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = func_returns_null,
      settings = tiny_settings())
    Message
      i Search space size: 1 models
      
      -- Generating full grid search job list... -------------------------------------
    Condition
      Error:
      ! Error: custom model_creation_func must return a single character string (model path).

---

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = func_returns_bad_path,
      settings = tiny_settings())
    Message
      i Search space size: 1 models
      
      -- Generating full grid search job list... -------------------------------------
    Condition
      Error:
      ! Error: custom model_creation_func returned non-existent path: path/that/does/not/exist.txt

---

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = func_returns_multiple,
      settings = tiny_settings())
    Message
      i Search space size: 1 models
      
      -- Generating full grid search job list... -------------------------------------
    Condition
      Error:
      ! Error: custom model_creation_func must return a single character string (model path).

