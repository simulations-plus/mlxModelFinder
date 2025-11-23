# library='custom' requires model_creation_func

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search")
    Condition
      Error in `calculateSearchSpace()`:
      ! object 'space' not found

# library='custom' rejects decision_tree algorithm

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "decision_tree", model_creation_func = valid_custom_func)
    Condition
      Error in `calculateSearchSpace()`:
      ! object 'space' not found

# library='custom' with iiv=TRUE requires param_mapping

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = valid_custom_func,
      settings = list(iiv = TRUE, plot = FALSE, max_models = 1))
    Condition
      Error in `calculateSearchSpace()`:
      ! object 'nb_params' not found

# param_mapping must have correct structure

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = valid_custom_func,
      param_mapping = list(common_params = c("V", "Cl")), settings = list(iiv = TRUE,
        plot = FALSE, max_models = 1))
    Condition
      Error in `calculateSearchSpace()`:
      ! object 'nb_params' not found

---

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = valid_custom_func,
      param_mapping = list(filter_params = list()), settings = list(iiv = TRUE, plot = FALSE,
        max_models = 1))
    Condition
      Error in `calculateSearchSpace()`:
      ! object 'nb_params' not found

---

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = valid_custom_func,
      param_mapping = list(common_params = list("V", "Cl"), filter_params = list()),
      settings = list(iiv = TRUE, plot = FALSE, max_models = 1))
    Condition
      Error in `calculateSearchSpace()`:
      ! object 'nb_params' not found

---

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = valid_custom_func,
      param_mapping = list(common_params = c("V", "Cl"), filter_params = c("ka",
        "Tk0")), settings = list(iiv = TRUE, plot = FALSE, max_models = 1))
    Condition
      Error in `calculateSearchSpace()`:
      ! object 'nb_params' not found

# model_creation_func must have correct signature

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = bad_func_1arg)
    Condition
      Error in `calculateSearchSpace()`:
      ! object 'space' not found

---

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = bad_func_3args)
    Condition
      Error in `calculateSearchSpace()`:
      ! object 'space' not found

# model_creation_func must return valid path

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = func_returns_null,
      settings = tiny_settings())
    Condition
      Error in `calculateSearchSpace()`:
      ! object 'space' not found

---

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = func_returns_bad_path,
      settings = tiny_settings())
    Condition
      Error in `calculateSearchSpace()`:
      ! object 'space' not found

---

    Code
      findModel(project = demo_proj, library = "custom", filters = list(absorption = "firstOrder"),
      algorithm = "exhaustive_search", model_creation_func = func_returns_multiple,
      settings = tiny_settings())
    Condition
      Error in `calculateSearchSpace()`:
      ! object 'space' not found

