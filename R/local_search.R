runLocalSearch <- function(model, model_info, list_of_runs,
                           iiv = NULL, error = NULL, metric_func = NULL,
                           settings = NULL, prev_runs_list = NULL) {

  models_tested <- 0

  best_metric <- model_info$metric

  base_model_info <- model_info
  base_params <- iiv[which(!is.na(iiv))]

  repeat {
    params <- base_params
    model_df <- base_model_info
    improved <- FALSE

    free_params <- names(params)
    if (!is.null(settings$fixed_iiv)) {
      free_params <- setdiff(free_params, names(settings$fixed_iiv))
    }

    for (param in free_params) {

      new_params <- replace(params, param, !params[[param]])
      new_model_df <- replace(model_df, param, !params[[param]])
      new_model_df$metric <- NULL

      # Generate a unique hash for the current model configuration
      model_hash <- create_model_hash(new_model_df)

      # Check if this exact model has been run before
      if (model_hash %in% names(list_of_runs) || !is.null(settings$previous_runs) && model_hash %in% names(prev_runs_list)) {

        if (!is.null(settings$previous_runs) && model_hash %in% names(prev_runs_list)) {
          models_tested <- models_tested + 1
          list_of_runs[[model_hash]] <- prev_runs_list[[model_hash]]
          prev_runs_list[[model_hash]] <- NULL
        }

        metric <- list_of_runs[[model_hash]]$metric

        if (metric < best_metric) {
          best_metric <- metric
          base_model_info <- list_of_runs[[model_hash]]
          base_params <- new_params
          improved <- TRUE
        }

        next
      }

      # If model is new, run it
      runModel(model, iiv = new_params, error = error, obsIDToUse = settings$obsIDToUse, linearization = settings$linearization, param_distributions = settings$param_distributions, obs_mapping = settings$obs_mapping)

      models_tested <- models_tested + 1

      # Calculate the performance metric
      metric <- metric_func()

      # Store results
      new_model_info <- new_model_df
      new_model_info$metric <- metric
      list_of_runs[[model_hash]] <- new_model_info

      if (settings$save_mode == "all")
        lixoftConnectors::saveProject(file.path(
          settings$result_folder,
          "autoBuild",
          paste0("run_", models_tested, ".mlxtran")
        ))

      if (metric < best_metric) {
        best_metric <- metric
        base_model_info <- new_model_info
        base_params <- new_params
        improved <- TRUE

        if (settings$save_mode == "best")
          lixoftConnectors::saveProject(
            file.path(settings$result_folder, "autoBuild", "best_run.mlxtran")
          )
      }
    }

    # if (!improved) break
    break
  }

  return(list(
    list_of_runs = list_of_runs,
    best_model_info = base_model_info,
    models_tested = models_tested,
    prev_runs_list = prev_runs_list
  ))
}

