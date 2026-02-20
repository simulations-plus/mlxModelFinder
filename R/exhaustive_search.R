# --- REMEMBER TO RUN THIS FIRST ---
# library(future)
# library(future.apply)
# library(progressr)
# library(future.batchtools) # (If using Slurm)
#
# plan(multisession, workers = 4)
# # OR
# # plan(batchtools_slurm, resources = list(walltime = "1:00:00", ...))
# ---

runExhaustiveSearch <- function(library, filters, requiredFilters, settings, metric_func, initial_func,
                                model_creation_func = NULL,  # NEW
                                param_mapping = NULL) {      # NEW

  project <- lixoftConnectors::getProjectSettings()$project

  # --- 1. "SPLIT" (Generate Job List) ---
  # First, create the flat list of all possible jobs to run.

  cli::cli_h1("Generating full grid search job list...")

  all_models <- expand.grid(filters, stringsAsFactors = FALSE)

  if ("transitCompartments" %in% all_models$absorption)
    all_models[all_models$absorption == "transitCompartments", ]$delay <- NA

  all_models <- unique(all_models)

  if (settings$error) {
    errors <- data.frame(error = c("constant", "proportional", "combined1", "combined2"))
    all_models <- merge(all_models, errors)
  }

  job_list <- list()

  # Determine parameter list based on library type
  if (library == "custom" && !is.null(param_mapping)) {
    # For custom models, get all possible parameters from mapping
    param_list <- unique(c(
      param_mapping$common_params,
      unlist(param_mapping$filter_params)
    ))
  } else {
    # Use standard PK parameter list
    param_list <- c("ka", "Tk0", "Tk0s", "Tlag", "Ktr", "Mtt",
                    "V", "V1", "Cl", "Q", "Q2", "V2", "Q3", "V3",
                    "k", "k12", "k21", "k13", "k31",
                    "Vm", "Km", "F")
  }

  # Create a blank data frame for storing IIV info
  param_cols <- as.data.frame(lapply(param_list, function(x) NA))
  names(param_cols) <- param_list

  # Loop sequentially to build the job list
  for (i in seq_len(nrow(all_models))) {

    model_row <- all_models[i, , drop = FALSE]
    struct_cols <- setdiff(names(model_row), "error")

    # Use custom or standard getModelName
    # For custom libraries, requiredFilters might be NULL
    if (library == "custom" || is.null(requiredFilters) || length(requiredFilters) == 0) {
      model_filters <- as.list(model_row[, struct_cols, drop = FALSE])
    } else {
      model_filters <- append(model_row[, struct_cols], requiredFilters)
    }
    
    model_name <- getModelName(
      library, 
      filters = model_filters,
      custom_func = model_creation_func  # Pass through custom function
    )

    # This is the stateful part that must be run sequentially
    lixoftConnectors::setStructuralModel(model_name)

    error_model <- if (settings$error) model_row[, "error", drop = FALSE] else NULL

    # Get IIV combinations
    if (settings$iiv) {
      # For custom models with param_mapping, determine which params are in this model
      if (library == "custom" && !is.null(param_mapping)) {
        # Get parameters for this specific filter combination
        model_params <- getParamListFromMapping(
          model_filters,
          param_mapping
        )

        # Get actual parameters from Monolix
        available_params <- lixoftConnectors::getIndividualParameterModel()$name

        # Intersect to get valid IIV parameters for this model
        valid_iiv_params <- intersect(model_params, available_params)

        if (length(valid_iiv_params) > 0) {
          variability_df <- as.data.frame(as.list(
            lixoftConnectors::getIndividualParameterModel()$variability$id
          ))

          # Only select columns that actually exist in the data frame
          existing_cols <- intersect(valid_iiv_params, names(variability_df))

          if (length(existing_cols) > 0) {
            default_iiv <- variability_df[, existing_cols, drop = FALSE]
            iiv_combos <- apply_fixed_iiv(default_iiv, settings$fixed_iiv)
          } else {
            iiv_combos <- data.frame(row.names = 1)
          }
        } else {
          iiv_combos <- data.frame(row.names = 1)
        }
      } else {
        # Standard PK logic (existing code)
        default_iiv <- as.data.frame(as.list(
          lixoftConnectors::getIndividualParameterModel()$variability$id
        ))
        if (ncol(default_iiv) > 0) {
          iiv_combos <- apply_fixed_iiv(default_iiv, settings$fixed_iiv)
        } else {
          iiv_combos <- data.frame(row.names = 1)
        }
      }
    } else {
      iiv_combos <- as.data.frame(as.list(
        lixoftConnectors::getIndividualParameterModel()$variability$id
      ))
      # Apply fixed_iiv even when not searching IIV
      if (!is.null(settings$fixed_iiv)) {
        for (p in names(settings$fixed_iiv)) {
          if (p %in% names(iiv_combos)) {
            iiv_combos[[p]] <- settings$fixed_iiv[[p]]
          }
        }
      }
    }

    # Add each IIV combination as a separate job
    for (j in 1:nrow(iiv_combos)) {
      iiv_model <- iiv_combos[j, , drop = FALSE]

      # Store the base info (struct model + error)
      if (!settings$iiv) {
        model_info_base <- model_row
      } else {
        # Add IIV columns, populating only the relevant ones
        model_info_base <- cbind(model_row, param_cols)
        # Only set values for parameters that exist in this IIV combo
        if(ncol(iiv_model) > 0) {
          model_info_base[, names(iiv_model)] <- iiv_model
        }
      }

      # Create the job packet
      job <- list(
        model_name = model_name,
        iiv_model = iiv_model,
        error_model = error_model,
        model_info_base = model_info_base, # This has all info for the final table
        job_index = length(job_list) + 1
      )

      job_list[[length(job_list) + 1]] <- job
    }
  }

  cli::cli_alert_success("Job list complete. {length(job_list)} total models to run.")

  # --- 2. "APPLY" (Define and Run Worker) ---
  # This "pure" function will be run in parallel on each worker.
  # It takes one job packet and returns a single data.frame row.

  run_exhaustive_worker <- function(job, project) {

    lixoftConnectors::initializeLixoftConnectors()
    lixoftConnectors::loadProject(project)

    # Run the model (this is the time-consuming part)
    runModel(job$model_name,
             iiv = job$iiv_model,
             error = job$error_model,
             obsIDToUse = settings$obsIDToUse,
             linearization = settings$linearization,
             initial_func = initial_func,
             param_distributions = settings$param_distributions)

    metric <- metric_func()

    # Save project if save_mode == "all"
    if (identical(settings$save_mode, "all")) {
      lixoftConnectors::saveProject(
        file.path(settings$result_folder, "autoBuild",
                  paste0("run_", job$job_index, ".mlxtran"))
      )
    }

    # Populate the result row
    model_info <- job$model_info_base
    model_info$metric <- metric

    if (!is.null(settings$table_func)) {
      result_row <- cbind(model_info, settings$table_func())
    } else {
      result_row <- model_info
    }

    return(result_row)
  }

  # --- Run all jobs in parallel ---
  cli::cli_h1("Running exhaustive search...")

  all_results_list <- NULL

  p <- progressr::progressor(steps = length(job_list))

  all_results_list <- future.apply::future_lapply(
    job_list,
    function(job) {
      res <- run_exhaustive_worker(job, project)
      p() # Update progress
      return(res)
    },
    future.seed = TRUE
  )

  cli::cli_alert_success("All models complete.")

  # --- 3. "COMBINE" (Process Results) ---
  cli::cli_h1("Combining results and saving...")

  results_df <- do.call(rbind, all_results_list)
  rownames(results_df) <- NULL

  # Sort to find the best
  results_df <- results_df[order(results_df$metric), ]
  best_model_df <- results_df[1, , drop = FALSE]

  # --- Save final CSV (replaces the inefficient in-loop save) ---
  if (!is.null(settings$output)) {
    tryCatch({
      write.csv(results_df, file = settings$output, row.names = FALSE)
      cli::cli_alert_success("Full results saved to {settings$output}")
    }, error = function(e) {
      cli::cli_alert_warning("Could not save results CSV. Is the file open?")
    })
  }

  # --- Save the "best" model project file (re-run once) ---
  if (settings$save_mode %in% c("best", "all") && nrow(best_model_df) > 0) {
    cli::cli_alert_info("Re-running best model to save project file...")

    # Re-build the components for the best model
    struct_cols <- names(filters)
    best_struct_info <- best_model_df[, struct_cols, drop = FALSE]
    best_model_name <- getModelName(
      library,  # Use the library parameter, not hardcoded "pk"
      filters = append(as.list(best_struct_info), requiredFilters),
      custom_func = model_creation_func  # Pass through custom function
    )

    best_error_model <- if (settings$error) best_model_df[, "error", drop = FALSE] else NULL

    if (settings$iiv) {
      iiv_cols_in_result <- intersect(param_list, names(best_model_df))
      best_iiv_model <- best_model_df[, iiv_cols_in_result, drop = FALSE]
      # Drop NAs (parameters not in this model)
      best_iiv_model <- best_iiv_model[, colSums(!is.na(best_iiv_model)) > 0, drop = FALSE]
    } else {
      lixoftConnectors::setStructuralModel(best_model_name)
      best_iiv_model <- as.data.frame(as.list(lixoftConnectors::getIndividualParameterModel()$variability$id))
      if (!is.null(settings$fixed_iiv)) {
        for (p in names(settings$fixed_iiv)) {
          if (p %in% names(best_iiv_model)) {
            best_iiv_model[[p]] <- settings$fixed_iiv[[p]]
          }
        }
      }
    }

    # Run it one last time
    runModel(best_model_name,
             iiv = best_iiv_model,
             error = best_error_model,
             obsIDToUse = settings$obsIDToUse,
             linearization = settings$linearization,
             initial_func = initial_func,
             param_distributions = settings$param_distributions)

    # And save it
    lixoftConnectors::saveProject(
      file.path(settings$result_folder, "autoBuild", "best_run.mlxtran")
    )
    cli::cli_alert_success("Best model project file saved.")
  }

  cli::cli_alert_success("Done.")

  return(list(results = results_df,
              best_model = best_model_df))
}
