# --- REMEMBER TO RUN THIS FIRST ---
# library(future)
# library(future.apply)
# plan(multisession, workers = 4)
# ---

runACO <- function(library, filters, requiredFilters, settings, metric_func, initial_func,
                   model_creation_func = NULL,  # NEW
                   param_mapping = NULL) {      # NEW

  project <- lixoftConnectors::getProjectSettings()$project

  # --- Initialization ---
  cli::cli_h1("Initializing ACO Search")

  # Initialize probabilities for the structural model components
  probas_struct <- probabilities(filters)
  struct_cols <- colnames(probas_struct$filters)
  probas_iiv <- if (settings$iiv) iiv_probabilities(library, filters, requiredFilters, favor_true = settings$iiv_favor_true) else NULL
  iiv_cols <- if (settings$iiv && !is.null(probas_iiv)) colnames(probas_iiv$filters) else NULL
  probas_error <- if (settings$error) probabilities(list(error = c("constant", "proportional", "combined1", "combined2"))) else NULL

  # --- ACO State Variables ---
  iteration <- 0
  stagnation_count <- -1
  models_tested <- 0
  best_metric <- Inf
  list_of_runs <- list()

  # Pre-populate with previous runs, if provided
  if (!is.null(settings$previous_runs) && nrow(settings$previous_runs) > 0) {
    cli::cli_alert_info("Seeding run history with {nrow(settings$previous_runs)} previous results.")
    model_def_cols <- setdiff(colnames(settings$previous_runs), c("metric", "rank"))

    rows <- split(settings$previous_runs[, model_def_cols, drop = FALSE],
                  seq_len(nrow(settings$previous_runs)))
    hashes <- lapply(rows, create_model_hash)

    prev_runs_list <- split(settings$previous_runs, seq(nrow(settings$previous_runs)))
    names(prev_runs_list) <- hashes
  } else {
    prev_runs_list <- NULL
  }

  # --- History Tracking ---
  log_history <- function(iter) {
    history <- tidyProbabilities(probas_struct, iter)
    if (settings$error) history <- rbind(history, tidyProbabilities(probas_error, iter))
    if (settings$iiv) history <- rbind(history, tidyProbabilities(probas_iiv, iter))
    return(history)
  }
  history_list <- list(log_history(0))
  iteration_progress <- list()

  # --- Main ACO Loop ---
  repeat {
    iteration <- iteration + 1
    cli::cli_h2("Iteration {iteration}")

    # --- Generate Models (Ants) ---
    models <- models(probas_struct, settings$N, p_iiv = probas_iiv, p_err = probas_error)

    # Pre-processing of Sampled Models
    if (settings$iiv && iteration <= settings$initial_iiv_forced_iter) {
      models$iiv[] <- TRUE
    }

    # For certain absorption models, delay is not meaningful (only for standard PK)
    if (library == "pk" && !is.null(requiredFilters) && requiredFilters$administration %in% c("oral", "oralBolus")) {
      is_special_abs <- models$struct[, "absorption"] == "transitCompartments"
      if (sum(is_special_abs) > 0)
        models$struct[is_special_abs, "delay"] <- NA
    }

    # Combine all model components into a single data frame
    if (settings$iiv) {
      for (i in 1:settings$N) {
        # Prepare filters for getModelName
        if (library == "custom" || is.null(requiredFilters) || length(requiredFilters) == 0) {
          model_filters <- as.list(models$struct[i, ])
        } else {
          model_filters <- append(models$struct[i, ], requiredFilters)
        }
        
        model <- getModelName(
          library,  # Use library parameter
          filters = model_filters,
          custom_func = model_creation_func  # Pass custom function
        )
        lixoftConnectors::setStructuralModel(model)
        
        # Get parameter names from the model
        param_names <- lixoftConnectors::getIndividualParameterModel()$name
        
        # For custom models, also check against param_mapping
        if (library == "custom" && !is.null(param_mapping)) {
          expected_params <- getParamListFromMapping(
            as.list(models$struct[i, ]), 
            param_mapping
          )
          param_names <- intersect(param_names, expected_params)
        }
        
        models$iiv[i, !(colnames(models$iiv) %in% param_names)] <- NA
      }
    }

    models_df <- as.data.frame(models$struct)

    if (settings$iiv)
      models_df <- cbind(models_df, as.data.frame(models$iiv))

    if (settings$error)
      models_df <- cbind(models_df, data.frame(error = models$error[, "error"]))

    # --- Evaluate Models (Parallelized) ---

    # 1. DEFINE A "PURE" WORKER FUNCTION
    # (This is safe to define outside the loop as it doesn't modify state)
    run_single_model_parallel <- function(model_df_row, project) {
      lixoftConnectors::initializeLixoftConnectors()
      lixoftConnectors::loadProject(project)
      model_struct <- model_df_row[, struct_cols, drop = FALSE]
      
      model_name <- getModelName(
        library,  # Use library parameter
        filters = append(as.list(model_struct), requiredFilters),
        custom_func = model_creation_func  # Pass custom function
      )
      
      iiv_params <- if (settings$iiv) model_df_row[, iiv_cols, drop = FALSE] else NULL
      error_model <- if (settings$error) model_df_row[, "error", drop = FALSE] else NULL
      runModel(model_name, iiv = iiv_params, error = error_model, obsIDToUse = settings$obsIDToUse, linearization = settings$linearization, initial_func = initial_func)
      metric <- metric_func()
      result <- model_df_row
      result$metric <- metric
      if (!is.null(settings$table_func))
        result <- cbind(result, settings$table_func())
      # Read-only check against best_metric
      improved <- metric < best_metric
      return(list(result = result, improved = improved))
    }

    # 2. SPLIT: Pre-filter models that need to be run
    models_to_run_list <- list()
    iter_results_list <- list()

    cli::cli_alert_info("Checking {nrow(models_df)} generated models against history...")
    for (i in 1:nrow(models_df)) {
      model_df_row <- models_df[i, , drop = FALSE]
      model_hash <- create_model_hash(model_df_row)

      if (model_hash %in% names(list_of_runs) || model_hash %in% names(prev_runs_list)) {
        if (model_hash %in% names(prev_runs_list)) {
          # FIX: Use = instead of <<-
          list_of_runs[[model_hash]] <- prev_runs_list[[model_hash]]
          prev_runs_list[[model_hash]] <- NULL
        }

        result <- list_of_runs[[model_hash]]
        iter_results_list[[model_hash]] <- result

        if (result$metric < best_metric) {
          # FIX: Use = instead of <<-
          best_metric <- result$metric
          stagnation_count <- -1
        }

      } else {
        models_to_run_list[[model_hash]] <- model_df_row
      }
    }

    # 3. APPLY: Run new models in parallel
    if (length(models_to_run_list) > 0) {
      cli::cli_alert_info("Running {length(models_to_run_list)} new models...")

      parallel_results <- future.apply::future_lapply(
        models_to_run_list,
        function(row) {
          res <- run_single_model_parallel(row, project)
          return(res)
        },
        future.seed = TRUE
      )
      cli::cli_alert_success("Execution complete.")

      # 4. COMBINE: Integrate parallel results back into the main state

      for (model_hash in names(parallel_results)) {
        parallel_run <- parallel_results[[model_hash]]
        result <- parallel_run$result
        improved <- parallel_run$improved

        models_tested <- models_tested + 1

        if (improved) {
          best_metric <- result$metric
          stagnation_count <- -1
        }

        save_project(
          mode = settings$save_mode,
          path_all = file.path(settings$result_folder, "autoBuild"),
          path_best = file.path(settings$result_folder, "autoBuild"),
          idx = models_tested,
          improved = improved
        )

        list_of_runs[[model_hash]] <- result
        iter_results_list[[model_hash]] <- result
      }
    } else {
      cli::cli_alert_info("No new models to run in this iteration (all found in history).")
    }

    iter_results_list <- unname(Filter(Negate(is.null), iter_results_list))

    # --- Search with XGBoost if stagnating ---
    if (settings$xgboost && settings$iiv && settings$error && stagnation_count == 7) {
      cli::cli_alert_info("Stagnation detected. Searching with XGBoost...")

      list_of_runs_df <- do.call(rbind, list_of_runs)
      rownames(list_of_runs_df) <- NULL

      xgb_runs <- train_xgb_and_propose(
        list_of_runs = list_of_runs_df,
        probas_struct = probas_struct, probas_iiv = probas_iiv, probas_error = probas_error,
        N = settings$xgboost_model_nb, requiredFilters = requiredFilters
      )

      xgb_runs$pred_score <- NULL

      cli::cli_alert_info("Evaluating {nrow(xgb_runs)} models from XGBoost...")
      xgb_results_list <- vector("list", nrow(xgb_runs))
      xgb_models_to_run <- list()

      # 1. SPLIT (XGBoost)
      for (i in 1:nrow(xgb_runs)) {
        model_df_row <- xgb_runs[i, , drop = FALSE]
        model_hash <- create_model_hash(model_df_row)
        if (model_hash %in% names(list_of_runs) || model_hash %in% names(prev_runs_list)) {
          if (model_hash %in% names(prev_runs_list)) {
            list_of_runs[[model_hash]] <- prev_runs_list[[model_hash]]
            prev_runs_list[[model_hash]] <- NULL
          }
          result <- list_of_runs[[model_hash]]
          xgb_results_list[[i]] <- result
          if (result$metric < best_metric) {
            best_metric <- result$metric
            stagnation_count <- -1
          }
        } else {
          xgb_models_to_run[[model_hash]] <- model_df_row
        }
      }

      # 2. APPLY (XGBoost)
      if (length(xgb_models_to_run) > 0) {
        cli::cli_alert_info("Running {length(xgb_models_to_run)} new XGBoost models in parallel...")
        xgb_parallel_results <- future.apply::future_lapply(
          xgb_models_to_run,
          function(row) {
            res <- run_single_model_parallel(row, project)
            return(res)
          },
          future.seed = TRUE
        )

        # 3. COMBINE (XGBoost)
        for (model_hash in names(xgb_parallel_results)) {
          parallel_run <- xgb_parallel_results[[model_hash]]
          result <- parallel_run$result
          improved <- parallel_run$improved

          models_tested <- models_tested + 1
          if (improved) {
            best_metric <- result$metric
            stagnation_count <- -1
          }
          save_project(
            mode = settings$save_mode,
            path_all = file.path(settings$result_folder, "autoBuild"),
            path_best = file.path(settings$result_folder, "autoBuild"),
            idx = models_tested,
            improved = improved
          )
          list_of_runs[[model_hash]] <- result
          iter_results_list[[length(iter_results_list) + 1]] <- result
        }
      }

      existing_xgb_results <- Filter(Negate(is.null), xgb_results_list)
      if (length(existing_xgb_results) > 0) {
        iter_results_list <- c(iter_results_list, existing_xgb_results)
      }
    }

    # --- Add global best ---
    if (settings$elitism && iteration > 1)
      iter_results_list[[length(iter_results_list) + 1]] <- best_model_info[, colnames(best_model_info) != "rank"]

    # --- Do local IIV search ---
    # (This section remains unchanged, as it was already correct)
    if (settings$local_search && settings$iiv && iteration > settings$initial_iiv_forced_iter) {
      cli::cli_alert_info("Running local IIV search on the best model...")

      min_metric <- min(sapply(iter_results_list, function(x) x$metric))
      index_best <- which(sapply(iter_results_list, function(x) x$metric == min_metric))[1]

      model <- getModelName(
        library,  # Use library parameter instead of hardcoded "pk"
        filters = append(iter_results_list[[index_best]][, struct_cols], requiredFilters),
        custom_func = model_creation_func  # Pass custom function
      )

      iiv_params <- iter_results_list[[index_best]][, iiv_cols]
      error_model <- if (settings$error) iter_results_list[[index_best]][, "error"] else NULL

      if (is.null(settings$previous_runs)) prev_runs_list <- NULL

      res <- runLocalSearch(
        model = model,
        model_info = iter_results_list[[index_best]],
        list_of_runs = list_of_runs,
        iiv = iiv_params,
        error = error_model,
        settings = settings,
        metric_func = metric_func,
        prev_runs_list = prev_runs_list
      )

      # This is the crucial part: runLocalSearch returns *new* lists,
      # which replace the old ones. My fix (using =) makes this safe.
      list_of_runs <- res$list_of_runs
      prev_runs_list <- res$prev_runs_list
      models_tested <- models_tested + res$models_tested

      if (res$best_model_info$metric < min_metric) {
        best_metric <- res$best_model_info$metric
        stagnation_count <- -1
        best_model_info <- res$best_model_info
      }

      iter_results_list[[index_best]] <- res$best_model_info
    }

    # --- Post-Iteration Processing ---
    iter_models <- do.call(rbind, Filter(Negate(is.null), iter_results_list))
    rownames(iter_models) <- NULL

    if (!is.null(settings$output)) {
      full_results_df <- do.call(rbind, list_of_runs)
      prev_runs <- do.call(rbind, prev_runs_list)
      full_results_df <- rbind(full_results_df, prev_runs)
      rownames(full_results_df) <- NULL
      tryCatch({
        write.csv(full_results_df, file = settings$output, row.names = FALSE)
      }, error = function(e) cli::cli_alert_warning("Error when writing model list to a file: File is likely open in Excel"))
    }

    iter_models <- iter_models[order(iter_models$metric), ]
    iter_models$rank <- 1:nrow(iter_models)
    iter_models <- postprocessModels(iter_models)

    # --- Update Pheromones (Probabilities) ---
    probas_struct <- updateProbabilities(probas_struct, iter_models, rho = settings$rho, clip = settings$clip)
    if (settings$iiv && iteration > settings$initial_iiv_forced_iter) {
      probas_iiv <- updateProbabilities(probas_iiv, iter_models, rho = settings$rho, clip = settings$clip)
    }
    if (settings$error) {
      probas_error <- updateProbabilities(probas_error, iter_models, rho = settings$rho, clip = settings$clip)
    }

    # --- Reporting ---
    cli::cli_alert_info("Total number of models tested so far: {models_tested}")
    if (settings$elitism)
      cli::cli_alert_info("Best model found so far:")
    else
      cli::cli_alert_info("Best model found in iteration {iteration}:")

    best_model_info <- iter_models[1, ]
    ulid <- cli::cli_ul()
    for (colname in colnames(best_model_info)) {
      val <- best_model_info[[colname]]
      if (colname %in% c("rank") || is.na(val)) next
      if (colname == "metric")
        cli::cli_li("{.field {colname}}: {.value {round(val, 2)}}")
      else if (colname %in% colnames(models$iiv))
        cli::cli_li("{.field IIV {colname}}: {.value {val}}")
      else
        cli::cli_li("{.field {colname}}: {.value {val}}")
    }
    cli::cli_end(ulid)

    history_list[[iteration + 1]] <- log_history(iteration)

    if (settings$plot)
      print(plotProbabilities(history_list))

    iteration_progress[[iteration]] <- cbind(
      data.frame(iteration = iteration, models_tested = models_tested),
      best_model_info
    )

    # --- Check for Convergence ---
    stagnation_count <- stagnation_count + 1

    if (stagnation_count >= 8) {
      cli::cli_alert_success("Convergence criteria met. {stagnation_count} iterations without improvement.")
      break
    }
    if (models_tested > settings$max_models) {
      cli::cli_alert_warning("Maximum number of tested models reached.")
      break
    }
  }

  cli::cli_alert_success("ACO search complete.")

  final_results <- do.call(rbind, list_of_runs)
  rownames(final_results) <- NULL

  res <- final_results[order(final_results$metric), ]
  res <- res[, colSums(!is.na(res)) > 0]

  iteration_progress <- do.call(rbind, iteration_progress)
  rownames(iteration_progress) <- NULL
  iteration_progress$rank <- NULL
  iteration_progress <- iteration_progress[, colSums(!is.na(iteration_progress)) > 0]

  best_index <- rownames(res[which(res$metric == min(res$metric)), ])

  best_model_path <- if (settings$save_mode == "best")
    file.path(settings$result_folder, "autoBuild", "best_run.mlxtran")
  else if (settings$save_mode == "all")
    file.path(settings$result_folder, "autoBuild", paste0("run_", best_index, ".mlxtran"))
  else
    NULL

  return(list(
    results = res,
    best_model = res[1, ],
    best_model_path = best_model_path,
    iteration_progress = iteration_progress,
    prob_history = history_list
  ))
}
