runACO <- function(library, filters, requiredFilters, settings, metric_func) {

  # --- Initialization ---
  cli::cli_h1("Initializing ACO Search")

  # Initialize probabilities for the structural model components
  probas_struct <- probabilities(filters)
  struct_cols <- colnames(probas_struct$filters)

  # Initialize probabilities for interindividual variability (IIV)
  if (settings$iiv) {
    probas_iiv <- iiv_probabilities(library, favor_true = settings$iiv_favor_true)
    iiv_cols <- colnames(probas_iiv$filters)
  }

  # Initialize probabilities for error model
  if (settings$error) {
    probas_error <- probabilities(list(error = c("constant", "proportional", "combined1", "combined2")))
  }

  # --- ACO State Variables ---
  iteration <- 0
  stagnation_count <- -1
  models_tested <- 0
  n_ants <- settings$N
  best_metric <- Inf
  list_of_runs <- list()
  history_list <- list(tidyProbabilities(probas_struct, iteration = 0))
  if (settings$error)
    history_list[[1]] <- rbind(history_list[[1]], tidyProbabilities(probas_error, iteration))
  if (settings$iiv)
    history_list[[1]] <- rbind(history_list[[1]], tidyProbabilities(probas_iiv, iteration))

  iteration_progress <- list()

  # Pre-populate with previous runs, if provided
  if (!is.null(settings$previous_runs) && nrow(settings$previous_runs) > 0) {
    cli::cli_alert_info("Seeding run history with {nrow(settings$previous_runs)} previous results.")
    model_def_cols <- setdiff(colnames(settings$previous_runs), c("metric", "rank"))

    rows <- split(settings$previous_runs[, model_def_cols, drop = FALSE],
                  seq_len(nrow(settings$previous_runs)))
    hashes <- lapply(rows, create_model_hash)

    prev_runs_list <- split(settings$previous_runs, seq(nrow(settings$previous_runs)))
    names(prev_runs_list) <- hashes
  }

  # --- Main ACO Loop ---
  repeat {
    iteration <- iteration + 1
    cli::cli_h2("Iteration {iteration}")
    cli::cli_progress_bar(
      format = "Running models {cli::pb_bar} | {i}/{n_ants}",
      total = n_ants,
      clear = FALSE
    )

    # Sample structural models based on their probabilities
    p_iiv <- if (settings$iiv) probas_iiv else NULL
    p_err <- if (settings$error) probas_error else NULL
    models <- models(probas_struct, n_ants, p_iiv = p_iiv, p_err = p_err)

    # --- Pre-processing of Sampled Models ---
    # Force IIV exploration in early iterations
    if (settings$iiv && iteration <= settings$initial_iiv_forced_iter) {
      models$iiv[] <- TRUE
    }

    # For certain absorption models, delay is not a meaningful parameter
    if (library == "pk" && requiredFilters$administration %in% c("oral", "oralBolus")) {
      is_special_abs <- models$struct[, "absorption"] == "transitCompartments"
      if (sum(is_special_abs) > 0)
        models$struct[is_special_abs, "delay"] <- NA
    }

    # This list will hold the results for the current iteration
    iter_results_list <- vector("list", n_ants)

    # --- Compile models in a data.frame ---
    if (settings$iiv) {
      for (i in 1:n_ants) {
        model <- getModelName("pk", filters = append(models$struct[i, ], requiredFilters))
        lixoftConnectors::setStructuralModel(model)
        param_names <- lixoftConnectors::getIndividualParameterModel()$name
        models$iiv[i, !(colnames(models$iiv) %in% param_names)] <- NA
      }
    }

    models_df <- as.data.frame(models$struct)

    if (settings$iiv)
      models_df <- cbind(models_df, as.data.frame(models$iiv))

    if (settings$error)
      models_df <- cbind(models_df, data.frame(error = models$error[, "error"]))

    # --- Evaluate Models (Ants) ---
    for (i in 1:n_ants) {

      model <- getModelName("pk", filters = append(models$struct[i, ], requiredFilters))

      current_model_df <- models_df[i, ]

      # Generate a unique hash for the current model configuration
      model_hash <- create_model_hash(current_model_df)

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
          stagnation_count <- -1
        }

        iter_results_list[[i]] <- list_of_runs[[model_hash]]
        cli::cli_progress_update()
        next
      }

      # If model is new, run it
      iiv_params <- if (settings$iiv) as.data.frame(models$iiv)[i, ] else NULL
      error_model <- if (settings$error) as.data.frame(models$error)[i, ] else NULL

      runModel(model, iiv = iiv_params, error = error_model, obsIDToUse = settings$obsIDToUse, linearization = settings$linearization)

      models_tested <- models_tested + 1

      # Calculate the performance metric
      metric <- metric_func()

      if (settings$save_mode == "all")
        lixoftConnectors::saveProject(file.path(
          settings$result_folder,
          "autoBuild",
          paste0("run_", models_tested, ".mlxtran")
        ))

      if (metric < best_metric) {
        best_metric <- metric
        if (settings$save_mode == "best")
          lixoftConnectors::saveProject(
            file.path(settings$result_folder, "autoBuild", "best_run.mlxtran")
          )
        stagnation_count <- -1
      }

      # Store results
      model_info <- current_model_df
      model_info$metric <- metric

      # Add the new result to the main history and the current iteration's list
      list_of_runs[[model_hash]] <- model_info
      iter_results_list[[i]] <- model_info

      cli::cli_progress_update()
    }
    cli::cli_progress_done()

    # --- Search with XGBoost if stagnating ---
    if (settings$xgboost && settings$iiv && stagnation_count == 5) {
      cli::cli_alert_info("Trying to find a good model using XGBoost...")

      list_of_runs_df <- do.call(rbind, list_of_runs)
      rownames(list_of_runs_df) <- NULL
      runs <- train_xgb_and_propose(
        list_of_runs = list_of_runs_df,
        probas_struct = probas_struct,
        probas_iiv = probas_iiv,
        probas_error = probas_error,
        N = settings$xgboost_model_nb,
        requiredFilters = requiredFilters
      )

      for (i in 1:nrow(runs)) {
        current_model_df <- runs[i, ]
        current_model_df$pred_metric <- NULL
        iiv_params <- runs[i, iiv_cols]
        error_model <- if (settings$error) runs[i, "error"] else NULL

        model <- getModelName("pk", filters = append(runs[i, struct_cols], requiredFilters))

        model_hash <- create_model_hash(current_model_df)

        if (!is.null(settings$previous_runs) && model_hash %in% names(prev_runs_list)) {
          list_of_runs[[model_hash]] <- prev_runs_list[[model_hash]]
          prev_runs_list[[model_hash]] <- NULL

          metric <- list_of_runs[[model_hash]]$metric
          iter_results_list[[length(iter_results_list) + 1]] <- list_of_runs[[model_hash]]

          if (metric < best_metric) {
            best_metric <- metric
            stagnation_count <- -1
          }

          next
        }

        runModel(model, iiv = iiv_params, error = error_model, obsIDToUse = settings$obsIDToUse, linearization = settings$linearization)

        models_tested <- models_tested + 1

        # Calculate the performance metric
        metric <- metric_func()

        if (settings$save_mode == "all")
          lixoftConnectors::saveProject(file.path(
            settings$result_folder,
            "autoBuild",
            paste0("run_", models_tested, ".mlxtran")
          ))

        if (metric < best_metric) {
          best_metric <- metric
          if (settings$save_mode == "best")
            lixoftConnectors::saveProject(
              file.path(settings$result_folder, "autoBuild", "best_run.mlxtran")
            )
          stagnation_count <- -1
        }

        # Store results
        model_info <- current_model_df
        model_info$metric <- metric

        # Add the new result to the main history and the current iteration's list
        list_of_runs[[model_hash]] <- model_info
        iter_results_list[[length(iter_results_list) + 1]] <- model_info
      }
    }

    # --- Add global best ---
    if (settings$elitism && iteration > 1)
      iter_results_list[[length(iter_results_list) + 1]] <- best_model_info[, colnames(best_model_info) != "rank"]

    # --- Do local IIV search ---
    if (settings$local_search && settings$iiv && iteration > settings$initial_iiv_forced_iter) {
      cli::cli_alert_info("Running local IIV search on the best model...")

      min_metric <- min(sapply(iter_results_list, function(x) x$metric))
      index_best <- which(sapply(iter_results_list, function(x) x$metric == min_metric))[1]

      model <- getModelName("pk", filters = append(iter_results_list[[index_best]][, struct_cols], requiredFilters))

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
    # Efficiently combine the list of results into a single data frame
    iter_models <- do.call(rbind, Filter(Negate(is.null), iter_results_list))
    rownames(iter_models) <- NULL

    # Save current progress to a file if requested
    if (!is.null(settings$output)) {
      full_results_df <- do.call(rbind, list_of_runs)
      prev_runs <- do.call(rbind, prev_runs_list)
      full_results_df <- rbind(full_results_df, prev_runs)
      rownames(full_results_df) <- NULL
      tryCatch({
        write.csv(full_results_df, file = settings$output, row.names = FALSE)
      }, error = function(e) cli::cli_alert_warning("Error when writing model list to a file: File is likely open in Excel"))
    }

    # Sort the models from the current iteration by performance and rank them
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

    # Store and plot probability history
    history_list[[iteration + 1]] <- tidyProbabilities(probas_struct, iteration)
    if (settings$error)
      history_list[[iteration + 1]] <- rbind(history_list[[iteration + 1]], tidyProbabilities(probas_error, iteration))
    if (settings$iiv)
      history_list[[iteration + 1]] <- rbind(history_list[[iteration + 1]], tidyProbabilities(probas_iiv, iteration))

    if (settings$plot)
      print(plotProbabilities(history_list))

    # Store results per iteration
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

  return(list(
    runs = res,
    iteration_progress = iteration_progress,
    prob_history = history_list
  ))
}

