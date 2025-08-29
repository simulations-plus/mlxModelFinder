runTournament <- function(library, filters, requiredFilters, settings, metric_func) {

  # --- Initialization ---
  cli::cli_h1("Initializing Tournament Search")

  # Initialize probabilities for the structural model components
  probas_struct <- probabilities(filters)

  # Initialize probabilities for interindividual variability (IIV)
  if (settings$iiv) {
    probas_iiv <- iiv_probabilities(library)
  }

  # Initialize probabilities for error model
  if (settings$error) {
    probas_error <- probabilities(list(error = c("constant", "proportional", "combined1", "combined2")))
  }

  # --- Tournament State Variables ---
  iteration <- 0
  models_tested <- 0
  n_wins <- settings$N
  best_metric <- Inf
  list_of_runs <- list()
  history_list <- list(tidyProbabilities(probas_struct, iteration = 0))

  # Pre-populate with previous runs, if provided
  if (!is.null(settings$previous_runs) && nrow(settings$previous_runs) > 0) {
    cli::cli_alert_info("Seeding run history with {nrow(settings$previous_runs)} previous results.")
    model_def_cols <- setdiff(colnames(settings$previous_runs), c("metric", "rank"))

    hashes <- apply(settings$previous_runs[, model_def_cols, drop = FALSE], 1, create_model_hash)

    prev_runs_list <- split(settings$previous_runs, seq(nrow(settings$previous_runs)))
    names(prev_runs_list) <- hashes
    list_of_runs <- prev_runs_list
  }

  repeat {

    ## Check if all structural model features converged
    if (sum(probas_struct$probas >= settings$convergence_threshold) == length(colnames(probas_struct$probas)) && (!settings$iiv || sum(probas_iiv$probas >= settings$convergence_threshold) >= 3))
      break

    i <- 0
    iteration <- iteration + 1
    cli::cli_h2("Iteration {iteration}")
    cli::cli_progress_bar(
      paste("Iteration", iteration),
      format = "Running models {cli::pb_bar} {i}/{n_wins} wins",
      total = n_wins,
      clear = FALSE
    )

    winning_model <- NULL

    # This list will hold the results for the current iteration
    iter_results_list <- list()

    while (i < n_wins) {

      # Sample structural models based on their probabilities
      p_iiv <- if (settings$iiv) probas_iiv else NULL
      p_err <- if (settings$error) probas_error else NULL
      current_model <- models(probas_struct, 1, p_iiv = p_iiv, p_err = p_err)

      # Force IIV exploration in early iterations
      if (settings$iiv && iteration <= settings$initial_iiv_forced_iter) {
        current_model$iiv[] <- TRUE
      }

      model <- getModelName("pk", filters = append(current_model$struct, requiredFilters))
      lixoftConnectors::setStructuralModel(model)

      # Combine the structural, IIV, and error components into a single definition
      current_model_df <- as.data.frame(current_model$struct)

      if (settings$iiv) {
        param_names <- lixoftConnectors::getIndividualParameterModel()$name
        current_model$iiv[1, !(colnames(current_model$iiv) %in% param_names)] <- NA
        current_model_df <- cbind(current_model_df, as.data.frame(current_model$iiv))
      }

      if (settings$error) {
        current_model_df <- cbind(current_model_df, data.frame(error = current_model$error[[1, "error"]]))
      }

      # Generate a unique hash for the current model configuration
      model_hash <- create_model_hash(current_model_df)

      # Check if this exact model has been run before
      if (model_hash %in% names(list_of_runs)) {
        model_info <- iter_results_list[[i + 1]] <- list_of_runs[[model_hash]]
        metric <- model_info$metric

        if (is.null(winning_model))
          winning_model <- iter_results_list[[length(iter_results_list)]]

        if (winning_model$metric <= metric) {
          i <- i + 1
        } else {
          i <- 1
          winning_model <- model_info
        }
        cli::cli_progress_update(set = i)
        next
      }

      # If model is new, run it
      iiv_params <- if (settings$iiv) as.data.frame(current_model$iiv) else NULL
      error_model <- if (settings$error) as.data.frame(current_model$error) else NULL

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
      }

      # Store results
      model_info <- current_model_df
      model_info$metric <- metric

      # Add the new result to the main history and the current iteration's list
      list_of_runs[[model_hash]] <- model_info
      iter_results_list[[i + 1]] <- model_info

      if (is.null(winning_model))
        winning_model <- iter_results_list[[length(iter_results_list)]]

      if (winning_model$metric <= metric) {
        i <- i + 1
      } else {
        i <- 1
        winning_model <- model_info
      }
      cli::cli_progress_update(set = i)
    }

    cli::cli_progress_done()

    # --- Post-Iteration Processing ---
    # Efficiently combine the list of results into a single data frame
    iter_models <- do.call(rbind, Filter(Negate(is.null), iter_results_list))
    rownames(iter_models) <- NULL

    # Save current progress to a file if requested
    if (!is.null(settings$output)) {
      full_results_df <- do.call(rbind, list_of_runs)
      rownames(full_results_df) <- NULL
      tryCatch({
        write.csv(full_results_df, file = settings$output, row.names = FALSE)
      }, error = function(e) cli::cli_alert_warning("Error when writing model list to a file: File is likely open in Excel"))
    }

    # Sort the models from the current iteration by performance and rank them
    iter_models <- iter_models[order(iter_models$metric), ]
    iter_models$rank <- 1:nrow(iter_models)

    iter_models <- postprocessModels(iter_models)

    # --- Update Probabilities ---
    probas_struct <- updateProbabilities(probas_struct, iter_models, rho = settings$rho)
    if (settings$iiv && iteration > settings$initial_iiv_forced_iter) {
      probas_iiv <- updateProbabilities(probas_iiv, iter_models, rho = settings$rho)
    }
    if (settings$error) {
      probas_error <- updateProbabilities(probas_error, iter_models, rho = settings$rho)
    }

    # --- Reporting ---
    cli::cli_alert_info("Total number of models tested so far: {models_tested}")
    cli::cli_alert_info("Best model found in iteration {iteration}:")

    best_model_info <- iter_models[1, ]
    ulid <- cli::cli_ul()
    for (colname in colnames(best_model_info)) {
      val <- best_model_info[[colname]]
      if (colname %in% c("rank") || is.na(val)) next
      if (colname == "metric")
        cli::cli_li("{.field {colname}}: {.value {round(val, 2)}}")
      else if (nchar(colname) < 5)
        cli::cli_li("{.field IIV {colname}}: {.value {val}}")
      else
        cli::cli_li("{.field {colname}}: {.value {val}}")
    }
    cli::cli_end(ulid)

    # Store and plot probability history
    history_list[[iteration + 1]] <- tidyProbabilities(probas_struct, iteration)
    plotProbabilities(history_list)
  }

  cli::cli_alert_success("Done")

  final_results <- do.call(rbind, list_of_runs)
  rownames(final_results) <- NULL

  res <- final_results[order(final_results$metric), ]
  res <- res[, colSums(!is.na(res)) > 0]

  return(res)
}
