rbind_fill <- function(lst) {
  all_cols <- Reduce(union, lapply(lst, names))
  lst2 <- lapply(lst, function(x) {
    missing <- setdiff(all_cols, names(x))
    if (length(missing)) x[missing] <- NA
    x[all_cols]
  })
  do.call(rbind, lst2)
}


featureLevels <- function(library) {

  order <- list()
  if (library == "pk") {
    order$delay <- data.frame(
      delay = c("noDelay", "lagTime"),
      level = c(1, 2)
    )
    order$absorption_delay <- data.frame(
      delay = c("noDelay", "noDelay", "noDelay", "lagTime", "lagTime", "lagTime", "noDelay"),
      absorption = c("firstOrder", "zeroOrder", "sigmoid", "firstOrder", "zeroOrder", "sigmoid", "transitCompartments"),
      level = c(1, 1, 2, 2, 2, 3, 3)
    )
    order$distribution <- data.frame(
      distribution = c("1compartment", "2compartments", "3compartments"),
      level = c(1, 2, 3)
    )
    order$elimination <- data.frame(
      elimination = c("linear", "MichaelisMenten", "combined"),
      level = c(1, 2, 3)
    )
    order$bioavailability <- data.frame(
      bioavailability = c("false", "true"),
      level = c(1, 2)
    )
  }

  return(order)
}

getBaseModel <- function(library, filters) {
  levels_map <- featureLevels(library)

  base_model <- list()

  # Handle joint absorption + delay first if both provided
  if (all(c("absorption", "delay") %in% names(filters))) {
    combos <- levels_map$absorption_delay
    valid <- combos[
      (is.na(combos$delay) | combos$delay %in% filters$delay) &
        combos$absorption %in% filters$absorption,
      , drop = FALSE]

    if (nrow(valid) == 0) {
      cli::cli_abort("No valid combination of absorption and delay found in filters.")
    }

    best <- valid[which.min(valid$level), , drop = FALSE]
    base_model$delay <- best$delay
    base_model$absorption <- best$absorption

    simple_features <- setdiff(names(filters), c("absorption", "delay"))
  } else {
    simple_features <- names(filters)
  }

  # Fill the rest feature-by-feature
  common_names <- intersect(names(levels_map), simple_features)
  for (key in common_names) {
    df <- levels_map[[key]]
    # choose the minimal level among the allowed values
    rows <- df[df[[key]] %in% filters[[key]], , drop = FALSE]
    if (nrow(rows) > 0) {
      base_model[[key]] <- rows[[key]][which.min(rows$level)]
    }
  }

  as.data.frame(base_model, stringsAsFactors = FALSE)
}

textify <- function(df, run = NULL) {

  if (is.data.frame(df))
    vals <- unlist(df[1, ], use.names = FALSE)
  else
    vals <- df

  vals <- gsub("([a-z])([A-Z])", "\\1 \\2", vals)

  vals <- gsub("([0-9])([A-Za-z])", "\\1 \\2", vals)
  vals <- gsub("([A-Za-z])([0-9])", "\\1 \\2", vals)

  vals <- tolower(vals)

  text <- paste(vals, collapse = ", ")

  if (!is.null(run))
    text <- paste0(text, "\n", "(run_", run, ")")

  return(text)
}

runDecisionTree <- function(library, filters, requiredFilters, settings, metric_func, initial_func) {

  # Validate that custom library is not used with decision tree
  if (library == "custom") {
    stop("Error: library='custom' is not supported with algorithm='decision_tree'. Use 'exhaustive_search' or 'ACO' instead.", 
         call. = FALSE)
  }

  # --- Initialization ---
  cli::cli_h1("Initializing Decision Tree Search")

  # --- Decision Tree State Variables ---
  list_of_runs <- list()
  prev_runs_list <- NULL

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


  base_model_df <- getBaseModel(library, filters)
  model <- getModelName(library, filters = append(base_model_df, requiredFilters))

  models_tested <- 1

  current_model_df <- base_model_df

  if (settings$iiv) {
    model <- getModelName(library, filters = append(current_model_df, requiredFilters))
    lixoftConnectors::setStructuralModel(model)
    iiv_params <- as.data.frame(as.list(lixoftConnectors::getIndividualParameterModel()$variability[[1]]))
    current_model_df <- do.call(transform, c(list(current_model_df), iiv_params))
  }

  if (settings$error) {
    current_model_df <- transform(current_model_df, error = "combined1")
  }

  model_hash <- create_model_hash(current_model_df)

  if (!is.null(prev_runs_list) && model_hash %in% names(prev_runs_list)) {

    list_of_runs[[model_hash]] <- prev_runs_list[[model_hash]]
    prev_runs_list[[model_hash]] <- NULL
    metric <- list_of_runs[[model_hash]]$metric

  } else {

    runModel(model, obsIDToUse = settings$obsIDToUse, linearization = settings$linearization, initial_func = initial_func, param_distributions = settings$param_distributions)
    metric <- metric_func()

    if (!is.null(settings$table_func))
      table <- settings$table_func()

    save_project(settings$save_mode,
                 path_all  = file.path(settings$result_folder, "autoBuild"),
                 path_best = file.path(settings$result_folder, "autoBuild"),
                 idx = models_tested,
                 improved = TRUE)

  }

  # Store results
  base_metric <- metric


  if (!is.null(settings$table_func))
    list_of_runs[[model_hash]] <- cbind(transform(current_model_df, metric = metric), table)
  else
    list_of_runs[[model_hash]] <- transform(current_model_df, metric = metric)

  levels_map <- featureLevels(library)

  keys <- if (is.null(settings$keys)) sample(names(filters)) else settings$keys
  if (all(c("absorption", "delay") %in% keys)) {
    keys[keys == "absorption"] <- "absorption_delay"
    keys <- keys[keys != "delay"]
  }

  if (settings$plot) {
    base <- data.tree::Node$new(paste0(
      "Base model\n(",
      textify(base_model_df),
      ")"
    ))
    base$metric <- metric
    base$diff <- 0
    data.tree::SetNodeStyle(base)
  }

  path <- c()

  for (key in keys) {
    cli::cli_h3("Optimizing {key}")

    levels_df <- levels_map[[key]]
    feat_cols <- setdiff(names(levels_df), "level")

    base_value <- base_model_df[, feat_cols, drop = FALSE]
    base_value_string <- textify(base_value)

    if (settings$plot) {
      child <- Reduce(`[[`, c(list(base), path))$AddChild(base_value_string)
      child$metric <- base_metric
      child$diff <- 0
      path <- c(path, base_value_string)
    }

    for (i in 1:nrow(levels_df)) {

      candidate <- levels_df[i, feat_cols, drop = FALSE]
      base_value <- base_model_df[, feat_cols, drop = FALSE]

      # must satisfy user filters (respect NA in delay when transit)
      if (!all(mapply(function(col, val) is.na(val) || val %in% filters[[col]], col = feat_cols, val = candidate[1, ]))) {
        next
      }

      # only consider adjacent levels
      base_level <- merge(base_value, levels_df, by = feat_cols, all.x = TRUE)$level
      cand_level <- levels_df$level[i]
      if (!isTRUE(abs(cand_level - base_level) <= 1) || identical(candidate, base_value)) next

      current_model_df <- base_model_df
      current_model_df[, feat_cols] <- candidate

      # transitCompartments => NA delay
      if (identical(library, "pk") && "absorption" %in% names(current_model_df)) {
        if (isTRUE(current_model_df$absorption[1] == "transitCompartments")) current_model_df$delay <- NA
      }

      struct_model_df <- current_model_df

      if (settings$iiv) {
        model <- getModelName(library, filters = append(struct_model_df, requiredFilters))
        lixoftConnectors::setStructuralModel(model)
        iiv_params <- as.data.frame(as.list(lixoftConnectors::getIndividualParameterModel()$variability[[1]]))
        current_model_df <- do.call(transform, c(list(current_model_df), iiv_params))
      }

      if (settings$error) {
        current_model_df <- transform(current_model_df, error = "combined1")
      }

      candidate_string <- textify(candidate)

      if (settings$plot)
        child <- Reduce(`[[`, c(list(base), path[-length(path)]))$AddChild(candidate_string)

      # Generate a unique hash for the current model configuration
      model_hash <- create_model_hash(current_model_df)

      # Check if this exact model has been run before
      if (model_hash %in% names(list_of_runs) || !is.null(prev_runs_list) && model_hash %in% names(prev_runs_list)) {

        if (!is.null(prev_runs_list) && model_hash %in% names(prev_runs_list)) {
          models_tested <- models_tested + 1
          list_of_runs[[model_hash]] <- prev_runs_list[[model_hash]]
          prev_runs_list[[model_hash]] <- NULL
        }

        metric <- list_of_runs[[model_hash]]$metric

        if (settings$plot) {
          child$metric <- metric
          child$diff <- metric - Reduce(`[[`, c(list(base), path[-length(path)]))$metric
        }

      } else {

        model <- getModelName(library, filters = append(struct_model_df, requiredFilters))
        runModel(model, obsIDToUse = settings$obsIDToUse, linearization = settings$linearization, initial_func = initial_func, param_distributions = settings$param_distributions)

        metric <- metric_func()
        models_tested <- models_tested + 1

        save_project(settings$save_mode,
                     path_all  = file.path(settings$result_folder, "autoBuild"),
                     path_best = file.path(settings$result_folder, "autoBuild"),
                     idx = models_tested,
                     improved = metric < base_metric)

        if (!is.null(settings$table_func))
          list_of_runs[[model_hash]] <- cbind(transform(current_model_df, metric = metric), settings$table_func())
        else
          list_of_runs[[model_hash]] <- transform(current_model_df, metric = metric)

      }

      if (settings$plot) {
        child$metric <- metric
        child$diff <- metric - Reduce(`[[`, c(list(base), path[-length(path)]))$metric
      }

      if (base_metric > metric) {
        base_metric <- metric
        base_model_df <- struct_model_df
        path[length(path)] <- candidate_string
      }

    }

    # Save current progress to a file if requested
    if (!is.null(settings$output)) {
      full_results_df <- do.call(rbind, list_of_runs)
      rownames(full_results_df) <- NULL
      tryCatch({
        write.csv(full_results_df, file = settings$output, row.names = FALSE)
      }, error = function(e) cli::cli_alert_warning("Error when writing model list to a file: File is likely open in Excel"))
    }

    cli::cli_alert_success(
      paste0("Done. Selected: ",
             paste(unlist(base_model_df[, feat_cols, drop = FALSE]), collapse = ", ")
      )
    )
  }

  model <- getModelName(library, filters = append(base_model_df, requiredFilters))
  lixoftConnectors::setStructuralModel(model)

  variability <- NA
  if (settings$iiv) {
    cli::cli_h3("Optimizing IIV")

    level <- names(lixoftConnectors::getIndividualParameterModel()$variability)[[1]]

    repeat {
      improved <- FALSE
      params <- lixoftConnectors::getIndividualParameterModel()$variability[[1]]
      params_wo_var <- params[params]
      # Remove fixed IIV params from exploration (TRUE-locked stay on, FALSE-locked already off)
      if (!is.null(settings$fixed_iiv)) {
        params_wo_var <- params_wo_var[!names(params_wo_var) %in% names(settings$fixed_iiv)]
      }
      base_params <- params

      if (settings$plot) {
        child <- Reduce(`[[`, c(list(base), path))$AddChild("No change")
        child$metric <- base_metric
        child$diff <- 0
        path <- c(path, "No change")
      }

      for (param in names(params_wo_var)) {
        iiv_params <- as.data.frame(as.list(params))
        iiv_params[, param] <- FALSE
        runModel(model, iiv = iiv_params, obsIDToUse = settings$obsIDToUse,
                 linearization = settings$linearization, initial_func = initial_func, param_distributions = settings$param_distributions)

        metric <- metric_func()
        models_tested <- models_tested + 1

        save_project(settings$save_mode,
                     path_all  = file.path(settings$result_folder, "autoBuild"),
                     path_best = file.path(settings$result_folder, "autoBuild"),
                     idx = models_tested,
                     improved = metric < base_metric)

        candidate_string <- paste0("No IIV on ", param)

        if (settings$plot) {
          child <- Reduce(`[[`, c(list(base), path[-length(path)]))$AddChild(candidate_string)
          child$metric <- metric
          child$diff <- metric - Reduce(`[[`, c(list(base), path[-length(path)]))$metric
        }

        if (metric < base_metric) {
          base_metric <- metric
          base_params <- replace(params, param, FALSE)
          improved <- TRUE
          path[length(path)] <- candidate_string
        }

        current_model_df <- do.call(transform, c(list(base_model_df), iiv_params))
        if (settings$error) {
          current_model_df <- transform(current_model_df, error = "combined1")
        }

        model_hash <- create_model_hash(current_model_df)
        if (!is.null(settings$table_func))
          list_of_runs[[model_hash]] <- cbind(transform(current_model_df, metric = metric), settings$table_func())
        else
          list_of_runs[[model_hash]] <- transform(current_model_df, metric = metric)
      }

      variability <- base_params

      if (!improved) break
    }

    lixoftConnectors::setIndividualParameterVariability(id = as.list(variability))
    cli::cli_alert_success("Done")
  }

  error <- NA
  if (settings$error) {
    cli::cli_h3("Optimizing error model")

    errors <- c("constant", "proportional", "combined1", "combined2")
    base_error <- error_model <- lixoftConnectors::getContinuousObservationModel()$errorModel[[1]]
    obs_name <- names(lixoftConnectors::getContinuousObservationModel()$errorModel)[[1]]

    if (settings$plot) {
      child <- Reduce(`[[`, c(list(base), path))$AddChild(base_error)
      child$metric <- base_metric
      child$diff <- 0
      path <- c(path, base_error)
    }

    for (error in errors) {
      if (error == error_model)
        next

      error_df <- data.frame(error = error)
      runModel(model, error = error_df, obsIDToUse = settings$obsIDToUse,
               linearization = settings$linearization, initial_func = initial_func, param_distributions = settings$param_distributions)

      metric <- metric_func()
      models_tested <- models_tested + 1

      candidate_string <- error

      if (settings$plot) {
        child <- Reduce(`[[`, c(list(base), path[-length(path)]))$AddChild(candidate_string)
        child$metric <- metric
        child$diff <- metric - Reduce(`[[`, c(list(base), path[-length(path)]))$metric
      }

      save_project(settings$save_mode,
                   path_all  = file.path(settings$result_folder, "autoBuild"),
                   path_best = file.path(settings$result_folder, "autoBuild"),
                   idx = models_tested,
                   improved = metric < base_metric)

      if (metric < base_metric) {
        base_metric <- metric
        base_error <- error
        path[length(path)] <- candidate_string
      }

      current_model_df <- base_model_df
      if (settings$iiv) {
        iiv_params <- as.data.frame(as.list(lixoftConnectors::getIndividualParameterModel()$variability[[1]]))
        current_model_df <- do.call(transform, c(list(current_model_df), iiv_params))
      }
      current_model_df <- transform(current_model_df, error = error)
      model_hash <- create_model_hash(current_model_df)
      if (!is.null(settings$table_func))
        list_of_runs[[model_hash]] <- cbind(transform(current_model_df, metric = metric), settings$table_func())
      else
        list_of_runs[[model_hash]] <- transform(current_model_df, metric = metric)
    }

    error <- base_error

    cli::cli_alert_success("Done")
  }

  if (settings$plot) {
    lapply(seq_along(path), function(i) {
      x <- Reduce(`[[`, c(list(base), path[seq_len(i)]))
      data.tree::SetNodeStyle(x,
                              fillcolor = "#e6f7e6",
                              style = "filled,rounded",
                              fontcolor = "black",
                              inherit = FALSE
      )
    })
    data.tree::SetEdgeStyle(base, label = function(node) round(node$diff, 2))
    data.tree::SetNodeStyle(base, shape = "box", style = "rounded")
    print(plot(base))
  }

  final_results <- rbind_fill(list_of_runs)
  rownames(final_results) <- NULL
  cols <- names(final_results)
  new_order <- c(setdiff(cols, c("error", "metric")), intersect(cols, c("error", "metric")))
  final_results <- final_results[, new_order]

  best_index <- rownames(final_results[which(final_results$metric == min(final_results$metric)), ])

  best_model_path <- if (settings$save_mode == "best")
    file.path(settings$result_folder, "autoBuild", "best_run.mlxtran")
  else if (settings$save_mode == "all")
    file.path(settings$result_folder, "autoBuild", paste0("run_", best_index, ".mlxtran"))
  else
    NULL

  return(list(
    results = final_results[order(final_results$metric), ],
    best_model = final_results[order(final_results$metric), ][1, ],
    best_model_path = best_model_path
  ))
}
