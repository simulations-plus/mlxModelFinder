featureLevels <- function(library) {

  order <- list()
  if (library == "pk") {
    order$delay <- data.frame(
      delay = c("noDelay", "lagTime"),
      level = c(1, 2)
    )
    order$absorption_delay <- data.frame(
      delay = c("noDelay", "noDelay", "noDelay", "lagTime", "lagTime", "lagTime", "transitCompartments"),
      absorption = c("firstOrder", "zeroOrder", "sigmoid", "firstOrder", "zeroOrder", "sigmoid", "firstOrder"),
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

  order <- featureLevels(library)

  base_model <- list()

  if (all(c("absorption", "delay") %in% names(filters))) {
    valid_combinations <- order$absorption_delay[
      order$absorption_delay$delay %in% filters$delay &
        order$absorption_delay$absorption %in% filters$absorption,
    ]
    if (nrow(valid_combinations) > 0) {
      best_combination <- valid_combinations[which.min(valid_combinations$level), ]
      base_model$delay <- best_combination$delay
      base_model$absorption <- best_combination$absorption
    } else {
      cli::cli_abort("No valid combination of absorption and delay found in filters.")
    }
    simple_features <- setdiff(names(filters), c("absorption", "delay"))
  } else {
    simple_features <- names(filters)
  }

  common_names <- intersect(names(order), simple_features)

  for (key in common_names) {
    levels <- order[[key]]
    valid_rows <- levels[levels[[key]] %in% filters[[key]], ]
    if(nrow(valid_rows) > 0) {
      base_model[[key]] <- valid_rows[which.min(valid_rows$level), key]
    }
  }

  return(as.data.frame(base_model))
}

runDecisionTree <- function(library, filters, requiredFilters, settings, metric_func) {

  # --- Initialization ---
  cli::cli_h1("Initializing Decision Tree Search")

  # --- Decision Tree State Variables ---
  list_of_runs <- list()

  # Pre-populate with previous runs, if provided
  if (!is.null(settings$previous_runs) && nrow(settings$previous_runs) > 0) {
    cli::cli_alert_info("Seeding run history with {nrow(settings$previous_runs)} previous results.")
    model_def_cols <- setdiff(colnames(settings$previous_runs), c("metric", "rank"))

    hashes <- apply(settings$previous_runs[, model_def_cols, drop = FALSE], 1, create_model_hash)

    prev_runs_list <- split(settings$previous_runs, seq(nrow(settings$previous_runs)))
    names(prev_runs_list) <- hashes
    list_of_runs <- prev_runs_list
  }

  base_model_df <- getBaseModel(library, filters)
  model_hash <- create_model_hash(base_model_df)
  model <- getModelName("pk", filters = append(base_model_df, requiredFilters))

  runModel(model, obsIDToUse = settings$obsIDToUse)

  # Calculate the performance metric
  metric <- metric_func()

  models_tested <- 1

  if (settings$save_mode == "all")
    lixoftConnectors::saveProject(file.path(
      settings$result_folder,
      "autoBuild",
      paste0("run_", models_tested, ".mlxtran")
    ))

  if (settings$save_mode == "best")
    lixoftConnectors::saveProject(
      file.path(settings$result_folder, "autoBuild", "best_run.mlxtran")
    )

  # Store results
  model_info <- base_model_df
  model_info$metric <- base_metric <- metric

  # Add the new result to the main history and the current iteration's list
  list_of_runs[[model_hash]] <- model_info

  order <- featureLevels(library)

  # random permutation of filters column names
  if (is.null(settings$keys))
    keys <- sample(names(filters))
  else
    keys <- settings$keys

  if (all(c("absorption", "delay") %in% keys)) {
    keys[keys == "absorption"] <- "absorption_delay"
    keys <- keys[keys != "delay"]
  }

  base <- data.tree::Node$new("Initial model")
  base$metric <- metric
  base$diff <- 0

  path <- c()

  for (key in keys) {
    cli::cli_h3("Optimizing {key}")

    levels <- order[[key]]

    base_value <- base_model_df[, colnames(base_model_df) %in% names(levels)]
    base_value_string <- paste(unlist(base_value), collapse = "_")

    child <- Reduce(`[[`, c(list(base), path))$AddChild(base_value_string)
    child$metric <- base_metric
    child$diff <- 0
    path <- c(path, base_value_string)

    for (i in 1:nrow(levels)) {

      features <- setdiff(names(levels), "level")
      value <- levels[i, features, drop = FALSE]

      all_match <- all(sapply(features, function(col) {
        val <- value[[col]]
        val %in% filters[[col]]
      }))

      if (!all_match) next

      base_value <- base_model_df[, features, drop = FALSE]
      base_level <- merge(base_value, levels, by = features, all.x = TRUE)$level

      level <- levels[i, "level"]

      if (identical(value, base_value) || abs(level - base_level) > 1)
        next

      current_model_df <- base_model_df
      current_model_df[, features] <- value

      value_string <- paste(unlist(value), collapse = "_")

      child <- Reduce(`[[`, c(list(base), path[-length(path)]))$AddChild(value_string)

      if (library == "pk") {
        is_special_delay <- current_model_df[, "delay"] == "transitCompartments"
        if (is_special_delay) current_model_df[is_special_delay, "absorption"] <- NA
      }

      # Generate a unique hash for the current model configuration
      model_hash <- create_model_hash(current_model_df)

      # Check if this exact model has been run before
      if (model_hash %in% names(list_of_runs)) {
        metric <- list_of_runs[[model_hash]]$metric
        child$metric <- metric
        child$diff <- metric - Reduce(`[[`, c(list(base), path[-length(path)]))$metric
        next
      }

      model <- getModelName("pk", filters = append(current_model_df, requiredFilters))
      runModel(model, obsIDToUse = settings$obsIDToUse)

      # Calculate the performance metric
      metric <- metric_func()

      models_tested <- models_tested + 1

      if (settings$save_mode == "all")
        lixoftConnectors::saveProject(file.path(
          settings$result_folder,
          "autoBuild",
          paste0("run_", models_tested, ".mlxtran")
        ))

      if (metric < base_metric) {
        if (settings$save_mode == "best")
          lixoftConnectors::saveProject(
            file.path(settings$result_folder, "autoBuild", "best_run.mlxtran")
          )
      }

      # Store results
      model_info <- current_model_df
      model_info$metric <- metric

      child$metric <- metric
      child$diff <- metric - Reduce(`[[`, c(list(base), path[-length(path)]))$metric

      # Add the new result to the main history and the current iteration's list
      list_of_runs[[model_hash]] <- model_info

      if (base_metric > metric) {
        base_metric <- metric
        base_model_df <- current_model_df
        path[length(path)] <- value_string
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

    cli::cli_alert_success("Done")
  }

  nodes <- data.tree::Traverse(base)
  data.tree::Do(nodes,
     function(node) {
       data.tree::SetNodeStyle(node,
                    fillcolor = ifelse(node$diff <= 0, "#e6f7e6", "white"),
                    style = "filled",
                    fontcolor = "black",
                    inherit = FALSE
       )
     }
  )
  data.tree::SetEdgeStyle(base, label = function(node) ifelse(is.null(node$diff), "0", round(node$diff, 2)))
  print(plot(base))

  final_results <- do.call(rbind, list_of_runs)
  rownames(final_results) <- NULL

  model <- getModelName("pk", filters = append(base_model_df, requiredFilters))
  lixoftConnectors::setStructuralModel(model)

  variability <- NA
  if (settings$iiv) {
    cli::cli_h3("Optimizing IIV")

    project_path <- gsub("\\\\", "/", tempfile(pattern = "project", fileext = ".mlxtran"))
    lixoftConnectors::saveProject(project_path)
    project <- Rsmlx::buildVar(project = project_path, print = FALSE)
    variability <- project$variability.model
    lixoftConnectors::loadProject(project$project)

    cli::cli_alert_success("Done")
  }

  error <- NA
  if (settings$error) {
    cli::cli_h3("Optimizing error model")

    project_path <- gsub("\\\\", "/", tempfile(pattern = "project", fileext = ".mlxtran"))
    lixoftConnectors::saveProject(project_path)
    suppressMessages(
      project <- Rsmlx::buildmlx(project = project_path, model = "residualError", print = FALSE)
    )
    lixoftConnectors::loadProject(project$project)
    error <- project$error.model

    cli::cli_alert_success("Done")
  }

  if (settings$save_mode == "best")
    lixoftConnectors::saveProject(
      file.path(settings$result_folder, "autoBuild", "best_run.mlxtran")
    )

  return(list(
    results = final_results[order(final_results$metric), ],
    best_model = final_results[order(final_results$metric), ][1, ],
    variability = variability,
    error = error
  ))
}
