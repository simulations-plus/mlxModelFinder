runExhaustiveSearch <- function(library, filters, requiredFilters, settings, metric_func, initial_func) {

  all_models <- expand.grid(filters, stringsAsFactors = FALSE)

  if ("transitCompartments" %in% all_models$absorption)
    all_models[all_models$absorption == "transitCompartments", ]$delay <- NA

  all_models <- unique(all_models)
  best_metric <- Inf

  if (settings$error) {
    errors <- data.frame(error = c("constant", "proportional", "combined1", "combined2"))
    all_models <- merge(all_models, errors)
  }

  res <- list()

  space <- calculateSearchSpace(library, filters, settings$error, settings$iiv)
  cli::cli_progress_bar(
    format = "Running models {cli::pb_bar} | {cli::pb_percent} | ETA: {cli::pb_eta}",
    total = space,
    clear = FALSE
  )

  for (i in seq_len(nrow(all_models))) {

    model <- all_models[i, ]
    model <- getModelName("pk", filters = append(model[, !(names(model) %in% c("error"))], requiredFilters))

    lixoftConnectors::setStructuralModel(model)

    error_model <- if (settings$error) as.data.frame(all_models$error)[i, ] else NULL

    if (settings$iiv) {
      default_iiv <- as.data.frame(as.list(lixoftConnectors::getIndividualParameterModel()$variability$id))
      iiv_combos <- expand.grid(lapply(default_iiv, function(x) c(TRUE, FALSE)))

      param_list <- c("ka", "Tk0", "Tk0s", "Tlag", "Ktr", "Mtt",
                      "V", "V1", "Cl", "Q", "Q2", "V2", "Q3", "V3",
                      "k", "k12", "k21", "k13", "k31",
                      "Vm", "Km", "F")

      param_cols <- as.data.frame(lapply(param_list, function(x) NA))
      names(param_cols) <- param_list
    } else {
      iiv_combos <- as.data.frame(as.list(lixoftConnectors::getIndividualParameterModel()$variability$id))
    }

    for (j in 1:nrow(iiv_combos)) {
      iiv_model <- iiv_combos[j, ]

      runModel(model, iiv = iiv_model, error = error_model, obsIDToUse = settings$obsIDToUse, linearization = settings$linearization, initial_func = initial_func)
      metric <- metric_func()

      if (metric < best_metric) {
        best_metric <- metric
        if (settings$save_mode == "best")
          lixoftConnectors::saveProject(
            file.path(settings$result_folder, "autoBuild", "best_run.mlxtran")
          )
      }

      if (!settings$iiv) {
        model_info <- all_models[i, ]
      } else {
        model_info <- cbind(all_models[i, ], param_cols)
        model_info[, names(iiv_model)] <- iiv_model
      }

      model_info$metric <- metric
      if (!is.null(settings$table_func))
        res[[length(res) + 1]] <- cbind(model_info, settings$table_func())
      else
        res[[length(res) + 1]] <- model_info


      if (!is.null(settings$output)) {
        results_df <- do.call(rbind, res)
        rownames(results_df) <- NULL

        write.csv(results_df, file = settings$output, row.names = FALSE)
      }

      cli::cli_progress_update()
    }
  }

  cli::cli_progress_done()
  cli::cli_alert_success("Done")

  results_df <- do.call(rbind, res)
  rownames(results_df) <- NULL

  return(results_df)
}
