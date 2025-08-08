runExhaustiveSearch <- function(library, filters, requiredFilters, settings, metric_func) {

  all_models <- expand.grid(filters, stringsAsFactors = FALSE)

  if ("transitCompartments" %in% all_models$delay)
    all_models[all_models$delay == "transitCompartments", ]$absorption <- NA

  all_models <- unique(all_models)

  if (settings$error) {
    errors <- data.frame(error = c("constant", "proportional", "combined1", "combined2"))
    all_models <- merge(all_models, errors)
  }

  res <- all_models
  res$metric <- rep(NA, nrow(res))

  cli_progress_bar(
    format = "Running models {cli::pb_bar} | {cli::pb_percent}",
    total = nrow(all_models),
    clear = FALSE
  )

  for (i in seq_len(nrow(all_models))) {
    model <- all_models[i, ]
    model <- getModelName("pk", filters = append(model[, !(names(model) %in% c("error"))], requiredFilters))

    error_model <- if (settings$error) as.data.frame(all_models$error)[i, ] else NULL

    runModel(model, error = error_model, obsIDToUse = settings$obsIDToUse)
    metric <- metric_func()
    res[i, ]$metric <- metric

    if (!is.null(settings$output))
      write.csv(res, file = settings$output, row.names = FALSE)

    cli_progress_update()
  }

  cli_progress_done()
  cli_alert_success("Done")
  return(res)
}
