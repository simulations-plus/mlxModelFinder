getAllFilters <- function(library, requiredFilters) {

  if (library == "pk" && requiredFilters$administration == "oralBolus")
    return(list(
      delay = c("noDelay", "lagTime", "transitCompartments"),
      absorption = c("zeroOrder", "firstOrder", "sigmoid"),
      distribution = c("1compartment", "2compartments", "3compartments"),
      elimination = c("linear", "MichaelisMenten", "combined"),
      bioavailability = c("true", "false")
    ))

  if (library == "pk" && requiredFilters$administration == "oral")
    return(list(
      delay = c("noDelay", "lagTime", "transitCompartments"),
      absorption = c("zeroOrder", "firstOrder", "sigmoid"),
      distribution = c("1compartment", "2compartments", "3compartments"),
      elimination = c("linear", "MichaelisMenten", "combined")
    ))

  if (library == "pk")
    return(list(
      delay = c("noDelay"),
      distribution = c("1compartment", "2compartments", "3compartments"),
      elimination = c("linear", "MichaelisMenten", "combined")
    ))
}

updateFilters <- function(library, filters, requiredFilters) {

  all_filters <- getAllFilters(library, requiredFilters)
  updated_filters <- list()
  for (key in names(all_filters)) {
    if (key %in% names(filters))
      updated_filters <- append(updated_filters, filters[key])
    else
      updated_filters <- append(updated_filters, all_filters[key])
  }

  return(updated_filters)
}

getFilterMatrix <- function(filters) {
  return(sapply(filters, '[', seq(max(sapply(filters, length)))))
}

initializeWeightMatrix <- function(filter_matrix) {
  matrix <- apply(filter_matrix, 1:2, function(x)
    as.numeric(!is.na(x)))

  return(matrix)
}

calculateMetricRSEs <- function() {
  BICc <- lixoftConnectors::getEstimatedLogLikelihood()$importanceSampling[4]
  RSEs <- lixoftConnectors::getEstimatedStandardErrors()$stochasticApproximation$rse
  if (is.null(RSEs)) {
    RSEs <- rep(NaN, length(lixoftConnectors::getIndividualParameterModel()$name))
  }
  obs_name <- lixoftConnectors::getObservationInformation()$name
  n_tot <- length(lixoftConnectors::getObservationInformation()[[obs_name]][[obs_name]])
  if (is.null(BICc)) return(1e12)
  return(BICc + length(RSEs[is.nan(RSEs) | RSEs > 50]) * log(n_tot))
}

calculateMetricBICc <- function() {
  BICc <- lixoftConnectors::getEstimatedLogLikelihood()$importanceSampling[4]
  if (is.null(BICc)) return(1e12)
  return(BICc)
}

renameParams <- function(input, library) {
  if (library == "pk") {
    map <- list(V = "V1", Q = "Q2")
    input <- ifelse(input %in% names(map), unlist(map[input]), input)
  }
  return(input)
}

runModel <- function(model, iiv = NULL, error = NULL, initial_func = NULL, obsIDToUse = NULL) {
  error_model <- lixoftConnectors::getContinuousObservationModel()$errorModel
  lixoftConnectors::setStructuralModel(model)

  if (!is.null(obsIDToUse)) {
    model_output <- lixoftConnectors::getMapping()$mapping[[1]]$modelOutput
    lixoftConnectors::setMapping(list(list(obsId = obsIDToUse, modelOutput = model_output)))
    if (obsIDToUse %in% names(error_model))
      lixoftConnectors::setErrorModel(error_model)
  }

  project_path <- gsub("\\\\", "/", tempfile(pattern = "project", fileext = ".mlxtran"))
  lixoftConnectors::saveProject(project_path)
  if (is.null(initial_func))
    path <- initialValues(project_path)
  lixoftConnectors::loadProject(path)

  if (!is.null(iiv)) {
    param_names <- lixoftConnectors::getIndividualParameterModel()$name
    if ("V1" %in% names(iiv) && "V" %in% param_names) {
      names(iiv) <- ifelse(names(iiv) == "V1", "V", names(iiv))
    }
    if ("Q2" %in% names(iiv) && "Q" %in% param_names) {
      names(iiv) <- ifelse(names(iiv) == "Q2", "Q", names(iiv))
    }
    iiv <- unlist(iiv[, !is.na(iiv[1, ])])
    lixoftConnectors::setIndividualParameterModel(list(variability = list(id = iiv)))
  }

  if (!is.null(error)) {
    error <- as.list(error)
    names(error) <- lixoftConnectors::getObservationInformation()$name[[1]]
    lixoftConnectors::setErrorModel(error)
  }

  scenario <- lixoftConnectors::getScenario()
  scenario$tasks <- c(populationParameterEstimation = TRUE,
                      conditionalDistributionSampling = TRUE,
                      conditionalModeEstimation = FALSE,
                      standardErrorEstimation = TRUE,
                      logLikelihoodEstimation = TRUE)
  scenario$linearization <- FALSE
  lixoftConnectors::setScenario(scenario)
  lixoftConnectors::runScenario()
}

calculateSearchSpace <- function(library, filters, error = FALSE, iiv = FALSE) {
  if (library == "pk") {
    absorption_delays <- sum(filters$delay %in% c("noDelay", "lagTime")) * max(length(filters$absorption), 1)
    direct_delays <- sum(filters$delay == "transitCompartments")
    space <- (absorption_delays + direct_delays) * length(filters$distribution) * length(filters$elimination)
  }
  if (iiv) {
    if (library == "pk") {
      nb_params <- list(
        delay = c(0, 1, 3),
        absorption = c(1, 1, 2),
        distribution = c(1, 3, 5),
        elimination = c(1, 2, 3),
        bioavailability = c(1, 0)
      )

      all_models <- expand.grid(filters, stringsAsFactors = FALSE)

      if ("transitCompartments" %in% all_models$delay)
        all_models[all_models$delay == "transitCompartments", ]$absorption <- NA
    }

    common_names <- intersect(names(nb_params), names(filters))
    param_maps <- mapply(setNames, nb_params[common_names], filters[common_names], SIMPLIFY = FALSE)

    param_values <- mapply(function(column_data, map) {
      map[column_data]
    }, all_models, param_maps)
    all_models$total_params <- rowSums(param_values, na.rm = TRUE)
    space <- sum(2^all_models$total_params)
  }
  if (error) space <- space * 4
  return(space)
}

tidyProbabilities <- function(probas, iteration) {
  col_names <- colnames(probas$probas)
  list_of_cols <- vector("list", length(col_names))
  names(list_of_cols) <- col_names

  for (col_name in col_names) {
    list_of_cols[[col_name]] <- data.frame(
      Iteration = iteration,
      Category = col_name,
      Probability = probas$probas[, col_name],
      Meaning = probas$filters[, col_name]
    )
  }

  current_tidy_data <- do.call(rbind, list_of_cols)
  rownames(current_tidy_data) <- NULL

  return(current_tidy_data)
}

plotProbabilities <- function(history_list, iteration) {

  plot_data <- do.call(rbind, history_list)
  plot_data <- subset(plot_data, !is.na(Meaning))

  out <- by(data = plot_data, INDICES = plot_data$Category, FUN = function(p) {
    p <- droplevels(p)
    label_data <- subset(p, Iteration == iteration)
    p <- ggplot(
      data = p,
      aes(x = Iteration, y = Probability, group = Meaning)
    ) +
      facet_wrap(~Category, ncol = 1) +
      geom_line(color = "grey60") +
      geom_point(data = label_data, aes(color = Meaning), size = 3) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_x_continuous(
        limits = c(0, iteration),
        breaks = seq(0, iteration, 1),
        minor_breaks = NULL
      ) +
      labs(
        x = "Iteration Number",
        y = "Probability"
      ) +
      theme(legend.position = "right")
  })
  grid <- do.call(grid.arrange, c(out, ncol=2))
  return(grid)
}

create_model_hash <- function(df_row) {
  df_row <- df_row[, sort(names(df_row)), drop = FALSE]
  paste(names(df_row), df_row, sep = "=", collapse = "&")
}
