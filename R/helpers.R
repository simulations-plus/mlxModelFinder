getAllFilters <- function(library, requiredFilters) {

  if (library == "pk" && requiredFilters$administration == "oralBolus")
    return(list(
      delay = c("noDelay", "lagTime"),
      absorption = c("zeroOrder", "firstOrder", "sigmoid", "transitCompartments"),
      distribution = c("1compartment", "2compartments", "3compartments"),
      elimination = c("linear", "MichaelisMenten", "combined"),
      bioavailability = c("true", "false")
    ))

  if (library == "pk" && requiredFilters$administration == "oral")
    return(list(
      delay = c("noDelay", "lagTime"),
      absorption = c("zeroOrder", "firstOrder", "sigmoid", "transitCompartments"),
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

calculateMetric <- function(RSEs, linearization) {
  method <- ifelse(linearization, "linearization", "importanceSampling")

  BICc <- lixoftConnectors::getEstimatedLogLikelihood()[[method]][4]
  if (is.null(BICc))
    return(1e12)

  if (!RSEs) return(BICc)

  method <- ifelse(linearization, "linearization", "stochasticApproximation")

  RSEs <- lixoftConnectors::getEstimatedStandardErrors()[[method]]$rse
  if (is.null(RSEs)) {
    RSEs <- rep(NaN, length(lixoftConnectors::getIndividualParameterModel()$name))
  }
  obs_name <- lixoftConnectors::getObservationInformation()$name
  n_tot <- length(lixoftConnectors::getObservationInformation()[[obs_name]][[obs_name]])

  return(BICc + length(RSEs[is.nan(RSEs) | RSEs > 50]) * log(n_tot))
}

calculateMetricBICc <- function() {
  BICc <- lixoftConnectors::getEstimatedLogLikelihood()$importanceSampling[4]
  if (is.null(BICc)) return(1e12)
  return(BICc)
}

runModel <- function(model, iiv = NULL, error = NULL, initial_func = NULL,
                     obsIDToUse = NULL, linearization = FALSE) {
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
                      conditionalDistributionSampling = !linearization,
                      conditionalModeEstimation = linearization,
                      standardErrorEstimation = TRUE,
                      logLikelihoodEstimation = TRUE)
  scenario$linearization <- linearization
  lixoftConnectors::setScenario(scenario)

  # Run scenario and handle an annoying error message that does not abort the task
  abort_msg <- "[ERROR] Error in function boost::math::cdf(const normal_distribution<double>&, double)"

  withCallingHandlers(
    lixoftConnectors::runScenario(),
    message = function(m) {
      if (grepl(abort_msg, conditionMessage(m), fixed = TRUE)) {
        lixoftConnectors:::.processRequest("monolix", "abort", list(), "asynchronous")
        scenario$tasks <- c(populationParameterEstimation = TRUE,
                            conditionalDistributionSampling = !linearization,
                            conditionalModeEstimation = linearization,
                            standardErrorEstimation = FALSE,
                            logLikelihoodEstimation = TRUE)
        lixoftConnectors::setScenario(scenario)
        lixoftConnectors::runScenario()
        invokeRestart("muffleMessage")
      }
    }
  )
}

calculateSearchSpace <- function(library, filters, error = FALSE, iiv = FALSE) {
  if (library == "pk") {
    absorption_delays <- max(sum(filters$absorption != "transitCompartments"), 1) * length(filters$delay)
    direct_abs <- sum(filters$absorption == "transitCompartments")
    space <- (absorption_delays + direct_abs) * length(filters$distribution) * length(filters$elimination)
  }
  if (iiv) {
    if (library == "pk") {
      nb_params <- list(
        delay = c(0, 1),
        absorption = c(1, 1, 2, 3),
        distribution = c(1, 3, 5),
        elimination = c(1, 2, 3),
        bioavailability = c(1, 0)
      )

      all_models <- expand.grid(filters, stringsAsFactors = FALSE)

      if ("transitCompartments" %in% all_models$absorption)
        all_models[all_models$absorption == "transitCompartments", ]$delay <- NA
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

plotProbabilities <- function(history_list, ncol = 2, panel_height = 360) {
  # history_list -> data.frame (base R only)
  df <- do.call(
    rbind,
    lapply(history_list, function(x) as.data.frame(x, stringsAsFactors = FALSE))
  )
  if (!nrow(df)) stop("history_list produced an empty data frame.")
  df <- df[!is.na(df$Meaning), , drop = FALSE]

  # ensure numeric, ordered lines
  df$Iteration   <- as.numeric(df$Iteration)
  df$Probability <- as.numeric(df$Probability)

  make_one_panel <- function(df_cat, title) {
    if (length(unique(df_cat$Meaning)) <= 1L) return(NULL)

    p <- plotly::plot_ly(height = panel_height)
    m_vals <- unique(df_cat$Meaning)

    for (m in m_vals) {
      sub <- df_cat[df_cat$Meaning == m, , drop = FALSE]
      sub <- sub[order(sub$Iteration), , drop = FALSE]

      p <- plotly::add_lines(
        p, data = sub,
        x = ~Iteration, y = ~Probability,
        name = m, showlegend = TRUE,
        hovertemplate = paste0(
          "Category: ", title, "<br>",
          "Meaning: ", m, "<br>",
          "Iteration: %{x}<br>",
          "Probability: %{y:.2f}<extra></extra>"
        )
      )
    }

    xmax <- max(df_cat$Iteration, na.rm = TRUE)
    p <- plotly::layout(
      p,
      title = list(text = title, x = 0.02),
      xaxis = list(title = "Iteration Number", rangemode = "tozero", range = c(0, xmax)),
      yaxis = list(title = "Probability", range = c(0, 1)),
      legend = list(orientation = "v", x = 1, xanchor = "right", y = 1, yanchor = "top")
    )
    plotly::config(p, displaylogo = FALSE, toImageButtonOptions = list(scale = 2))
  }

  cats <- unique(df$Category)
  panels <- lapply(cats, function(cat) {
    make_one_panel(df[df$Category == cat, , drop = FALSE], as.character(cat))
  })

  # drop NULLs (categories with a single series)
  keep <- !vapply(panels, is.null, logical(1))
  panels <- panels[keep]

  if (!length(panels)) {
    return(htmltools::browsable(
      htmltools::div("No panels to draw (all categories have a single series).")
    ))
  }

  # Build a responsive 2-column grid with htmltools
  grid_children <- lapply(panels, function(p) {
    htmltools::div(style = sprintf("height:%dpx; padding:4px;", panel_height), p)
  })

  container <- htmltools::div(
    style = sprintf("
      display:grid;
      grid-template-columns: repeat(%d, minmax(0,1fr));
      gap: 8px;
      align-items: start;
    ", ncol),
    grid_children
  )

  return(htmltools::browsable(container))
}

create_model_hash <- function(df_row) {
  df_row <- df_row[, sort(names(df_row)), drop = FALSE]
  paste(names(df_row), df_row, sep = "=", collapse = "&")
}
