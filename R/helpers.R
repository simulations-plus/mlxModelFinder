getAllFilters <- function(library, requiredFilters, custom_filters = NULL) {
  # If custom filters provided, return them directly
  if (library == "custom" && !is.null(custom_filters)) {
    return(custom_filters)
  }

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
  
  # For custom library, just return filters as-is
  if (library == "custom") {
    return(filters)
  }

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
                     obsIDToUse = NULL, linearization = FALSE,
                     param_distributions = NULL) {

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
  else
    path <- initial_func(project_path)
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

  if (!is.null(param_distributions)) {
    current_params <- lixoftConnectors::getIndividualParameterModel()$name
    dists_to_set <- param_distributions[names(param_distributions) %in% current_params]
    if (length(dists_to_set) > 0) {
      do.call(lixoftConnectors::setIndividualParameterDistribution, dists_to_set)
    }
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

calculateSearchSpace <- function(library, filters, error = FALSE, iiv = FALSE, param_mapping = NULL, fixed_iiv = NULL, requiredFilters = NULL) {
  if (library == "pk") {
    absorption_delays <- max(sum(filters$absorption != "transitCompartments"), 1) * length(filters$delay)
    direct_abs <- sum(filters$absorption == "transitCompartments")
    space <- (absorption_delays + direct_abs) * length(filters$distribution) * length(filters$elimination)
  } else if (library == "custom") {
    # For custom libraries, calculate space based on filter combinations
    if (length(filters) > 0) {
      space <- prod(sapply(filters, length))
    } else {
      space <- 1
    }
  } else {
    # Unknown library
    space <- 1
  }
  if (iiv) {
    if (library == "pk") {

      nb_params <- list(
        delay = data.frame(
          delay = c("noDelay", "lagTime"),
          level = c(0, 1)
        ),
        absorption = data.frame(
          absorption = c("firstOrder", "zeroOrder", "sigmoid", "transitCompartments"),
          level = c(1, 1, 2, 3)
        ),
        distribution = data.frame(
          distribution = c("1compartment", "2compartments", "3compartments"),
          level = c(1, 3, 5)
        ),
        elimination = data.frame(
          elimination = c("linear", "MichaelisMenten", "combined"),
          level = c(1, 2, 3)
        ),
        bioavailability = data.frame(
          bioavailability = c("false", "true"),
          level = c(0, 1)
        )
      )

      all_models <- expand.grid(filters, stringsAsFactors = FALSE)

      if ("transitCompartments" %in% all_models$absorption)
        all_models[all_models$absorption == "transitCompartments", ]$delay <- NA
      
      common_names <- intersect(names(nb_params), colnames(all_models))

      param_maps <- mapply(function(df, feat) {
        setNames(df$level, as.character(df[[feat]]))
      }, nb_params[common_names], common_names, SIMPLIFY = FALSE)

      param_values <- mapply(function(column_data, map) {
        unname(map[as.character(column_data)])
      }, all_models[common_names], param_maps, SIMPLIFY = FALSE)

      all_models$total_params <- rowSums(as.data.frame(param_values), na.rm = TRUE)

      if (!is.null(fixed_iiv) && !is.null(requiredFilters)) {
        # Count free params per model using actual parameter names
        space <- 0
        for (i in 1:nrow(all_models)) {
          model_filters <- as.list(all_models[i, , drop = FALSE])
          params <- get_pk_parameters(model_filters, requiredFilters)
          n_free <- length(setdiff(params, names(fixed_iiv)))
          space <- space + 2^n_free
        }
      } else {
        space <- sum(2^all_models$total_params)
      }
    } else if (library == "custom") {
      # For custom libraries, use param_mapping if available
      if (!is.null(param_mapping)) {
        # Calculate accurate IIV space based on param_mapping
        # Generate all structural model combinations
        all_models <- expand.grid(filters, stringsAsFactors = FALSE)
        
        # For each structural model, count parameters
        total_iiv_space <- 0
        for (i in 1:nrow(all_models)) {
          model_filters <- as.list(all_models[i, , drop = FALSE])
          
          # Get parameters for this model variant
          params <- getParamListFromMapping(model_filters, param_mapping)
          free_params <- if (!is.null(fixed_iiv)) setdiff(params, names(fixed_iiv)) else params
          n_free <- length(free_params)

          # Each free parameter can have IIV on or off: 2^n_free combinations
          total_iiv_space <- total_iiv_space + 2^n_free
        }
        
        space <- total_iiv_space
      } else {
        # No param_mapping - use conservative estimate
        # Assume average of 3 parameters per model, giving 2^3 = 8 IIV combinations
        space <- space * 8
      }
    }
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
          "Sampling Probability: %{y:.2f}<extra></extra>"
        )
      )
    }

    xmax <- max(df_cat$Iteration, na.rm = TRUE)
    p <- plotly::layout(
      p,
      title = list(text = title, x = 0.02),
      xaxis = list(title = "Iteration Number", rangemode = "tozero", range = c(0, xmax)),
      yaxis = list(title = "Sampling Probability", range = c(0, 1)),
      legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.2, yanchor = "top")
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

save_project <- function(mode, path_all, path_best, idx, improved) {

  project <- lixoftConnectors::getProjectSettings()$project

  if (identical(mode, "all")) {
    lixoftConnectors::saveProject(file.path(path_all, paste0("run_", idx, ".mlxtran")))
  } else if (identical(mode, "best") && isTRUE(improved)) {
    lixoftConnectors::saveProject(file.path(path_best, "best_run.mlxtran"))
  }

  lixoftConnectors::saveProject(project)
}

# Get parameter list for IIV based on filters and param_mapping
getParamListFromMapping <- function(filters, param_mapping) {
  if (is.null(param_mapping)) {
    return(character(0))
  }
  
  params <- character(0)
  
  # Add common parameters (always present)
  if (!is.null(param_mapping$common_params)) {
    params <- c(params, param_mapping$common_params)
  }
  
  # Add filter-specific parameters
  if (!is.null(param_mapping$filter_params)) {
    for (filter_name in names(filters)) {
      filter_value <- filters[[filter_name]]
      
      # Check if this filter category has parameter mappings
      if (filter_name %in% names(param_mapping$filter_params)) {
        category_mapping <- param_mapping$filter_params[[filter_name]]
        
        # Check if this specific filter value has parameters
        if (filter_value %in% names(category_mapping)) {
          params <- c(params, category_mapping[[filter_value]])
        }
      }
    }
  }
  
  return(unique(params))
}

# Generate IIV combinations, respecting fixed_iiv constraints.
# default_iiv: 1-row data frame with TRUE/FALSE columns for each param
# fixed_iiv: named logical list (or NULL) locking params to TRUE/FALSE
apply_fixed_iiv <- function(default_iiv, fixed_iiv = NULL) {
  iiv_cols <- names(default_iiv)

  # Determine which params are fixed vs free
  if (!is.null(fixed_iiv)) {
    fixed_params <- intersect(names(fixed_iiv), iiv_cols)
  } else {
    fixed_params <- character(0)
  }
  free_params <- setdiff(iiv_cols, fixed_params)

  # Generate combinations for free params only
  if (length(free_params) > 0) {
    iiv_combos <- expand.grid(
      setNames(rep(list(c(TRUE, FALSE)), length(free_params)), free_params)
    )
  } else {
    iiv_combos <- data.frame(row.names = 1)
  }

  # Add fixed params as constant columns
  for (p in fixed_params) {
    iiv_combos[[p]] <- fixed_iiv[[p]]
  }

  iiv_combos
}
