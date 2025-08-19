validate_project <- function(project) {
  if (!is.character(project) || length(project) != 1) {
    stop("Error: 'project' must be a single string.", call. = FALSE)
  }
  if (!file.exists(project)) {
    stop("Error: The file specified by 'project' does not exist.", call. = FALSE)
  }
  if (!grepl("\\.mlxtran$", project)) {
    stop("Error: The 'project' file must have the .mlxtran extension.", call. = FALSE)
  }
  return(TRUE)
}

validate_library <- function(library) {
  allowed_values <- c("pk", "tmdd", "pkpd")
  if (!is.character(library) || length(library) != 1 || !(library %in% allowed_values)) {
    stop("Error: 'library' must be one of 'pk', 'tmdd', or 'pkpd'.", call. = FALSE)
  }
  return(TRUE)
}

validate_required_filters <- function(requiredFilters) {
  if (!is.list(requiredFilters) || !all(c("administration", "parametrization") %in% names(requiredFilters))) {
    stop("Error: 'requiredFilters' must be a list with headers 'administration' and 'parametrization'.", call. = FALSE)
  }
  administration_values <- c("bolus", "infusion", "oral", "oralBolus")
  parametrization_values <- c("rate", "clearance")
  if (!(requiredFilters$administration %in% administration_values)) {
    stop("Error: 'administration' must be one of 'bolus', 'infusion', 'oral', or 'oralBolus'.", call. = FALSE)
  }
  if (!(requiredFilters$parametrization %in% parametrization_values)) {
    stop("Error: 'parametrization' must be 'rate' or 'clearance'.", call. = FALSE)
  }
  return(TRUE)
}

validate_filters <- function(filters, requiredFilters) {
  allowed_names <- c("delay", "absorption", "distribution", "elimination", "bioavailability")
  if (!is.list(filters) || !all(names(filters) %in% allowed_names)) {
    stop("Error: 'filters' must be a list with keys from: 'delay', 'absorption', 'distribution', 'elimination', 'bioavailability'.", call. = FALSE)
  }

  delay_values <- c("noDelay", "lagTime")
  absorption_values <- c("zeroOrder", "firstOrder", "sigmoid", "transitCompartments")
  distribution_values <- c("1compartment", "2compartments", "3compartments")
  elimination_values <- c("linear", "MichaelisMenten", "combined")

  if (!is.null(filters$delay) && !all(filters$delay %in% delay_values)) {
    stop("Error: 'delay' must contain only values 'noDelay', or 'lagTime'.", call. = FALSE)
  }
  if (!is.null(filters$absorption) && !all(filters$absorption %in% absorption_values)) {
    stop("Error: 'absorption' must contain only values 'zeroOrder', 'firstOrder', 'sigmoid', or 'transitCompartments'.", call. = FALSE)
  }
  if (!is.null(filters$distribution) && !all(filters$distribution %in% distribution_values)) {
    stop("Error: 'distribution' must contain only values '1compartment', '2compartments', or '3compartments'.", call. = FALSE)
  }
  if (!is.null(filters$elimination) && !all(filters$elimination %in% elimination_values)) {
    stop("Error: 'elimination' must contain only values 'linear', 'MichaelisMenten', or 'combined'.", call. = FALSE)
  }

  if (!is.null(filters$bioavailability)) {
    if (requiredFilters$administration != "oralBolus") {
      stop("Error: 'bioavailability' filter is only applicable when administration is 'oralBolus'.", call. = FALSE)
    }
    if (!is.character(filters$bioavailability) || !all(filters$bioavailability %in% c("true", "false"))) {
      stop("Error: 'bioavailability' must contain only string values 'true' or 'false'.", call. = FALSE)
    }
  }

  return(TRUE)
}
validate_algorithm <- function(algorithm) {
  allowed_algorithms <- c("ACO", "tournament", "exhaustive_search", "decision_tree")
  if (!is.character(algorithm) || length(algorithm) != 1 || !(algorithm %in% allowed_algorithms)) {
    stop("Error: 'algorithm' must be one of 'ACO', 'tournament', 'exhaustive_search', or 'decision_tree'.", call. = FALSE)
  }
  return(TRUE)
}

validate_settings <- function(settings) {
  if (!is.list(settings)) {
    stop("Error: 'settings' must be a list.", call. = FALSE)
  }

  if (!is.null(settings$N) && (!is.numeric(settings$N) || length(settings$N) != 1 || settings$N <= 0)) {
    stop("Error: 'settings$N' must be a single positive number.", call. = FALSE)
  }
  if (!is.null(settings$iiv) && (!is.logical(settings$iiv) || length(settings$iiv) != 1)) {
    stop("Error: 'settings$iiv' must be a single boolean value (TRUE/FALSE).", call. = FALSE)
  }
  if (!is.null(settings$error) && (!is.logical(settings$error) || length(settings$error) != 1)) {
    stop("Error: 'settings$error' must be a single boolean value (TRUE/FALSE).", call. = FALSE)
  }
  if (!is.null(settings$seed) && (!is.numeric(settings$seed) || length(settings$seed) != 1)) {
    stop("Error: 'settings$seed' must be a single number.", call. = FALSE)
  }
  if (!is.null(settings$previous_runs) && (!is.character(settings$previous_runs) || length(settings$previous_runs) != 1)) {
    stop("Error: 'settings$previous_runs' must be a single string (file path).", call. = FALSE)
  }
  if (!is.null(settings$output) && (!is.character(settings$output) || length(settings$output) != 1)) {
    stop("Error: 'settings$output' must be a single string (file path).", call. = FALSE)
  }
  if (!is.null(settings$convergence_threshold) && (!is.numeric(settings$convergence_threshold) || length(settings$convergence_threshold) != 1 || settings$convergence_threshold <= 0 || settings$convergence_threshold > 1)) {
    stop("Error: 'settings$convergence_threshold' must be a single number between 0 and 1.", call. = FALSE)
  }
  if (!is.null(settings$rho) && (!is.numeric(settings$rho) || length(settings$rho) != 1 || settings$rho <= 0 || settings$rho > 1)) {
    stop("Error: 'settings$rho' must be a single number between 0 and 1.", call. = FALSE)
  }
  if (!is.null(settings$initial_iiv_forced_iter) && (!is.numeric(settings$initial_iiv_forced_iter) || length(settings$initial_iiv_forced_iter) != 1 || settings$initial_iiv_forced_iter < 0)) {
    stop("Error: 'settings$initial_iiv_forced_iter' must be a single non-negative number.", call. = FALSE)
  }

  return(TRUE)
}

validate_no_argument_function <- function(func, arg_name) {
  if (!is.function(func)) {
    stop(paste0("Error: The argument '", arg_name, "' must be a function."), call. = FALSE)
  }
  if (length(formals(func)) != 0) {
    if (arg_name == "metric_func") {
      stop("Error: The 'metric_func' function must take no arguments. It must retrieve the results using lixoftConnectors functions.", call. = FALSE)
    } else {
      stop(paste0("Error: The '", arg_name, "' function must take no arguments."), call. = FALSE)
    }
  }
  return(TRUE)
}

validate_obsidtouse <- function(settings) {
  obs_ids <- lixoftConnectors::getData()$observationNames
  if (length(obs_ids) > 1) {
    if (is.null(settings$obsIDToUse)) {
      stop("settings$obsIDToUse must be provided when multiple observation types are present.", call. = FALSE)
    }

    if (!is.character(settings$obsIDToUse) || length(settings$obsIDToUse) != 1) {
      stop("settings$obsIDToUse must be a single string.", call. = FALSE)
    }

    if (!(settings$obsIDToUse %in% obs_ids)) {
      stop(sprintf(
        "Invalid obsIDToUse: '%s'. It must match one of the available observation types: %s",
        settings$obsIDToUse,
        paste(shQuote(obs_ids), collapse = ", ")
      ), call. = FALSE)
    }
  }

  return(TRUE)
}

.check_pkgs <- function(pkgs, what) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    pkg_word <- if (length(missing) == 1L) "Package" else "Packages"
    verb     <- if (length(missing) == 1L) "is" else "are"
    warning(
      sprintf("%s %s %s not installed, so the %s will not be plotted.",
              pkg_word, paste(sprintf("'%s'", missing), collapse = " and "), verb, what),
      call. = FALSE
    )
    return(FALSE)
  }
  TRUE
}
