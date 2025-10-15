.is_scalar_logical <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
.is_scalar_number  <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
.is_scalar_string  <- function(x) is.character(x) && length(x) == 1L && !is.na(x)

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
  allowed_values <- c("pk")  # keep future values when supported
  if (!is.character(library) || length(library) != 1 || !(library %in% allowed_values)) {
    stop("Error: 'library' must be 'pk'.", call. = FALSE)
  }
  TRUE
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

validate_settings <- function(settings, algorithm = NULL) {
  if (!is.list(settings)) {
    stop("Error: 'settings' must be a list.", call. = FALSE)
  }

  # Allowed keys (sync with default_settings)
  allowed <- c(
    "N","iiv","error","seed","previous_runs","output","convergence_threshold",
    "initial_iiv_forced_iter","rho","obsIDToUse","plot","linearization","save_mode",
    "clip","keys","elitism","max_models","iiv_favor_true","xgboost","xgboost_model_nb",
    "local_search", "table_func"
  )

  unknown <- setdiff(names(settings), allowed)
  if (length(unknown)) {
    stop(sprintf("Error: Unknown 'settings' key(s): %s",
                 paste(sprintf("'%s'", unknown), collapse = ", ")), call. = FALSE)
  }

  # --- Generic type/range checks --------------------------------------------
  if (!is.null(settings$N) && !(.is_scalar_number(settings$N) && settings$N > 0))
    stop("Error: 'settings$N' must be a single positive number.", call. = FALSE)

  if (!is.null(settings$iiv) && !.is_scalar_logical(settings$iiv))
    stop("Error: 'settings$iiv' must be TRUE/FALSE.", call. = FALSE)

  if (!is.null(settings$error) && !.is_scalar_logical(settings$error))
    stop("Error: 'settings$error' must be TRUE/FALSE.", call. = FALSE)

  if (!is.null(settings$seed) && !.is_scalar_number(settings$seed))
    stop("Error: 'settings$seed' must be a single number.", call. = FALSE)

  if (!is.null(settings$previous_runs)) {
    if (!.is_scalar_string(settings$previous_runs))
      stop("Error: 'settings$previous_runs' must be a single string (path to .csv).", call. = FALSE)
    if (!file.exists(settings$previous_runs))
      stop("Error: file specified in 'settings$previous_runs' does not exist.", call. = FALSE)
    if (!grepl("\\.csv$", settings$previous_runs, ignore.case = TRUE))
      stop("Error: 'settings$previous_runs' must be a .csv file.", call. = FALSE)
  }

  if (!is.null(settings$output)) {
    if (!.is_scalar_string(settings$output))
      stop("Error: 'settings$output' must be a single string (path to .csv).", call. = FALSE)
    if (!grepl("\\.csv$", settings$output, ignore.case = TRUE))
      stop("Error: 'settings$output' must end with .csv.", call. = FALSE)
    out_dir <- dirname(normalizePath(settings$output, mustWork = FALSE))
    if (!dir.exists(out_dir))
      stop(sprintf("Error: directory for 'settings$output' does not exist: %s", out_dir), call. = FALSE)
  }

  if (!is.null(settings$convergence_threshold) &&
      !(.is_scalar_number(settings$convergence_threshold) &&
        settings$convergence_threshold > 0 && settings$convergence_threshold <= 1))
    stop("Error: 'settings$convergence_threshold' must be in (0, 1].", call. = FALSE)

  if (!is.null(settings$initial_iiv_forced_iter) &&
      !(.is_scalar_number(settings$initial_iiv_forced_iter) && settings$initial_iiv_forced_iter >= 0))
    stop("Error: 'settings$initial_iiv_forced_iter' must be a non-negative number.", call. = FALSE)

  if (!is.null(settings$rho) &&
      !(.is_scalar_number(settings$rho) && settings$rho > 0 && settings$rho <= 1))
    stop("Error: 'settings$rho' must be in (0, 1].", call. = FALSE)

  if (!is.null(settings$obsIDToUse) && !.is_scalar_string(settings$obsIDToUse))
    stop("Error: 'settings$obsIDToUse' must be a single string.", call. = FALSE)

  if (!is.null(settings$plot) && !.is_scalar_logical(settings$plot))
    stop("Error: 'settings$plot' must be TRUE/FALSE.", call. = FALSE)

  if (!is.null(settings$linearization) && !.is_scalar_logical(settings$linearization))
    stop("Error: 'settings$linearization' must be TRUE/FALSE.", call. = FALSE)

  if (!is.null(settings$save_mode)) {
    allowed_sm <- c("none","best","all")
    if (!.is_scalar_string(settings$save_mode) || !(settings$save_mode %in% allowed_sm))
      stop("Error: 'settings$save_mode' must be one of 'none', 'best', 'all'.", call. = FALSE)
  }

  if (!is.null(settings$clip) && !.is_scalar_logical(settings$clip))
    stop("Error: 'settings$clip' must be TRUE/FALSE.", call. = FALSE)

  if (!is.null(settings$keys) && !(is.character(settings$keys) && length(settings$keys) >= 1L))
    stop("Error: 'settings$keys' must be a non-empty character vector.", call. = FALSE)

  if (!is.null(settings$elitism) && !.is_scalar_logical(settings$elitism))
    stop("Error: 'settings$elitism' must be TRUE/FALSE.", call. = FALSE)

  if (!is.null(settings$max_models) && !(.is_scalar_number(settings$max_models) && settings$max_models >= 1))
    stop("Error: 'settings$max_models' must be a positive number.", call. = FALSE)

  if (!is.null(settings$iiv_favor_true) && !.is_scalar_logical(settings$iiv_favor_true))
    stop("Error: 'settings$iiv_favor_true' must be TRUE/FALSE.", call. = FALSE)

  if (!is.null(settings$xgboost) && !.is_scalar_logical(settings$xgboost))
    stop("Error: 'settings$xgboost' must be TRUE/FALSE.", call. = FALSE)

  if (!is.null(settings$xgboost_model_nb) &&
      !(.is_scalar_number(settings$xgboost_model_nb) && settings$xgboost_model_nb >= 1))
    stop("Error: 'settings$xgboost_model_nb' must be a positive number.", call. = FALSE)

  if (!is.null(settings$local_search) && !.is_scalar_logical(settings$local_search))
    stop("Error: 'settings$local_search' must be TRUE/FALSE.", call. = FALSE)

  # --- Algorithm-specific enforcement ---------------------------------------
  if (!is.null(algorithm)) {
    algo <- tolower(algorithm)

    # Settings only meaningful for stochastic search (ACO/Tournament)
    stochastic_only <- c("convergence_threshold","rho","initial_iiv_forced_iter","clip",
                         "xgboost","xgboost_model_nb","elitism","local_search")

    if (algo %in% c("decision_tree","exhaustive_search")) {
      present <- intersect(names(settings), stochastic_only)
      if (length(present)) {
        stop(sprintf(
          "Error: setting(s) %s are only valid for 'ACO' or 'tournament'.",
          paste(sprintf("'%s'", present), collapse = ", ")
        ), call. = FALSE)
      }
    }

    # 'keys' used by decision tree only
    if (!is.null(settings$keys) && !(algo %in% c("decision_tree"))) {
      stop("Error: 'settings$keys' is only valid with algorithm 'decision_tree'.", call. = FALSE)
    }

    # ACO-specific
    if (algo == "aco") {
      # ok: convergence_threshold, rho, initial_iiv_forced_iter, clip, xgboost*, elitism, local_search
      # nothing special to reject here
    }

    # Tournament-specific
    if (algo == "tournament") {
      # ok: rho, initial_iiv_forced_iter, clip, xgboost*
      # reject ACO-only fields if you consider them ACO-only (e.g., 'convergence_threshold', 'elitism', 'local_search')
      aco_only <- intersect(names(settings), c("convergence_threshold","elitism","local_search"))
      if (length(aco_only)) {
        stop(sprintf(
          "Error: setting(s) %s are only valid for 'ACO'.",
          paste(sprintf("'%s'", aco_only), collapse = ", ")
        ), call. = FALSE)
      }
    }
  }

  TRUE
}

validate_metric_func_no_args <- function(func) {
  if (!is.function(func)) {
    stop("Error: 'metric_func' must be a function.", call. = FALSE)
  }
  if (length(formals(func)) != 0) {
    stop("Error: 'metric_func' must take no arguments and read results via lixoftConnectors.", call. = FALSE)
  }
  TRUE
}

validate_initial_func_one_arg <- function(func) {
  if (!is.function(func)) {
    stop("Error: 'initial_func' must be a function.", call. = FALSE)
  }
  if (length(formals(func)) != 1) {
    stop("Error: 'initial_func' must accept exactly one argument: the project path.", call. = FALSE)
  }
  TRUE
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
      sprintf("%s %s %s not installed, so the %s will not be used.",
              pkg_word, paste(sprintf("'%s'", missing), collapse = " and "), verb, what),
      call. = FALSE
    )
    return(FALSE)
  }
  TRUE
}
