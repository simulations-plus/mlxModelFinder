get_pk_parameters <- function(filters, requiredFilters) {
  params <- character(0)

  # --- Administration and Bioavailability ---
  if (!is.null(filters$bioavailability) &&
      "true" %in% filters$bioavailability &&
      requiredFilters$administration == "oralBolus") {
    params <- c(params, "F")
  }

  # --- Delay ---
  if (!is.null(filters$delay) && "lagTime" %in% filters$delay) {
    params <- c(params, "Tlag")
  }

  # --- Absorption ---
  if (!is.null(filters$absorption)) {
    for (absorption_model in filters$absorption) {
      params <- switch(
        absorption_model,
        "zeroOrder" = c(params, "Tk0"),
        "firstOrder" = c(params, "ka"),
        "sigmoid" = c(params, "ka", "Tk0s"),
        "transitCompartments" = c(params, "ka", "Mtt", "Ktr"),
        params
      )
    }
  }

  # --- Distribution and Elimination ---
  is_clearance_param <- requiredFilters$parametrization == "clearance"

  # --- Elimination ---
  if (!is.null(filters$elimination)) {
    for (elimination_model in filters$elimination) {
      if (elimination_model %in% c("linear", "combined")) {
        if (is_clearance_param) {
          params <- c(params, "Cl")
        } else {
          params <- c(params, "k")
        }
      }
      if (elimination_model %in% c("MichaelisMenten", "combined")) {
        params <- c(params, "Vm", "Km")
      }
    }
  }

  # --- Distribution ---
  if (!is.null(filters$distribution)) {
    for (dist_model in filters$distribution) {
      if (dist_model == "1compartment") {
        params <- c(params, "V")
      } else if (dist_model == "2compartments") {
        if (is_clearance_param) {
          params <- c(params, "Q", "V1", "V2")
        } else {
          params <- c(params, "k12", "k21", "V")
        }
      } else if (dist_model == "3compartments") {
        if (is_clearance_param) {
          params <- c(params, "Q2", "Q3", "V1", "V2", "V3")
        } else {
          params <- c(params, "k12", "k21", "k13", "k31", "V")
        }
      }
    }
  }

  return(unique(params))
}
probabilities <- function(x, weights_override = NULL) {
  filters <- sapply(x, `[`, seq(max(sapply(x, length))))
  if (is.vector(filters)) filters <- t(as.matrix(filters))

  if (is.null(weights_override)) {
    weights <- apply(filters, 1:2, function(x) as.numeric(!is.na(x)))
  } else {
    wmat <- sapply(weights_override, `[`, seq(nrow(filters)))
    wmat[is.na(filters)] <- NA
    weights <- wmat
  }

  probas <- apply(weights, 2, function(col) col / sum(col, na.rm = TRUE))
  if (is.vector(probas)) probas <- t(as.matrix(probas))

  list(probas = probas, filters = filters, weights = weights)
}

iiv_probabilities <- function(library, filters, requiredFilters, favor_true = FALSE) {
  if (library == "pk") {
    # param_list <- c("ka", "Tk0", "Tk0s", "Tlag", "Ktr", "Mtt",
    #                 "V", "V1", "Cl", "Q", "Q2", "V2", "Q3", "V3",
    #                 "k", "k12", "k21", "k13", "k31",
    #                 "Vm", "Km", "F")
    param_list <- get_pk_parameters(filters, requiredFilters)
    param_map <- lapply(param_list, function(x) c(TRUE, FALSE))
    names(param_map) <- param_list
  } else if (library == "custom") {
    # For custom libraries, IIV parameters are determined dynamically
    # Return empty probabilities - ACO will handle this differently
    return(NULL)
  } else {
    stop("Unsupported library", call. = FALSE)
  }

  if (!favor_true) {
    return(probabilities(param_map))
  }

  strong_true  <- c("Cl", "k", "V")  # TRUE:FALSE = 9:1
  medium_true  <- c("ka", "Tk0", "Tlag", "Tk0s")          # TRUE:FALSE = 3:1

  weight_map <- lapply(names(param_map), function(p) {
    if (p %in% strong_true) {
      c(9, 1)
    } else if (p %in% medium_true) {
      c(3, 1)
    } else {
      c(1, 1)
    }
  })
  names(weight_map) <- names(param_map)

  probabilities(param_map, weights_override = weight_map)
}


updateProbabilities <- function(obj, runs, rho = 0.4, clip = FALSE,
                                min_val = 0.05, max_val = 0.95, tol = 1e-3) {
  # decay previous weights
  obj$weights <- obj$weights * (1 - rho)

  keys <- colnames(obj$weights)
  for (key in keys) {
    for (i in 1:nrow(obj$weights)) {
      obj$weights[i, key] <- obj$weights[i, key] +
        sum(1 / runs[runs[, key] == obj$filters[i, key], "rank"], na.rm = TRUE)
    }
  }

  # recompute probabilities
  obj$probas <- apply(obj$weights, 2, function(x) {
    p <- x / sum(x, na.rm = TRUE)
    if (clip) {
      clip_probabilities(p, min_val = min_val, max_val = max_val, tol = tol)
    } else {
      p
    }
  })

  if (is.vector(obj$probas)) obj$probas <- t(as.matrix(obj$probas))

  return(obj)
}


clip_probabilities <- function(p, min_val = 0.05, max_val = 0.95, tol = 1e-3) {
  if (any(p == 1)) return(ifelse(p == 1, 1, 0))

  zero_idx <- which(p == 0)
  pos_idx  <- which(p > 0 & p < 1)

  if (length(pos_idx) == 0) return(p)

  repeat {
    p[pos_idx] <- pmin(pmax(p[pos_idx], min_val), max_val)

    scale_factor <- sum(p[pos_idx])
    if (scale_factor > 0) {
      p[pos_idx] <- p[pos_idx] / scale_factor
    }

    if (all(p[pos_idx] >= min_val - tol & p[pos_idx] <= max_val + tol)) break
  }

  p
}
