updateModelContent <- function(model_content = NULL,
                               pk_lines = NULL,
                               equation_lines = NULL,
                               input_params = NULL,
                               output_params = NULL) {
  
  if (is.null(model_content)) {
    lines <- c(
      "[LONGITUDINAL]",
      "input = {}",
      "",
      "PK:",
      "",
      "EQUATION:",
      "",
      "OUTPUT:",
      "output = {}"
    )
  } else {
    lines <- model_content
  }
  
  # --- Helper function to add parameters to a line like 'input = {p1, p2}' ---
  add_params_to_line <- function(current_lines, line_prefix, params_to_add) {
    if (is.null(params_to_add)) {
      return(current_lines)
    }
    
    line_idx <- grep(paste0("^", line_prefix), current_lines)
    if (length(line_idx) == 0) {
      warning(paste("Warning: Could not find line starting with '", line_prefix, "'"))
      return(current_lines)
    }
    
    target_line <- current_lines[line_idx]
    
    # Extract existing parameters from inside the curly braces {}
    existing_params_str <- gsub(".*\\{(.*)\\}.*", "\\1", target_line)
    existing_params <- trimws(strsplit(existing_params_str, ",")[[1]])
    existing_params <- existing_params[existing_params != ""]
    
    # Combine with new parameters, ensuring no duplicates
    all_params <- unique(c(existing_params, params_to_add))
    
    # Reconstruct the line
    new_params_str <- paste(all_params, collapse = ", ")
    new_line <- paste0(line_prefix, " = {", new_params_str, "}")
    
    # Replace the old line with the new one
    current_lines[line_idx] <- new_line
    return(current_lines)
  }
  
  # --- Helper function to add lines to a specific section ---
  add_lines_to_section <- function(current_lines, section_header, lines_to_add) {
    if (is.null(lines_to_add)) {
      return(current_lines)
    }
    
    section_idx <- which(current_lines == section_header)
    if (length(section_idx) == 0) {
      warning(paste("Warning: Section '", section_header, "' not found."))
      return(current_lines)
    }
    
    # Find where the section ends (either next section header or end of file)
    all_headers <- c("PK:", "EQUATION:", "OUTPUT:", "[LONGITUDINAL]")
    next_headers_indices <- which(current_lines %in% all_headers)
    
    # Find the header that comes immediately after our target section
    next_section_idx_pos <- which(next_headers_indices > section_idx)
    
    if (length(next_section_idx_pos) > 0) {
      end_idx <- next_headers_indices[min(next_section_idx_pos)] -1
    } else {
      end_idx <- length(current_lines)
    }
    
    # Insert the new lines at the end of the section
    prefix <- current_lines[1:end_idx]
    suffix <- if (end_idx < length(current_lines)) current_lines[(end_idx + 1):length(current_lines)] else character(0)
    
    return(c(prefix, lines_to_add, suffix))
  }
  
  # Sequentially apply the modifications
  lines <- add_params_to_line(lines, "input", input_params)
  lines <- add_params_to_line(lines, "output", output_params)
  lines <- add_lines_to_section(lines, "PK:", pk_lines)
  lines <- add_lines_to_section(lines, "EQUATION:", equation_lines)

  return(lines)
}

getModelName <- function(library, filters = list()) {
  clean_filters <- .preprocess_filters(library, filters)
  
  is_library_pk <- library == "pk" &&
    (is.null(clean_filters$absorption) || clean_filters$absorption != "sigmoid") &&
    (is.null(clean_filters$elimination) || clean_filters$elimination != "combined")
  
  if (is_library_pk) {
    return(lixoftConnectors::getLibraryModelName(library, clean_filters))
  }
  
  if (library == "pk") {
    updates <- list(
      pk_lines = "compartment(cmt=1, volume=V, concentration=Cc)",
      input_params = "V",
      output_params = "Cc"
    )

    admin_updates <- .handle_administration(clean_filters)

    dist_updates <- .handle_distribution(clean_filters)
    elim_updates <- .handle_elimination(clean_filters)

    updates$pk_lines <- c(updates$pk_lines, admin_updates$pk_lines, dist_updates$pk_lines, elim_updates$pk_lines)
    updates$input_params <- c(updates$input_params, admin_updates$input_params, dist_updates$input_params, elim_updates$input_params)

    model_content <- updateModelContent(
      pk_lines = updates$pk_lines,
      equation_lines = updates$equation_lines,
      input_params = unique(updates$input_params),
      output_params = updates$output_params
    )
    
    path <- tempfile("model", fileext = ".txt")
    writeLines(model_content, path)
    return(path)
  }
  
  warning(paste("Model generation for library '", library, "' is not fully implemented."))
  return(NULL)
}


# --- Helper Functions ---

.preprocess_filters <- function(library, filters) {
  new_filters <- filters
  
  if (library == "pk") {
    if (!is.null(filters$delay) && filters$delay == "transitCompartments") {
      new_filters$absorption <- NULL
    }
    if (!is.null(filters$elimination) && filters$elimination == "MichaelisMenten" &&
        !is.null(filters$distribution) && filters$distribution == "1compartment") {
      new_filters$parametrization <- NULL
    }
  }
  return(new_filters)
}

.handle_administration <- function(filters) {
  # Extract filter values for easier access
  adm <- filters$administration
  delay <- filters$delay
  absorption <- filters$absorption
  bioavailability <- filters$bioavailability
  
  # Main switch for the administration route
  switch(adm,
         "bolus" = ,
         "infusion" = {
           # --- IV Administration ---
           switch(delay,
                  "noDelay" = list(pk_lines = "iv(cmt=1)", input_params = NULL),
                  "lagTime" = list(pk_lines = "iv(cmt=1, Tlag)", input_params = "Tlag"),
                  list(pk_lines = NULL, input_params = NULL) # Default
           )
         },
         
         "oral" = {
           # --- Oral Administration ---
           # Switch based on the DELAY model
           switch(delay,
                  "transitCompartments" = {
                    list(pk_lines = "absorption(cmt=1, Ktr, Mtt, ka)", input_params = c("Ktr", "Mtt", "ka"))
                  },
                  "noDelay" = ,
                  "lagTime" = {
                    tlag_pk <- if (delay == "lagTime") ", Tlag" else ""
                    tlag_param <- if (delay == "lagTime") "Tlag" else NULL
                    
                    switch(absorption,
                           "firstOrder" = ,
                           "zeroOrder" = {
                             param <- if (absorption == "firstOrder") "ka" else "Tk0"
                             pk_line <- sprintf("absorption(cmt=1, %s%s)", param, tlag_pk)
                             list(pk_lines = pk_line, input_params = c(param, tlag_param))
                           },
                           "sigmoid" = {
                             pk_lines <- c("compartment(cmt=0)", "transfer(from=0, to=1, kt=ka)", "elimination(cmt=0, k=0)")
                             pk_line_absorp <- sprintf("absorption(cmt=0, Tk0%s)", tlag_pk)
                             list(pk_lines = c(pk_lines, pk_line_absorp), input_params = c("Tk0", "ka", tlag_param))
                           },
                           list(pk_lines = NULL, input_params = NULL) # Default for absorption
                    )
                  },
                  list(pk_lines = NULL, input_params = NULL) # Default for delay
           )
         },
         
         "oralBolus" = {
           # --- Combined IV and Oral Administration ---
           iv_base <- "iv(adm=1, cmt=1)"
           
           # Handle bioavailability suffix and input
           bioavail_suffix <- if (!is.null(bioavailability) && bioavailability == "true") ", p=F" else ""
           bioavail_input <- if (!is.null(bioavailability) && bioavailability == "true") "F" else NULL
           
           # Switch based on the DELAY model
           switch(delay,
                  "transitCompartments" = {
                    pk_line <- sprintf("absorption(adm=2, cmt=1, Ktr, Mtt, ka%s)", bioavail_suffix)
                    list(pk_lines = c(iv_base, pk_line), input_params = c("Ktr", "Mtt", "ka", bioavail_input))
                  },
                  "noDelay" = ,
                  "lagTime" = {
                    # For noDelay and lagTime, we then check the absorption model
                    tlag_pk <- if (delay == "lagTime") ", Tlag" else ""
                    tlag_param <- if (delay == "lagTime") "Tlag" else NULL
                    
                    switch(absorption,
                           "firstOrder" = ,
                           "zeroOrder" = {
                             param <- if (absorption == "firstOrder") "ka" else "Tk0"
                             pk_line <- sprintf("absorption(adm=2, cmt=1, %s%s%s)", param, tlag_pk, bioavail_suffix)
                             list(pk_lines = c(iv_base, pk_line), input_params = c(param, tlag_param, bioavail_input))
                           },
                           "sigmoid" = {
                             pk_lines <- c("compartment(cmt=0)", "transfer(from=0, to=1, kt=ka)", "elimination(cmt=0, k=0)")
                             pk_line_absorp <- sprintf("absorption(adm=2, cmt=0, Tk0%s%s)", tlag_pk, bioavail_suffix)
                             list(pk_lines = c(iv_base, pk_lines, pk_line_absorp), input_params = c("Tk0", "ka", tlag_param, bioavail_input))
                           },
                           list(pk_lines = NULL, input_params = NULL) # Default for absorption
                    )
                  },
                  list(pk_lines = NULL, input_params = NULL) # Default for delay
           )
         },
         
         # Default case for the main switch
         list(pk_lines = NULL, input_params = NULL)
  )
}


.handle_distribution <- function(filters) {
  updates <- list(pk_lines = NULL, input_params = NULL)
  dist <- filters$distribution
  param <- filters$parametrization
  
  if (dist == "2compartments") {
    if (param == "clearance") updates <- list(pk_lines = "peripheral(k12=Q/V, k21=Q/V2)", input_params = c("Q", "V2"))
    else updates <- list(pk_lines = "peripheral(k12, k21)", input_params = c("k12", "k21"))
  } else if (dist == "3compartments") {
    if (param == "clearance") updates <- list(pk_lines = c("peripheral(k12=Q2/V, k21=Q2/V2)", "peripheral(k13=Q3/V, k31=Q3/V3)"), input_params = c("Q2", "V2", "Q3", "V3"))
    else updates <- list(pk_lines = c("peripheral(k12, k21)", "peripheral(k13, k31)"), input_params = c("k12", "k21", "k13", "k31"))
  }
  return(updates)
}


.handle_elimination <- function(filters) {
  updates <- list(pk_lines = NULL, input_params = NULL)
  elim <- filters$elimination
  param <- filters$parametrization
  
  if (elim %in% c("linear", "combined")) {
    if (param == "clearance") {
      updates$pk_lines <- c(updates$pk_lines, "elimination(cmt=1, k=Cl/V)")
      updates$input_params <- c(updates$input_params, "Cl")
    } else {
      updates$pk_lines <- c(updates$pk_lines, "elimination(cmt=1, k)")
      updates$input_params <- c(updates$input_params, "k")
    }
  }

  if (elim %in% c("MichaelisMenten", "combined")) {
    updates$pk_lines <- c(updates$pk_lines, "elimination(cmt=1, Vm, Km)")
    updates$input_params <- c(updates$input_params, "Vm", "Km")
  }
  return(updates)
}
