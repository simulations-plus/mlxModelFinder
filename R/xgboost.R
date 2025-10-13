make_applicable <- function(struct_features, param_df, param_cols, requiredFilters = list()) {

  params_out <- as.data.frame(param_df)

  # Helper to create a logical vector for a given rule across all rows.
  # It checks both the candidate's own features and the fixed requiredFilters.
  is_rule <- function(feature) {
    required_vals <- as.character(unlist(requiredFilters))
    apply(struct_features, 1, function(row) {
      all_features_for_row <- c(as.character(row), required_vals)
      feature %in% all_features_for_row
    })
  }

  # --- Create logical vectors for each rule ---
  has_no_delay        <- is_rule("noDelay")
  is_first_order      <- is_rule("firstOrder")
  is_zero_order       <- is_rule("zeroOrder")
  is_sigmoid          <- is_rule("sigmoid")
  is_transit          <- is_rule("transitCompartments")
  is_1cpt             <- is_rule("1compartment")
  is_2cpt             <- is_rule("2compartments")
  is_3cpt             <- is_rule("3compartments")
  is_linear_elim      <- is_rule("linear")
  is_mm_elim          <- is_rule("MichaelisMenten")
  has_bioavailability <- is_rule("true")
  is_clearance_param  <- is_rule("clearance")
  is_rate_param       <- is_rule("rate")

  has_any_absorption  <- is_first_order | is_zero_order | is_sigmoid | is_transit

  # --- Vectorized application of rules ---
  set_na_if <- function(df, condition, cols) {
    valid_cols <- intersect(cols, names(df))
    if (length(valid_cols) > 0 && any(condition)) {
      df[condition, valid_cols] <- NA
    }
    return(df)
  }

  # Delay
  params_out <- set_na_if(params_out, has_no_delay, "Tlag")

  # Absorption
  params_out <- set_na_if(params_out, is_first_order, c("Tk0", "Ktr", "Mtt", "Tk0s"))
  params_out <- set_na_if(params_out, is_zero_order, c("ka", "Ktr", "Mtt", "Tk0s"))
  params_out <- set_na_if(params_out, is_sigmoid, c("Tk0", "Ktr", "Mtt"))
  params_out <- set_na_if(params_out, is_transit, c("Tk0", "Tk0s", "Tlag"))
  params_out <- set_na_if(params_out, !has_any_absorption, c("ka", "Tk0", "Tk0s", "Ktr", "Mtt"))

  # Distribution
  params_out <- set_na_if(params_out, is_1cpt, c("V1", "V2", "V3", "Q2", "Q", "Q3", "k12", "k21", "k31", "k13"))
  params_out <- set_na_if(params_out, is_2cpt, c("V", "Q2", "V3", "Q3", "k13", "k31"))
  params_out <- set_na_if(params_out, is_3cpt, c("V", "Q"))

  # Elimination
  params_out <- set_na_if(params_out, is_linear_elim, c("Km", "Vm"))
  params_out <- set_na_if(params_out, is_mm_elim, c("Cl", "k"))

  # Bioavailability
  params_out <- set_na_if(params_out, !has_bioavailability, "F")

  # Parametrization
  params_out <- set_na_if(params_out, is_clearance_param, c("k12", "k21", "k13", "k31", "k"))
  params_out <- set_na_if(params_out, is_rate_param, c("Q", "Q2", "Q3", "Cl"))

  # --- Finalize: Ensure all required columns exist and are correctly ordered ---
  missing_cols <- setdiff(param_cols, names(params_out))
  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      params_out[[col]] <- NA
    }
  }

  return(params_out[, param_cols, drop = FALSE])
}


build_per_column_grid <- function(choices_list, max_combos, n_cols) {
  # expand.grid creates a data frame from all combinations of the supplied vectors
  full_grid <- expand.grid(choices_list, stringsAsFactors = FALSE)

  if (is.finite(max_combos) && nrow(full_grid) > max_combos) {
    # If the grid is too large, take a random sample
    sampled_indices <- sample.int(nrow(full_grid), size = max_combos)
    full_grid <- full_grid[sampled_indices, , drop = FALSE]
  }
  return(full_grid)
}

build_candidates <- function(probas_list, runs_df, target_cols, max_combos) {
  filters_df <- as.data.frame(probas_list$filters)

  choices_list <- lapply(target_cols, function(col) {
    # Values observed in training for this feature (including NA)
    seen_raw <- as.character(runs_df[[col]])
    seen_non_na <- unique(seen_raw[!is.na(seen_raw)])

    # If training never had a non-NA value, keep a single NA choice
    if (length(seen_non_na) == 0L) {
      allowed <- NA_character_
    } else {
      allowed <- seen_non_na
    }

    # What filters requested
    req_raw <- as.character(filters_df[[col]])
    allow_any <- any(is.na(req_raw))
    requested <- unique(req_raw[!is.na(req_raw)])

    # Apply restriction but never return 0-length
    vals <- if (length(requested) == 0L || allow_any) {
      allowed
    } else {
      inter <- intersect(requested, allowed)
      if (length(inter) == 0L) {
        # no overlap with training → use training-seen values
        allowed
      } else {
        inter
      }
    }

    # Final guard: if still empty (paranoia), use NA
    if (length(vals) == 0L) NA_character_ else vals
  })
  names(choices_list) <- target_cols

  # Sanity: ensure no 0-lengths
  stopifnot(all(sapply(choices_list, length) > 0L))

  build_per_column_grid(choices_list, max_combos, length(target_cols))
}

encode_for_model <- function(df, all_feature_cols, param_cols, levels_list) {
  for (col in param_cols) {
    if (col %in% names(df)) {
      val <- df[[col]]
      df[[col]] <- ifelse(is.na(val), "__NA__", ifelse(val, "TRUE", "FALSE"))
    }
  }

  for (col in all_feature_cols) {
    if (col %in% names(df)) {
      col_chr <- as.character(df[[col]])
      col_chr[is.na(col_chr)] <- "__NA__"
      df[[col]] <- factor(col_chr, levels = levels_list[[col]])
    }
  }

  return(df)
}

train_xgb_and_propose <- function(
    list_of_runs,
    probas_struct,
    probas_iiv,
    probas_error,
    N = 10L,
    metric_col = "metric",
    struct_max_combos = 10000L,
    iiv_max_combos = 10000L,
    max_total_combos = 200000L,
    requiredFilters = list(),
    modelsToTest = NULL
) {

  # --- 1. SETUP ---
  struct_cols <- colnames(probas_struct$filters)
  param_cols <- colnames(probas_iiv$filters)
  error_col <- colnames(probas_error$filters)[1]
  all_feature_cols <- c(struct_cols, error_col, param_cols)

  list_of_runs <- list_of_runs[list_of_runs$metric < 1e16, ]

  # --- 2. TRAIN XGBOOST MODEL ---
  feature_levels <- lapply(list_of_runs[all_feature_cols], function(col) {
    vals <- as.character(col)
    vals[is.na(vals)] <- "__NA__"
    if (is.logical(col)) return(union(vals, c("TRUE", "FALSE", "__NA__")))
    unique(vals)
  })

  runs_encoded <- encode_for_model(list_of_runs, all_feature_cols, param_cols, feature_levels)

  num_levels <- sapply(runs_encoded[all_feature_cols], function(col) length(levels(col)))
  constant_cols <- names(num_levels[num_levels < 2])

  if (length(constant_cols) > 0) {
    warning("The following columns have only one level and will be removed from the model: ",
            paste(constant_cols, collapse = ", "))
  }

  valid_feature_cols <- setdiff(all_feature_cols, constant_cols)

  if (length(valid_feature_cols) == 0) {
    stop("Error: No features with more than one level found. Cannot build a model.")
  }
  formula <- as.formula(paste("~", paste(valid_feature_cols, collapse = " + "), "- 1"))

  X_train <- model.matrix(formula, data = runs_encoded)

  y_train_raw <- list_of_runs[[metric_col]]
  y_rank <- -y_train_raw

  xgb_params <- list(
    objective = "rank:pairwise", max_depth = 6, eta = 0.05,
    subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 1
  )

  dtrain <- xgboost::xgb.DMatrix(data = X_train, label = y_rank)
  xgboost::setinfo(dtrain, "group", nrow(X_train))

  # cv <- xgboost::xgb.cv(
  #   params = xgb_params, data = dtrain, nrounds = 1000, nfold = 5,
  #   metrics = "rmse", early_stopping_rounds = 50, verbose = 0
  # )
  #
  # best_nrounds <- if (!is.null(cv$best_iteration)) cv$best_iteration else 200

  best_nrounds <- 300L
  model <- xgboost::xgb.train(params = xgb_params, data = dtrain, nrounds = best_nrounds, verbose = 0)

  # --- 3. GENERATE OR SELECT CANDIDATE MODELS ---

  if (is.null(modelsToTest)) {
    # --- A: Generate candidates automatically ---
    struct_candidates <- build_candidates(probas_struct, list_of_runs, struct_cols, struct_max_combos)
    iiv_candidates <- build_candidates(probas_iiv, list_of_runs, param_cols, iiv_max_combos)
    error_candidates <- build_candidates(probas_error, list_of_runs, error_col, Inf)

    indices <- expand.grid(
      struct = seq_len(nrow(struct_candidates)),
      iiv = seq_len(nrow(iiv_candidates)),
      error = seq_len(nrow(error_candidates))
    )

    if (is.finite(max_total_combos) && nrow(indices) > max_total_combos) {
      indices <- indices[sample.int(nrow(indices), max_total_combos), ]
    }

    struct_part <- struct_candidates[indices$struct, , drop = FALSE]
    iiv_part <- iiv_candidates[indices$iiv, , drop = FALSE]
    error_part <- error_candidates[indices$error, , drop = FALSE]

    fixed_iivs_df <- make_applicable(
      struct_features = struct_part,
      param_df = iiv_part,
      param_cols = param_cols,
      requiredFilters = requiredFilters
    )

    candidates <- cbind(struct_part, error_part, fixed_iivs_df)

  } else {
    # --- B: Use user-provided models ---
    message("Using user-provided 'modelsToTest' data frame. Skipping candidate generation.")

    user_struct_cols <- intersect(c(struct_cols, error_col), names(modelsToTest))
    struct_part_user <- modelsToTest[, user_struct_cols, drop = FALSE]

    user_param_cols <- intersect(param_cols, names(modelsToTest))
    param_part_user <- modelsToTest[, user_param_cols, drop = FALSE]

    fixed_iivs_df <- make_applicable(
      struct_features = struct_part_user,
      param_df = param_part_user,
      param_cols = param_cols,
      requiredFilters = requiredFilters
    )

    candidates <- cbind(struct_part_user, fixed_iivs_df)
  }

  candidates <- candidates[!duplicated(candidates), ]

  # --- 4. FILTER OUT ALREADY EVALUATED MODELS ---
  if (nrow(candidates) > 0) {
    # FIX: Create a unique key for each model configuration, ensuring column order is identical
    create_key <- function(df) {
      apply(df, 1, function(row) {
        row[is.na(row)] <- "_NA_"
        paste(row, collapse = "|")
      })
    }

    # Force both dataframes into the same canonical column order before creating keys
    runs_subset_ordered <- list_of_runs[, all_feature_cols, drop = FALSE]
    candidates_ordered <- candidates[, all_feature_cols, drop = FALSE]

    runs_keys <- create_key(runs_subset_ordered)
    candidate_keys <- create_key(candidates_ordered)

    is_new <- !(candidate_keys %in% runs_keys)
    candidates <- candidates[is_new, , drop = FALSE]
  }

  if (nrow(candidates) == 0) {
    warning("No new, unique candidate models to test after filtering already evaluated models.")
    return(data.frame())
  }

  # --- 5. PREDICT ON CANDIDATES & RETURN TOP N ---
  candidates_encoded <- encode_for_model(candidates, all_feature_cols, param_cols, feature_levels)
  X_cand <- model.matrix(formula, data = candidates_encoded)

  missing_cols <- setdiff(colnames(X_train), colnames(X_cand))
  if (length(missing_cols) > 0) {
    add_matrix <- matrix(0, nrow = nrow(X_cand), ncol = length(missing_cols), dimnames = list(NULL, missing_cols))
    X_cand <- cbind(X_cand, add_matrix)
  }
  X_cand <- X_cand[, colnames(X_train), drop = FALSE]

  candidates$pred_score <- predict(model, newdata = X_cand)
  ordered_candidates <- candidates[order(-candidates$pred_score), ]

  top_n_count <- min(N, nrow(ordered_candidates))
  top_results <- ordered_candidates[seq_len(top_n_count), , drop = FALSE]

  rownames(top_results) <- NULL
  return(top_results)
}

