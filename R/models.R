models <- function(p, n_ants, p_err = NULL, p_iiv = NULL) {
  x <- sapply(
    colnames(p$probas),
    function(x) sample(p$filters[, x], n_ants, TRUE, p$probas[, x])
  )

  x <- list(struct = x)
  if (!is.null(p_iiv))
    x$iiv <- sapply(
      colnames(p_iiv$probas),
      function(x) sample(p_iiv$filters[, x], n_ants, TRUE, p_iiv$probas[, x])
    )

  if (!is.null(p_err))
    x$error <- sapply(
      colnames(p_err$probas),
      function(x) sample(p_err$filters[, x], n_ants, TRUE, p_err$probas[, x])
    )

  if (n_ants == 1)
    x <- lapply(x, function(x) as.data.frame(as.list(x)))

  return(x)
}

postprocessModels <- function(models, library) {
  if ("delay" %in% names(models)) {
    if ("transitCompartments" %in% models$absorption)
      models[models$absorption == "transitCompartments", ]$delay <- NA
  }

  return(models)
}
