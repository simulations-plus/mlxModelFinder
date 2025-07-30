probabilities <- function(x) {
  filters <- sapply(x, '[', seq(max(sapply(x, length))))
  weights <- apply(filters, 1:2, function(x) as.numeric(!is.na(x)))
  probas <- apply(weights, 2, function(x) x/sum(x, na.rm = TRUE))

  p <- list(probas = probas, filters = filters, weights = weights)
  return(p)
}

iiv_probabilities <- function(library) {
  if (library == "pk") {
    param_list <- c("ka", "Tk0", "Tlag", "Ktr", "Mtt", "V1", "Cl", "Q2", "V2", "Q3", "V3", "k", "k12", "k21", "k13", "k31", "Vm", "Km", "F")
    param_map <- lapply(param_list, function(x) c(TRUE, FALSE))
    names(param_map) <- param_list
  }
  return(probabilities(param_map))
}

updateProbabilities <- function(obj, runs, rho = 0.4) {
  obj$weights <- obj$weights * (1 - rho)
  keys <- colnames(obj$weights)
  for (key in keys)
    for (i in 1:nrow(obj$weights))
      obj$weights[i, key] <- obj$weights[i, key] + sum(1 / runs[runs[, key] == obj$filters[i, key], "rank"], na.rm = TRUE)
  obj$probas <- apply(obj$weights, 2, function(x) x/sum(x, na.rm = TRUE))
  return(obj)
}
