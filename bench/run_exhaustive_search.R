library(lixoftConnectors)

initializeLixoftConnectors()

devtools::load_all()

projects <- list.files("bench/projects/", pattern="\\.mlxtran$", full.names = TRUE)

settings <- list(
    findIIV = TRUE,
    findError = TRUE
  )

results <- NULL

for (project in projects) {
  results_file_name <- paste0("bench/exhaustive_search/", tools::file_path_sans_ext(basename(project)), ".csv")

  message(tools::file_path_sans_ext(basename(project)))

  new_path <- paste0("bench/exhaustive_projects/", tools::file_path_sans_ext(basename(project)), ".mlxtran")

  loadProject(project)
  saveProject(new_path)

  default_settings <- list(
    output = results_file_name
  )

  res <- findModel(
    project = new_path,
    library = "pk",
    requiredFilters = list(administration = "oral", parametrization = "clearance"),
    filters = list(absorption = c("firstOrder", "zeroOrder", "sigmoid"), elimination = "linear"),
    algorithm = "exhaustive_search",
    settings = append(settings, default_settings)
  )
}
