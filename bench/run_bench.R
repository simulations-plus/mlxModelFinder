devtools::load_all()

projects <- list.files("bench/projects/", pattern="\\.mlxtran$", full.names = TRUE)

projects <- projects[2:length(projects)]

seeds <- 1:5
subdirectory <- "elitism"
settings <- list(
  findIIV = TRUE,
  findError = TRUE,
  max_models = 100,
  elitism = TRUE,
  plot = FALSE
)


results <- NULL

for (project in projects) {
  for (seed in seeds) {
    results_file_name <- paste0("bench/run_lists/", tools::file_path_sans_ext(basename(project)), ".csv")

    message(tools::file_path_sans_ext(basename(project)))

    if (file.exists(results_file_name))
      previous_runs <- results_file_name
    else
      previous_runs <- NULL

    default_settings <- list(
      previous_runs = previous_runs,
      output = results_file_name,
      seed = seed
    )

    res <- findModel(
      project = project,
      library = "pk",
      requiredFilters = list(administration = "oral", parametrization = "clearance"),
      filters = list(absorption = c("firstOrder", "zeroOrder", "sigmoid"), elimination = "linear"),
      settings = append(settings, default_settings)
    )

    path <- paste0("bench/results/", subdirectory)

    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
    }

    write.csv(res$iteration_progress, paste0(path, "/", tools::file_path_sans_ext(basename(project)), "_", seed, ".csv"), row.names = FALSE)
  }
}
