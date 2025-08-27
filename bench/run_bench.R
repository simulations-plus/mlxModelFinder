devtools::load_all()

projects <- list.files("bench/projects/", pattern="\\.mlxtran$", full.names = TRUE)

seeds <- 1:5
subdirectories <- c("classic", "clip", "elitism", "elitism_clip", "favor_true", "favor_true_1", "favor_true_0", "favor_true_0_elitism", "favor_true_1_elitism")
all_settings <- list(
  list(
    findIIV = TRUE,
    findError = TRUE,
    elitism = FALSE,
    plot = FALSE,
    clip = FALSE,
    iiv_favor_true = FALSE
  ),
  list(
    findIIV = TRUE,
    findError = TRUE,
    elitism = FALSE,
    plot = FALSE,
    clip = TRUE,
    iiv_favor_true = FALSE
  ),
  list(
    findIIV = TRUE,
    findError = TRUE,
    elitism = TRUE,
    plot = FALSE,
    clip = FALSE,
    iiv_favor_true = FALSE
  ),
  list(
    findIIV = TRUE,
    findError = TRUE,
    elitism = TRUE,
    plot = FALSE,
    clip = TRUE,
    iiv_favor_true = FALSE
  ),
  list(
    findIIV = TRUE,
    findError = TRUE,
    elitism = FALSE,
    plot = FALSE,
    clip = FALSE,
    iiv_favor_true = TRUE
  ),
  list(
    findIIV = TRUE,
    findError = TRUE,
    elitism = FALSE,
    plot = FALSE,
    clip = FALSE,
    iiv_favor_true = TRUE,
    initial_iiv_forced_iter = 1
  ),
  list(
    findIIV = TRUE,
    findError = TRUE,
    elitism = FALSE,
    plot = FALSE,
    clip = FALSE,
    iiv_favor_true = TRUE,
    initial_iiv_forced_iter = 0
  ),
  list(
    findIIV = TRUE,
    findError = TRUE,
    elitism = TRUE,
    plot = FALSE,
    clip = FALSE,
    iiv_favor_true = TRUE,
    initial_iiv_forced_iter = 0
  ),
  list(
    findIIV = TRUE,
    findError = TRUE,
    elitism = TRUE,
    plot = FALSE,
    clip = FALSE,
    iiv_favor_true = TRUE,
    initial_iiv_forced_iter = 1
  )
)


results <- NULL

for (i in c(2, 4, 9)) {
  subdirectory <- subdirectories[i]
  settings <- all_settings[[i]]

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
}
