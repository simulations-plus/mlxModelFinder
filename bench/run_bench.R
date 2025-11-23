devtools::load_all()

projects <- list.files("bench/projects/", pattern="\\.mlxtran$", full.names = TRUE)
projects <- projects[c(3)]

seeds <- 1:10
subdirectories <- c("rho_06", "clip", "elitism", "elitism_clip", "favor_true",
                    "favor_true_1", "favor_true_0", "favor_true_0_elitism",
                    "favor_true_1_elitism", "favor_true_1_nostrong", "rho_02",
                    "rho_04", "rho_08", "local_search", "xgboost", "elitism_xgboost",
                    "elitism_xgboost_15", "elitism_xgboost_iiv1")
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
    iiv_favor_true = FALSE,
    rho = 0.4
  ),
  list(
    findIIV = TRUE,
    findError = TRUE,
    elitism = TRUE,
    plot = FALSE,
    clip = FALSE,
    iiv_favor_true = FALSE,
    rho = 0.4
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
    iiv_favor_true = FALSE,
    rho = 0.2
  ),
  list(
    findIIV = TRUE,
    findError = TRUE,
    elitism = FALSE,
    plot = FALSE,
    clip = FALSE,
    iiv_favor_true = FALSE,
    rho = 0.4
  ),
  list(
    findIIV = TRUE,
    findError = TRUE,
    elitism = FALSE,
    plot = FALSE,
    clip = FALSE,
    iiv_favor_true = FALSE,
    rho = 0.8
  ),
  list(
    findIIV = TRUE,
    findError = TRUE,
    elitism = FALSE,
    plot = FALSE,
    clip = FALSE,
    iiv_favor_true = FALSE,
    rho = 0.4,
    local_search = TRUE
  ),
  list(
    findIIV = TRUE,
    findError = TRUE,
    elitism = FALSE,
    plot = FALSE,
    clip = FALSE,
    iiv_favor_true = FALSE,
    rho = 0.4,
    xgboost = TRUE
  ),
  list(
    findIIV = TRUE,
    findError = TRUE,
    elitism = TRUE,
    plot = FALSE,
    clip = FALSE,
    iiv_favor_true = FALSE,
    rho = 0.4,
    xgboost = TRUE
  ),
  list(
    findIIV = TRUE,
    findError = TRUE,
    elitism = TRUE,
    plot = FALSE,
    clip = FALSE,
    iiv_favor_true = FALSE,
    rho = 0.4,
    xgboost = TRUE,
    xgboost_model_nb = 15
  ),
  list(
    findIIV = TRUE,
    findError = TRUE,
    elitism = TRUE,
    plot = FALSE,
    clip = FALSE,
    iiv_favor_true = FALSE,
    rho = 0.4,
    xgboost = TRUE,
    initial_iiv_forced_iter = 1
  )
)


results <- NULL

for (i in c(1:3, 11:16)) {
# for (i in 18) {
  subdirectory <- subdirectories[i]
  settings <- all_settings[[i]]

  for (project in projects) {
    for (seed in seeds) {
      results_file_name <- paste0("bench/exhaustive_search/", tools::file_path_sans_ext(basename(project)), ".csv")

      message(tools::file_path_sans_ext(basename(project)))

      if (file.exists(results_file_name))
        previous_runs <- results_file_name
      else
        previous_runs <- NULL

      default_settings <- list(
        previous_runs = previous_runs,
        # output = results_file_name,
        seed = seed
      )

      res <- findModel(
        project = project,
        library = "pk",
        requiredFilters = list(administration = "oral", parametrization = "clearance"),
        filters = list(absorption = c("firstOrder", "zeroOrder", "sigmoid"), elimination = "linear"),
        settings = append(settings, default_settings)
      )

      path <- paste0("bench/results/exhaustive_comp/", subdirectory)

      if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
      }

      write.csv(res$iteration_progress, paste0(path, "/", tools::file_path_sans_ext(basename(project)), "_", seed, ".csv"), row.names = FALSE)
    }
  }
}
