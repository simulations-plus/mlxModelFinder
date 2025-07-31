#' Run Model Selection
#'
#' @param project (required) Path to the Monolix project with the loaded data set (SAEM settings, seed will be taken from the project).
#' @param library (required) The model library to use. Options: "pk", "tmdd", "pkpd".
#' @param filters (optional) A list defining the structural model components to explore. For the "pk" library, the options are:
#' \itemize{
#'   \item \strong{delay}: "noDelay", "lagTime", "transitCompartments"
#'   \item \strong{absorption}: "zeroOrder", "firstOrder", "sigmoid"
#'   \item \strong{distribution}: "1compartment", "2compartments", "3compartments"
#'   \item \strong{elimination}: "linear", "MichaelisMenten", "combined"
#'   \item \strong{bioavailability} (if administration is "oralBolus"): "true", "false"
#' }
#' If a filter is not specified, all its possible options will be included in the search.
#' @param requiredFilters (required) A list of required structural model components. For PK models, options are: administration ("bolus", "infusion", "oral", "oralBolus") and parametrization ("rate", "clearance").
#' @param algorithm (optional) Algorithm to use. Options are "ACO", "tournament", "exhaustive_search", "decision_tree". If not provided, it's automatically selected: "ACO" for search spaces larger than 15,000 models, and "decision_tree" otherwise.
#' @param settings (optional) A list of settings for the algorithms.
#'   - N: Number of models to test per iteration (ants/wins for ACO/Tournament).
#'   - iiv: Boolean, whether to explore interindividual variability.
#'   - error: Boolean, whether to explore the residual error model.
#'   - seed: A random seed for reproducibility.
#'   - previous_runs: An optional data frame of previous model runs to learn from.
#'   - output: Optional file path to save results incrementally (with the csv extension).
#'   - convergence_threshold (ACO): The proportion of ants that must select the same model to trigger convergence (e.g., 0.75).
#'   - initial_iiv_forced_iter (ACO/Tournament): The number of initial iterations where IIV is forced on all parameters to encourage exploration (e.g., 3).
#'   - rho (ACO/Tournament): Evaporation (learning) rate.
#'   - save_mode: Character string indicating how model runs should be saved during training.
#'     Must be one of:
#'     \itemize{
#'       \item \code{"none"} (default): Do not save any model runs.
#'       \item \code{"best"}: Save only the best-performing model based.
#'       \item \code{"all"}: Save all model runs.
#'     }
#' @param metric_func (optional) A function that returns the performance metric for the last run model (based on lixoftConnectors). By default, BICc if no IIV, BICc + count(RSE>50%) * number_of_observations, if IIV.
#' @param initial_func (optional) A function to return initial estimates. The pkpopini function from Rsmlx by default.
#'
#' @return A data frame containing all models tested and their performance metrics.
#'
#' @importFrom stats aggregate lm optim setNames
#' @importFrom utils read.csv write.csv
#' @importFrom ggplot2 ggplot aes geom_line geom_point facet_wrap
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 labs theme
#' @importFrom cli cli_alert_success cli_progress_bar
#' @importFrom cli cli_progress_done cli_progress_update
#' @importFrom gridExtra grid.arrange
#'
#' @export
runStructuralModelBuilding <- function(project, library, requiredFilters,
                                       filters = NULL, algorithm = NULL,
                                       settings = NULL, metric_func = NULL,
                                       initial_func = NULL) {

  op = options()
  op$lixoft_notificationOptions$warnings = 1
  op$lixoft_notificationOptions$info = 1
  op$cli.progress_show_after = 0
  options(op)

  # --- Argument Checks ---
  validate_project(project)
  validate_library(library)
  validate_required_filters(requiredFilters)
  if (!is.null(algorithm)) {
    validate_algorithm(algorithm)
  }
  if (!is.null(filters))
    validate_filters(filters, requiredFilters)
  if (!is.null(settings))
    validate_settings(settings)
  if (!is.null(metric_func))
    validate_no_argument_function(metric_func, "metric_func")
  if (!is.null(initial_func))
    validate_no_argument_function(initial_func, "initial_func")

  # --- Default Settings ---
  default_settings <- list(
    N = ifelse(is.null(settings$iiv) || !settings$iiv, 10, 20),
    iiv = FALSE,
    error = FALSE,
    seed = 123456,
    previous_runs = NULL,
    output = NULL,
    convergence_threshold = 0.75,
    initial_iiv_forced_iter = 5,
    rho = ifelse(is.null(settings$iiv) || !settings$iiv, 0.4, 0.6),
    save_mode = "best",
    keys = c("distribution", "absorption_delay", "distribution", "elimination")
  )

  # Overwrite defaults with user-provided settings
  if (!is.null(settings)) {
    default_settings[names(settings)] <- settings
  }
  settings <- default_settings

  # --- Setup ---
  set.seed(settings$seed)
  filters <- updateFilters(library, filters, requiredFilters)
  space <- calculateSearchSpace(library, filters, settings$error, settings$iiv)
  cli::cli_alert_info("Search space size: {space} models")

  # --- Algorithm Selection ---
  if (is.null(algorithm)) {
    if (space > 15000) {
      algorithm <- "ACO"
    } else {
      algorithm <- "decision_tree"
    }
    cli::cli_alert_info("Algorithm automatically selected: {algorithm}")
  }

  # --- Selection of Metric Function ---
  if (is.null(metric_func)) {
    if (settings$iiv) {
      metric_func <- calculateMetricRSEs
    } else {
      metric_func <- calculateMetricBICc
    }
  }

  lixoftConnectors::initializeLixoftConnectors(force = TRUE)

  lixoftConnectors::loadProject(project)

  settings$result_folder <- lixoftConnectors::getProjectSettings()$directory

  if (settings$save_mode != "none") {
    dir.create(file.path(settings$result_folder, "autoBuild"), showWarnings = FALSE)
    lixoftConnectors::setProjectSettings(dataandmodelnexttoproject = TRUE)
  }

  if (!is.null(settings$previous_runs))
    settings$previous_runs <- read.csv(settings$previous_runs)

  switch(algorithm,
         ACO = {
           return(runACO(library, filters, requiredFilters, settings, metric_func))
         },
         tournament = {
           return(runTournament(library, filters, requiredFilters, settings, metric_func))
         },
         # random_forest = {
         #   return(runRandomForest(library, filters, requiredFilters, settings, metric_func))
         # },
         exhaustive_search = {
           return(runExhaustiveSearch(library, filters, requiredFilters, settings, metric_func))
         },
         decision_tree = {
           return(runDecisionTree(library, filters, requiredFilters, settings, metric_func))
         })
}

