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
#' @param algorithm (optional) Algorithm to use. Options are "ACO", "tournament", "exhaustive_search", "decision_tree". If not provided, it's automatically selected: "ACO" if searching for IIV, and "decision_tree" otherwise.
#' @param settings (optional) A list of settings for the algorithms.
#'   - N: Number of models to test per iteration (ants/wins for ACO/Tournament). Default depends on the search space size.
#'   - findIIV: Boolean, whether to explore interindividual variability. Default: FALSE.
#'   - findError: Boolean, whether to explore the residual error model. Default: FALSE.
#'   - seed: A random seed for reproducibility. Default: 123456.
#'   - previous_runs: An optional data frame of previous model runs to learn from.
#'   - output: Optional file path to save results incrementally (with the csv extension). If not provided, the results are not saved incrementally.
#'   - convergence_threshold (ACO): The proportion of ants that must select the same model to trigger convergence (e.g., 0.75). Default: 0.75.
#'   - initial_iiv_forced_iter (ACO/Tournament): The number of initial iterations where IIV is forced on all parameters to encourage exploration (e.g., 3). Default: 3.
#'   - rho (ACO/Tournament): Evaporation (learning) rate. Default depends on the search space size.
#'   - obsIDToUse: observation ID to use (required if multiple observation IDs).
#'   - plot: plot the progress (TRUE/FALSE). Default: TRUE.
#'   - save_mode: Character string indicating how model runs should be saved during training.
#'     Must be one of:
#'     \itemize{
#'       \item \code{"none"}: Do not save any model runs.
#'       \item \code{"best"} (default): Save only the best-performing model based.
#'       \item \code{"all"}: Save all model runs.
#'     }
#'     Models are saved in the autoBuild subdirectory of the result folder of the project.
#' @param metric_func (optional) A function that returns the performance metric for the last run model (based on lixoftConnectors). By default, BICc if no IIV, BICc + count(RSE>50%) * number_of_observations, if IIV.
#' @param initial_func (optional) A function to return initial estimates. The pkpopini function from Rsmlx by default.
#'
#' @return A data frame containing all models tested and their performance metrics.
#' @examples
#' \dontrun{
#' # Example 1: Default search with oral clearance models and additional filters
#' findModel(
#'   project = file.path(getDemoPath(), "1.creating_and_using_models",
#'                       "1.1.libraries_of_models", "theophylline_project.mlxtran"),
#'   library = "pk",
#'   requiredFilters = list(administration = "oral", parametrization = "clearance"),
#'   filters = list(delay = c("noDelay", "lagTime"), elimination = "linear"),
#'   settings = list(iiv = FALSE, error = TRUE, save_mode = "best")
#' )
#'
#' # Example 2: Oral bolus rate-parametrized models on an IV/Oral project
#' findModel(
#'   project = file.path(getDemoPath(), "6.PK_models", "6.2.multiple_routes",
#'                       "ivOral1_project.mlxtran"),
#'   library = "pk",
#'   requiredFilters = list(administration = "oralBolus", parametrization = "rate"),
#'   filters = list(delay = c("noDelay", "lagTime"), elimination = "linear"),
#'   settings = list(iiv = FALSE, error = FALSE)
#' )
#'
#' # Example 3: Bolus route with IIV and residual error estimation
#' findModel(
#'   project = file.path(getDemoPath(), "6.PK_models", "6.1.single_route",
#'                       "bolusMixed_project.mlxtran"),
#'   library = "pk",
#'   requiredFilters = list(administration = "bolus", parametrization = "rate"),
#'   settings = list(iiv = TRUE, error = TRUE, output = "list_of_runs.csv")
#' )
#'
#' # Example 4: Exhaustive search over oral clearance models
#' findModel(
#'   project = file.path(getDemoPath(), "1.creating_and_using_models",
#'                       "1.1.libraries_of_models", "theophylline_project.mlxtran"),
#'   library = "pk",
#'   requiredFilters = list(administration = "oral", parametrization = "clearance"),
#'   algorithm = "exhaustive_search",
#'   settings = list(output = "list_of_runs.csv")
#' )
#'
#' # Example 5: Decision tree search strategy
#' findModel(
#'   project = file.path(getDemoPath(), "1.creating_and_using_models",
#'                       "1.1.libraries_of_models", "theophylline_project.mlxtran"),
#'   library = "pk",
#'   requiredFilters = list(administration = "oral", parametrization = "clearance"),
#'   algorithm = "decision_tree"
#' )
#' }
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
findModel <- function(project, library, requiredFilters, filters = NULL,
                      algorithm = NULL, settings = NULL, metric_func = NULL,
                      initial_func = NULL) {

  op = options()
  op$lixoft_notificationOptions$warnings = 1
  op$lixoft_notificationOptions$info = 1
  op$cli.progress_show_after = 0
  options(op)

  names(settings)[names(settings) == "findIIV"] <- "iiv"
  names(settings)[names(settings) == "findError"] <- "error"

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
  if (requiredFilters$administration %in% c("oral", "oralBolus"))
    keys <- c("distribution", "absorption_delay", "distribution", "elimination")
  else
    keys <- c("delay", "distribution", "elimination")

  default_settings <- list(
    N = ifelse(is.null(settings$iiv) || !settings$iiv, 10, 20),
    iiv = FALSE,
    error = FALSE,
    seed = 123456,
    previous_runs = NULL,
    output = NULL,
    convergence_threshold = 0.75,
    initial_iiv_forced_iter = 3,
    rho = ifelse(is.null(settings$iiv) || !settings$iiv, 0.4, 0.6),
    save_mode = "best",
    keys = keys,
    plot = TRUE
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

  # --- Do not plot diagrams if DiagrammeR not installed ---
  if (settings$plot && algorithm == "decision_tree") {
    settings$plot <- isTRUE(requireNamespace("DiagrammeR", quietly = TRUE))
    if (!settings$plot)
      warning("Package 'DiagrammeR' is not installed, so the diagram will not be plotted.", call. = FALSE)
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

  # Load and clean project
  lixoftConnectors::loadProject(project)

  validate_obsidtouse(settings)

  init_var <- lixoftConnectors::getIndividualParameterModel()$variability
  init_var <- lapply(init_var, function(x) {
    x[] <- TRUE
    return(x)
  })

  init_cor <- lixoftConnectors::getIndividualParameterModel()$correlationBlocks
  init_cor <- lapply(init_cor, function(x) list())

  init_cov <- lixoftConnectors::getIndividualParameterModel()$covariateModel
  init_cov <- lapply(init_cov, function(x) {
    x[] <- FALSE
    return(x)
  })

  lixoftConnectors::setIndividualParameterModel(list(
    variability = init_var,
    correlationBlocks = init_cor,
    covariateModel = init_cov
  ))

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

