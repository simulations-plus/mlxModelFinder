#' Run Model Selection (PK, Monolix)
#'
#' Automates structural model selection for PK models in Monolix via
#' `lixoftConnectors`. The function constructs a search space from
#' `requiredFilters` (route/parametrization) and optional `filters`
#' (absorption, delay, distribution, elimination, bioavailability), then
#' explores that space using a chosen algorithm and returns a tidy data frame
#' of all tested models with their performance metrics and run metadata.
#'
#' If `algorithm` is not provided, it is selected automatically:
#' **"ACO"** when searching for IIV (`settings$iiv = TRUE`), otherwise
#' **"decision_tree"**.
#'
#' @param project (required) Path to the **Monolix** project (`.mlxtran`) that already
#'   contains the dataset and SAEM settings. Estimation settings/seed are taken
#'   from the project.
#' @param library (required) Model library to use. Options: `"pk"` (currently supported).
#' @param filters (optional) Named list defining **optional** structural components
#'   to explore. For the `"pk"` library:
#'   \itemize{
#'     \item \strong{delay}: `"noDelay"`, `"lagTime"`
#'     \item \strong{absorption}: `"zeroOrder"`, `"firstOrder"`, `"sigmoid"`, `"transitCompartments"`
#'     \item \strong{distribution}: `"1compartment"`, `"2compartments"`, `"3compartments"`
#'     \item \strong{elimination}: `"linear"`, `"MichaelisMenten"`, `"combined"`
#'     \item \strong{bioavailability} (only if \code{requiredFilters$administration == "oralBolus"}):
#'           `"true"`, `"false"`
#'   }
#'   If a component is omitted, **all** of its options are included in the search.
#' @param requiredFilters (required) Named list of **required** structural components:
#'   \itemize{
#'     \item \strong{administration}: `"bolus"`, `"infusion"`, `"oral"`, `"oralBolus"`
#'     \item \strong{parametrization}: `"rate"`, `"clearance"`
#'   }
#' @param algorithm (optional) Search algorithm. One of
#'   `"ACO"`, `"tournament"`, `"exhaustive_search"`, `"decision_tree"`.
#'   If `NULL`, defaults to `"ACO"` when `settings$iiv = TRUE`, otherwise `"decision_tree"`.
#' @param settings (optional) Named list of search/algorithm settings. Unless stated,
#'   defaults are chosen adaptively (see below). Supported entries:
#'   \itemize{
#'     \item \strong{N}: Integer; models per iteration (ants for ACO; wins/cohort size for Tournament).
#'           Default: `10` (no IIV) or `20` (with IIV).
#'     \item \strong{findIIV}: Logical; explore interindividual variability (IIV). Default: `FALSE`.
#'     \item \strong{findError}: Logical; explore residual error model. Default: `FALSE`.
#'     \item \strong{seed}: Integer; RNG seed for reproducibility. Default: `123456`.
#'     \item \strong{previous_runs}: Optional path to a CSV of past runs to warm-start/learn from. Default: `NULL`.
#'     \item \strong{output}: Optional CSV path; if provided, results are appended incrementally. Default: `NULL`.
#'     \item \strong{initial_iiv_forced_iter} (ACO/Tournament): Integer; number of initial iterations
#'           where IIV is forced on all parameters to encourage exploration. Default: `3`.
#'     \item \strong{rho} (ACO/Tournament): Learning/evaporation rate. Default: `0.4` (no IIV) or `0.6` (with IIV).
#'     \item \strong{obsIDToUse}: Character/integer; which observation ID to use (required if the project has multiple).
#'     \item \strong{plot}: Logical; if `TRUE` and optional deps are available, render progress
#'           (DiagrammeR/data.tree for Decision Tree; plotly/htmltools for ACO/Tournament). Default: `TRUE`.
#'     \item \strong{linearization}: Logical; if `TRUE`, use linearization for log-likelihood in the default metric.
#'           Default: `FALSE`.
#'     \item \strong{save_mode}: How to persist Monolix runs during training. One of:
#'       \itemize{
#'         \item \code{"none"}: do not save runs
#'         \item \code{"best"} (default): save only the current best-performing run
#'         \item \code{"all"}: save every run
#'       }
#'       Models (if saved) are written under \code{autoBuild/} in the project's result directory.
#'     \item \strong{clip}: Logical; clip sampling probabilities during stochastic search. Default: `FALSE`.
#'     \item \strong{keys}: Character vector; internal ordering of structural “keys” to update
#'           during decision tree (auto-derived from administration).
#'     \item \strong{elitism}: Logical; keep the current best candidate across iterations (ACO).
#'           Default: `TRUE`.
#'     \item \strong{max_models}: Integer cap on the total number of models to evaluate. Default: `1e6`.
#'     \item \strong{iiv_favor_true}: Logical; bias exploration toward keeping IIV Default: `FALSE`.
#'     \item \strong{xgboost}: Logical; enable XGBoost-based guidance for ACO. Default: `FALSE`. Limited testing performed.
#'     \item \strong{xgboost_model_nb}: Integer; number of boosting models to maintain if \code{xgboost = TRUE}. Default: `5`.
#'     \item \strong{local_search}: Logical; enable optional local IIV refinement around good candidates (ACO). Default: `FALSE`.
#'   }
#' @param metric_func (optional) **Zero-argument** function that returns a single numeric
#'   score for the **last run model** (lower is better). Use `lixoftConnectors` inside
#'   the function to read results (e.g., information criteria). **Default**:
#'   \emph{BICc} when `settings$iiv = FALSE`; \emph{BICc + penalty} when `settings$iiv = TRUE`,
#'   where the penalty is \code{count(RSE > 50\%) * log(n_observations)}. This default was chosen
#'   after thorough testing for robustness and precision.
#' @param initial_func (optional) Function that **takes one argument** (the \code{project}
#'   path) and **returns a project path** (same or new) after setting initial values via
#'   `lixoftConnectors`. Use this to inject custom population initials before each run.
#'   If `NULL`, the package uses its internal/default initializer.
#'
#' @return A data frame where each row corresponds to a tested model, including:
#'   structural choices, IIV/error settings, convergence status, metric value, timing,
#'   and (if enabled) paths to saved Monolix runs.
#'
#' @details
#' The function initializes `lixoftConnectors`, loads the provided project, validates inputs,
#' constructs the search space from `requiredFilters` + `filters`, and then runs the selected
#' algorithm. Progress and diagrams are only plotted if `settings$plot = TRUE` and the relevant
#' suggested packages are installed.
#'
#' @examples
#' \dontrun{
#' # Example 1: Default search with oral clearance models and additional filters
#' findModel(
#'   project = file.path(getDemoPath(), "1.creating_and_using_models",
#'                       "1.1.libraries_of_models", "theophylline_project.mlxtran"),
#'   library = "pk",
#'   requiredFilters = list(administration = "oral", parametrization = "clearance"),
#'   filters = list(absorption = c("zeroOrder", "firstOrder"), elimination = "linear"),
#'   settings = list(findIIV = FALSE, findError = TRUE)
#' )
#'
#' # Example 2: Oral bolus rate-parametrized models on an IV/Oral project w/o plotting progress
#' findModel(
#'   project = file.path(getDemoPath(), "6.PK_models", "6.2.multiple_routes",
#'                       "ivOral1_project.mlxtran"),
#'   library = "pk",
#'   requiredFilters = list(administration = "oralBolus", parametrization = "rate"),
#'   filters = list(absorption = c("zeroOrder", "firstOrder"),
#'                  elimination = "linear", bioavailability = "true"),
#'   settings = list(findIIV = FALSE, findError = FALSE, plot = FALSE)
#' )
#'
#' # Example 3: Bolus route with IIV and residual error estimation, saving progress
#' findModel(
#'   project = file.path(getDemoPath(), "6.PK_models", "6.1.single_route",
#'                       "bolusMixed_project.mlxtran"),
#'   library = "pk",
#'   requiredFilters = list(administration = "bolus", parametrization = "rate"),
#'   settings = list(findIIV = TRUE, findError = TRUE, output = "list_of_runs.csv")
#' )
#'
#' # Example 4: Exhaustive search over oral clearance models
#' findModel(
#'   project = file.path(getDemoPath(), "1.creating_and_using_models",
#'                       "1.1.libraries_of_models", "theophylline_project.mlxtran"),
#'   library = "pk",
#'   requiredFilters = list(administration = "oral", parametrization = "clearance"),
#'   algorithm = "exhaustive_search"
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
#' @importFrom stats aggregate lm optim setNames as.formula model.matrix predict
#' @importFrom utils read.csv write.csv
#' @importFrom cli cli_alert_success cli_progress_bar
#' @importFrom cli cli_progress_done cli_progress_update
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

  if (!is.null(settings)) {
    names(settings)[names(settings) == "findIIV"] <- "iiv"
    names(settings)[names(settings) == "findError"] <- "error"
  }

  # --- Argument Checks ---
  validate_project(project)
  validate_library(library)
  validate_required_filters(requiredFilters)
  if (!is.null(algorithm)) validate_algorithm(algorithm)
  if (!is.null(filters))   validate_filters(filters, requiredFilters)

  if (!is.null(metric_func))  validate_metric_func_no_args(metric_func)
  if (!is.null(initial_func)) validate_initial_func_one_arg(initial_func)

  # --- Default Settings ---
  keys <- switch(requiredFilters$administration,
     "oral" = c("distribution", "absorption_delay", "distribution", "elimination"),
     "oralBolus" = c("bioavailability", "distribution", "absorption_delay", "distribution", "elimination"),
     c("distribution", "elimination")
  )

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
    plot = TRUE,
    linearization = FALSE,
    elitism = TRUE,
    max_models = 1000000,
    iiv_favor_true = FALSE,
    clip = FALSE,
    xgboost = FALSE,
    xgboost_model_nb = 5,
    local_search = FALSE
  )

  # --- Setup ---
  set.seed(settings$seed)
  filters <- updateFilters(library, filters, requiredFilters)
  space <- calculateSearchSpace(library, filters, settings$error, settings$iiv)
  cli::cli_alert_info("Search space size: {space} models")

  # --- Algorithm Selection ---
  if (is.null(algorithm)) {
    algorithm <- if (isTRUE(settings$iiv)) "ACO" else "decision_tree"
    cli::cli_alert_info("Algorithm automatically selected: {algorithm}")
  }

  validate_settings(settings, algorithm)

  # --- Override default settings ---
  if (!is.null(settings)) {
    default_settings[names(settings)] <- settings
  }
  settings <- default_settings

  # --- Do not plot diagrams / progress if packages not installed ---
  if (isTRUE(settings$plot) && identical(algorithm, "decision_tree")) {
    settings$plot <- .check_pkgs(c("DiagrammeR", "data.tree"), what = "diagram")
  }

  if (isTRUE(settings$plot) && algorithm %in% c("ACO", "tournament")) {
    settings$plot <- .check_pkgs(c("plotly", "htmltools"), what = "progress plots")
  }

  if (isTRUE(settings$xgboost) && algorithm %in% c("ACO", "tournament")) {
    settings$xgboost <- .check_pkgs(c("xgboost"), what = "XGBoost")
  }

  # --- Selection of Metric Function ---
  if (is.null(metric_func))
    metric_func <- function() calculateMetric(RSEs = settings$iiv,
                                              linearization = settings$linearization)

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
           return(runACO(library, filters, requiredFilters, settings, metric_func, initial_func))
         },
         tournament = {
           return(runTournament(library, filters, requiredFilters, settings, metric_func, initial_func))
         },
         exhaustive_search = {
           return(runExhaustiveSearch(library, filters, requiredFilters, settings, metric_func, initial_func))
         },
         decision_tree = {
           return(runDecisionTree(library, filters, requiredFilters, settings, metric_func, initial_func))
         })
}

