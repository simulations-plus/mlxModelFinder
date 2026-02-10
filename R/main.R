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
#' @param library (required) Model library to use. Options: `"pk"` (built-in PK library),
#'   `"custom"` (user-defined models via `model_creation_func`).
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
#' @param requiredFilters (optional for `library = "custom"`, required for `library = "pk"`)
#'   Named list of **required** structural components for the `"pk"` library:
#'   \itemize{
#'     \item \strong{administration}: `"bolus"`, `"infusion"`, `"oral"`, `"oralBolus"`
#'     \item \strong{parametrization}: `"rate"`, `"clearance"`
#'   }
#'   Not used when `library = "custom"`.
#' @param algorithm (optional) Search algorithm. One of
#'   `"ACO"`, `"exhaustive_search"`, `"decision_tree"`.
#'   If `NULL`, defaults to `"ACO"` when `settings$iiv = TRUE`, otherwise `"decision_tree"`.
#' @param settings (optional) Named list of search/algorithm settings. Unless stated,
#'   defaults are chosen adaptively (see below). Supported entries:
#'   \itemize{
#'     \item \strong{N}: Integer; models per iteration (ants for ACO).
#'           Default: `10` (no IIV) or `20` (with IIV).
#'     \item \strong{findIIV}: Logical; explore interindividual variability (IIV). Default: `FALSE`.
#'     \item \strong{findError}: Logical; explore residual error model. Default: `FALSE`.
#'     \item \strong{seed}: Integer; RNG seed for reproducibility. Default: `123456`.
#'     \item \strong{previous_runs}: Optional path to a CSV of past runs to warm-start/learn from. Default: `NULL`.
#'     \item \strong{output}: Optional CSV path; if provided, results are appended incrementally. Default: `NULL`.
#'     \item \strong{initial_iiv_forced_iter} (ACO): Integer; number of initial iterations
#'           where IIV is forced on all parameters to encourage exploration. Default: `3`.
#'     \item \strong{rho} (ACO): Learning/evaporation rate. Default: `0.4`.
#'     \item \strong{obsIDToUse}: Character/integer; which observation ID to use (required if the project has multiple).
#'     \item \strong{plot}: Logical; if `TRUE` and optional deps are available, render progress
#'           (DiagrammeR/data.tree for Decision Tree; plotly/htmltools for ACO). Default: `TRUE`.
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
#'     \item \strong{xgboost}: Logical; enable XGBoost-based guidance for ACO. Default: `TRUE`.
#'     \item \strong{xgboost_model_nb}: Integer; number of boosting models to maintain if \code{xgboost = TRUE}. Default: `5`.
#'     \item \strong{table_func}: function that allows you to output additional information about each run. See online doc for more information.
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
#' @param model_creation_func (optional) **Custom model creation function** for use with
#'   `library = "custom"`. Must be a function that accepts exactly two arguments:
#'   `library` (character) and `filters` (named list), and returns a character string
#'   containing the **absolute path** to a valid Monolix model file. This allows you to
#'   use your own model library instead of the built-in PK library. When `library = "custom"`,
#'   this parameter is **required**. Only compatible with `algorithm = "exhaustive_search"`
#'   or `algorithm = "ACO"` (not compatible with `"decision_tree"`).
#' @param param_mapping (optional) **Parameter mapping structure** for custom models with IIV.
#'   Required when `library = "custom"` and `settings$iiv = TRUE`. Must be a named list with:
#'   \itemize{
#'     \item \strong{common_params}: Character vector of parameter names present in all models
#'     \item \strong{filter_params}: Named list where each name corresponds to a filter category
#'           (e.g., "absorption"), and each value is a named list mapping filter values
#'           (e.g., "firstOrder") to character vectors of parameter names specific to that filter value
#'   }
#'   This structure tells the algorithm which parameters are available for IIV exploration
#'   in each model variant. See vignette for detailed examples.
#'
#' @return A list that contains:
#'   \itemize{
#'     \item **results**: a data frame where each row corresponds to a tested model, including:
#'   structural choices, IIV/error settings, metric value,
#'     \item **best_model**: a data frame with one row corresponding to the best model,
#'     \item **best_path**: path to the saved Monolix runs corresponding to the best model (if save_mode is `all` or `best`),
#'     \item **iteration_progress**: progress per iteration for ACO,
#'     \item **prob_history**: list of sampling probabilities at the end of each iteration for ACO.
#'   }
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
#' res <- findModel(
#'   project = file.path(getDemoPath(), "1.creating_and_using_models",
#'                       "1.1.libraries_of_models", "theophylline_project.mlxtran"),
#'   library = "pk",
#'   requiredFilters = list(administration = "oral", parametrization = "clearance"),
#'   filters = list(absorption = c("zeroOrder", "firstOrder"), elimination = "linear"),
#'   settings = list(findIIV = FALSE, findError = TRUE)
#' )
#'
#' # Example 2: Oral bolus rate-parametrized models on an IV/Oral project w/o plotting progress
#' res <- findModel(
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
#' res <- findModel(
#'   project = file.path(getDemoPath(), "6.PK_models", "6.1.single_route",
#'                       "bolusMixed_project.mlxtran"),
#'   library = "pk",
#'   requiredFilters = list(administration = "bolus", parametrization = "rate"),
#'   settings = list(findIIV = TRUE, findError = TRUE, output = "list_of_runs.csv")
#' )
#'
#' # Example 4: Exhaustive search over oral clearance models
#' res <- findModel(
#'   project = file.path(getDemoPath(), "1.creating_and_using_models",
#'                       "1.1.libraries_of_models", "theophylline_project.mlxtran"),
#'   library = "pk",
#'   requiredFilters = list(administration = "oral", parametrization = "clearance"),
#'   algorithm = "exhaustive_search"
#' )
#'
#' # Example 5: Decision tree search strategy even though searching for IIV
#' res <- findModel(
#'   project = file.path(getDemoPath(), "1.creating_and_using_models",
#'                       "1.1.libraries_of_models", "theophylline_project.mlxtran"),
#'   library = "pk",
#'   requiredFilters = list(administration = "oral", parametrization = "clearance"),
#'   algorithm = "decision_tree",
#'   settings = list(findIIV = TRUE)
#' )
#'
#' # Example 6: Custom model library with exhaustive search
#' # Define a custom model creation function
#' my_custom_func <- function(library, filters) {
#'   # Return path to appropriate model file based on filters
#'   if (filters$absorption == "firstOrder") {
#'     return("/path/to/my/first_order_model.txt")
#'   } else {
#'     return("/path/to/my/zero_order_model.txt")
#'   }
#' }
#'
#' res <- findModel(
#'   project = "my_project.mlxtran",
#'   library = "custom",
#'   filters = list(absorption = c("firstOrder", "zeroOrder")),
#'   algorithm = "exhaustive_search",
#'   model_creation_func = my_custom_func,
#'   settings = list(iiv = FALSE, error = FALSE)
#' )
#' }
#'
#' @importFrom stats aggregate lm optim setNames as.formula model.matrix predict
#' @importFrom utils read.csv write.csv
#' @importFrom cli cli_alert_success cli_progress_bar
#' @importFrom cli cli_progress_done cli_progress_update
#'
#' @export
findModel <- function(project, library, requiredFilters = NULL, filters = NULL,
                      algorithm = NULL, settings = NULL, metric_func = NULL,
                      initial_func = NULL, model_creation_func = NULL,
                      param_mapping = NULL) {

  op = options()
  op$lixoft_notificationOptions$warnings = 1
  op$lixoft_notificationOptions$info = 1
  op$lixoft_notificationOptions$errors = 1
  op$cli.progress_show_after = 0
  options(op)

  if (!is.null(settings)) {
    names(settings)[names(settings) == "findIIV"] <- "iiv"
    names(settings)[names(settings) == "findError"] <- "error"
  }

  # --- Argument Checks ---
  validate_project(project)
  validate_library(library)
  
  # For custom library, requiredFilters can be empty
  if (library != "custom") {
    validate_required_filters(requiredFilters)
  }
  
  if (!is.null(algorithm)) validate_algorithm(algorithm)
  if (!is.null(filters))   validate_filters(filters, requiredFilters, library)

  if (!is.null(metric_func))  validate_metric_func_no_args(metric_func)
  if (!is.null(initial_func)) validate_initial_func_one_arg(initial_func)

  # --- Default Settings ---
  # For custom library, keys are not needed (no default filters)
  if (library != "custom") {
    keys <- switch(requiredFilters$administration,
       "oral" = c("distribution", "absorption_delay", "distribution", "elimination"),
       "oralBolus" = c("bioavailability", "distribution", "absorption_delay", "distribution", "elimination"),
       c("distribution", "elimination")
    )
  } else {
    keys <- NULL
  }

  if (is.null(settings)) settings <- list(iiv = FALSE, error = FALSE)
  if (is.null(settings$iiv)) settings$iiv <- FALSE
  if (is.null(settings$error)) settings$error <- FALSE

  default_settings <- list(
    N = ifelse(is.null(settings$iiv) || !settings$iiv, 10, 20),
    iiv = FALSE,
    error = FALSE,
    seed = 123456,
    previous_runs = NULL,
    output = NULL,
    convergence_threshold = 0.75,
    initial_iiv_forced_iter = 3,
    rho = 0.4,
    save_mode = "best",
    keys = keys,
    plot = TRUE,
    linearization = FALSE,
    elitism = TRUE,
    max_models = 1000000,
    iiv_favor_true = FALSE,
    clip = FALSE,
    xgboost = TRUE,
    xgboost_model_nb = 5,
    stagnation_xgboost = 7,
    stagnation_stop = 8,
    local_search = FALSE,
    table_func = NULL
  )

  # --- Setup ---
  set.seed(settings$seed)
  filters <- updateFilters(library, filters, requiredFilters)
  space <- calculateSearchSpace(library, filters, settings$error, settings$iiv, param_mapping)
  cli::cli_alert_info("Search space size: {space} models")

  # --- Algorithm Selection ---
  if (is.null(algorithm)) {
    algorithm <- if (isTRUE(settings$iiv)) "ACO" else "decision_tree"
    cli::cli_alert_info("Algorithm automatically selected: {algorithm}")
  }

  validate_settings(settings, algorithm)
  
  # Validate custom model setup
  validate_custom_model_setup(library, algorithm, model_creation_func, param_mapping, settings)

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

  init_cov <- lixoftConnectors::getIndividualParameterModel()$covariateModel
  init_cov <- lapply(init_cov, function(x) {
    x[] <- FALSE
    return(x)
  })

  init_var <- lixoftConnectors::getIndividualParameterModel()$variability
  init_var <- lapply(names(init_var), function(nm) {
    x <- init_var[[nm]]
    if (nm == "id") {
      x[] <- TRUE
    } else {
      x[] <- FALSE
    }
    return(x)
  })
  names(init_var) <- names(lixoftConnectors::getIndividualParameterModel()$variability)

  init_cor <- lixoftConnectors::getIndividualParameterModel()$correlationBlocks
  init_cor <- lapply(init_cor, function(x) list())

  lixoftConnectors::setIndividualParameterModel(list(
    correlationBlocks = init_cor,
    covariateModel = init_cov
  ))
  lixoftConnectors::setIndividualParameterModel(list(
    variability = init_var
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
           return(runACO(library, filters, requiredFilters, settings, metric_func, initial_func,
                         model_creation_func, param_mapping))
         },
         tournament = {
           return(runTournament(library, filters, requiredFilters, settings, metric_func, initial_func))
         },
         exhaustive_search = {
           return(runExhaustiveSearch(library, filters, requiredFilters, settings, metric_func, initial_func,
                                      model_creation_func, param_mapping))
         },
         decision_tree = {
           return(runDecisionTree(library, filters, requiredFilters, settings, metric_func, initial_func))
         })
}

