# mlxModelFinder 2.0.0

## New Features

### Custom Model Libraries
- **Custom model support** (`library = "custom"`): use your own model library instead of the built-in PK library. Provide a `model_creation_func(library, filters)` that returns a model file path or a Monolix `lib:` reference. Compatible with `"exhaustive_search"` and `"ACO"` algorithms.
- **Monolix library integration**: `model_creation_func` can return `lib:` references from `getLibraryModelName()`, enabling search over any Monolix library (PKPD, PD, TMDD, TGI, etc.) without writing model files.
- **`param_mapping`** argument: defines which parameters are available for IIV in each model variant via `common_params` and `filter_params`.

### IIV and Parameter Control
- **`fixed_iiv`** argument: lock IIV on (`TRUE`) or off (`FALSE`) for specific parameters, reducing the search space. Works with all algorithms and with both `iiv = TRUE` and `iiv = FALSE`.
- **`param_distributions`** argument: override the default log-normal distribution for specific parameters (supports `"normal"`, `"logNormal"`, `"logitNormal"` with optional limits).

### Parallel Execution
- **Parallel model evaluation** for ACO and exhaustive search using `future` / `future.apply`. Workers run models concurrently across local cores or HPC clusters.

### Search and Algorithm Improvements
- **`settings$tasks`**: configure which Monolix scenario tasks to run (e.g. disable `standardErrorEstimation` for faster exploration).
- **`settings$table_func`**: custom function to output additional information about each run.
