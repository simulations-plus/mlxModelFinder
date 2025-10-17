# Fail fast if demo project is missing in the installed pkg
demo_proj <- system.file("extdata", "theophylline_project.mlxtran",
                         package = "mlxModelFinder")

mlx_available <- function() {
  if (!file.exists(demo_proj)) return(FALSE)
  if (!requireNamespace("lixoftConnectors", quietly = TRUE)) return(FALSE)
  TRUE
}

skip_if_no_monolix <- function() {
  skip_on_cran()
  if (!mlx_available()) skip("Monolix or demo project not available.")
}

# Small, deterministic settings for integration tests
tiny_settings <- function(overrides = list()) {
  base <- list(
    plot = FALSE, seed = 123456, N = 1, max_models = 1, save_mode = "none"
  )
  utils::modifyList(base, overrides, keep.null = TRUE)
}

# Helper: make a tiny search space so “exhaustive” is still quick
tiny_pk_filters <- function() {
  list(
    absorption   = c("zeroOrder", "firstOrder"),
    distribution = "1compartment",
    elimination  = "linear",
    delay        = "noDelay"
  )
}
