# tests/testthat/helper-custom-models.R

# Create a simple test model file for custom model testing
create_test_model_file <- function() {
  model_content <- '[LONGITUDINAL]
input = {V, Cl, ka}

PK:
compartment(cmt=1, volume=V, concentration=Cc)
oral(cmt=1, ka)
elimination(cmt=1, k=Cl/V)

OUTPUT:
output = {Cc}
'
  
  temp_file <- tempfile(pattern = "test_model_", fileext = ".txt")
  writeLines(model_content, temp_file)
  return(temp_file)
}

# Create a simple custom model function for testing
create_simple_custom_func <- function() {
  # Return a function that creates the model file on-demand
  # This avoids closure issues with parallel processing
  function(library, filters) {
    model_content <- '[LONGITUDINAL]
input = {V, Cl, ka}

PK:
compartment(cmt=1, volume=V, concentration=Cc)
oral(cmt=1, ka)
elimination(cmt=1, k=Cl/V)

OUTPUT:
output = {Cc}
'
    temp_file <- tempfile(pattern = "test_model_", fileext = ".txt")
    writeLines(model_content, temp_file)
    return(temp_file)
  }
}

# Create a custom model function that varies by filter
create_filter_based_custom_func <- function() {
  # Return a function that creates model files on-demand
  # This avoids closure issues with parallel processing
  function(library, filters) {
    if (!is.null(filters$absorption) && filters$absorption == "zeroOrder") {
      model_content <- '[LONGITUDINAL]
input = {V, Cl, Tk0}

PK:
compartment(cmt=1, volume=V, concentration=Cc)
oral(cmt=1, Tk0=Tk0, p=0)
elimination(cmt=1, k=Cl/V)

OUTPUT:
output = {Cc}
'
    } else {
      # Default: firstOrder
      model_content <- '[LONGITUDINAL]
input = {V, Cl, ka}

PK:
compartment(cmt=1, volume=V, concentration=Cc)
oral(cmt=1, ka)
elimination(cmt=1, k=Cl/V)

OUTPUT:
output = {Cc}
'
    }
    
    temp_file <- tempfile(pattern = "test_model_", fileext = ".txt")
    writeLines(model_content, temp_file)
    return(temp_file)
  }
}

# Create param_mapping for simple 1-compartment models
create_simple_param_mapping <- function() {
  list(
    common_params = c("V", "Cl"),
    filter_params = list(
      absorption = list(
        firstOrder = c("ka"),
        zeroOrder = c("Tk0")
      )
    )
  )
}

# Create param_mapping for more complex models
create_complex_param_mapping <- function() {
  list(
    common_params = c("V", "Cl"),
    filter_params = list(
      absorption = list(
        firstOrder = c("ka"),
        zeroOrder = c("Tk0"),
        transitCompartments = c("Ktr", "Mtt")
      ),
      distribution = list(
        "1compartment" = character(0),
        "2compartments" = c("Q", "V2"),
        "3compartments" = c("Q", "V2", "Q2", "V3")
      )
    )
  )
}
