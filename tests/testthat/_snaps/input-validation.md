# errors if required args are missing

    Code
      findModel()
    Condition
      Error in `findModel()`:
      ! argument "project" is missing, with no default

---

    Code
      findModel(project = NULL)
    Condition
      Error:
      ! Error: 'project' must be a single string.

---

    Code
      findModel(project = proj_path)
    Condition
      Error in `findModel()`:
      ! argument "library" is missing, with no default

---

    Code
      findModel(project = proj_path, library = "pk")
    Condition
      Error in `findModel()`:
      ! argument "requiredFilters" is missing, with no default

# errors on unsupported library

    Code
      findModel(project = proj_path, library = "pd", requiredFilters = list(
        administration = "oral", parametrization = "clearance"))
    Condition
      Error:
      ! Error: 'library' must be 'pk'.

# requiredFilters must include administration + parametrization with valid values

    Code
      findModel(project = proj_path, library = "pk", requiredFilters = list())
    Condition
      Error:
      ! Error: 'requiredFilters' must be a list with headers 'administration' and 'parametrization'.

---

    Code
      findModel(project = proj_path, library = "pk", requiredFilters = list(
        administration = "oral"))
    Condition
      Error:
      ! Error: 'requiredFilters' must be a list with headers 'administration' and 'parametrization'.

---

    Code
      findModel(project = proj_path, library = "pk", requiredFilters = list(
        parametrization = "rate"))
    Condition
      Error:
      ! Error: 'requiredFilters' must be a list with headers 'administration' and 'parametrization'.

---

    Code
      findModel(project = proj_path, library = "pk", requiredFilters = list(
        administration = "iv_push", parametrization = "clearance"))
    Condition
      Error:
      ! Error: 'administration' must be one of 'bolus', 'infusion', 'oral', or 'oralBolus'.

---

    Code
      findModel(project = proj_path, library = "pk", requiredFilters = list(
        administration = "oral", parametrization = "foo"))
    Condition
      Error:
      ! Error: 'parametrization' must be 'rate' or 'clearance'.

# optional filters validate known choices and bioavailability rule

    Code
      findModel(project = proj_path, library = "pk", requiredFilters = list(
        administration = "oral", parametrization = "rate"), filters = list(
        absorption = c("firstOrder", "bogus")))
    Condition
      Error:
      ! Error: 'absorption' must contain only values 'zeroOrder', 'firstOrder', 'sigmoid', or 'transitCompartments'.

---

    Code
      findModel(project = proj_path, library = "pk", requiredFilters = list(
        administration = "oral", parametrization = "rate"), filters = list(
        bioavailability = "true"))
    Condition
      Error:
      ! Error: 'bioavailability' filter is only applicable when administration is 'oralBolus'.

# algorithm argument validates choices

    Code
      findModel(project = proj_path, library = "pk", requiredFilters = list(
        administration = "oral", parametrization = "rate"), algorithm = "genetic")
    Condition
      Error:
      ! Error: 'algorithm' must be one of 'ACO', 'tournament', 'exhaustive_search', or 'decision_tree'.

