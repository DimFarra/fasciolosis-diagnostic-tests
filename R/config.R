# Central configuration for the final BLCM workflow.

blcm_config <- function() {
  list(
    models = c("Herd", "Slaughterhouse"),
    scenarios = c(
      "baseline", "widened", "uniform_all", "beta22_tests",
      "uniform_tests", "uniform_prevalence", "baseline_independent",
      "widened_independent", "uniform_all_independent",
      "beta22_tests_independent", "uniform_tests_independent", "uniform_prevalence_independent",
      "liver_sp90", "liver_sp90_independent"
    ),
    herd_transfer_scenario = "baseline_independent",
    primary_scenarios = c(
      "baseline_independent", "widened_independent",
      "beta22_tests_independent", "uniform_tests_independent",
      "uniform_all_independent"
    ),
    optional_unreported_scenarios = "uniform_prevalence_independent",
    dependence_scenarios = c("baseline", "widened", "beta22_tests", "uniform_tests", "uniform_all", "uniform_prevalence"),
    chains = 4L, adapt = 5000L, burnin = 10000L,
    sample_block = 100000L, max_sample = 2000000L, thin = 5L,
    max_rhat = 1.01, min_ess = 1500, base_seed = 910260L
  )
}

valid_scenarios <- function(model, config = blcm_config()) {
  scenarios <- config$scenarios
  if (identical(model, "Herd")) scenarios <- scenarios[!grepl("liver", scenarios)]
  scenarios
}


# Manuscript snapshot: 11 herd fits and 13 slaughterhouse fits.
reported_scenarios <- function(model, config = blcm_config()) {
  setdiff(valid_scenarios(model, config), config$optional_unreported_scenarios)
}

