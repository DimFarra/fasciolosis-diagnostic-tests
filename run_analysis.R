# Entry point for the Bayesian latent class analysis.
# Preferred interface: Rscript scripts/run_models.R --model Herd --scenario baseline_independent

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
project_dir <- if (sys.nframe() == 0 && length(script_argument)) {
  dirname(normalizePath(sub("^--file=", "", script_argument[1]), mustWork = TRUE))
} else {
  normalizePath(".", mustWork = TRUE)
}
source(file.path(project_dir, "R", "load_modules.R")); load_blcm_modules(project_dir)

if (sys.nframe() == 0) {
  args <- commandArgs(TRUE)
  config <- blcm_config()
  models <- config$models
  include_dependence <- "--include-dependence" %in% args
  failed <- FALSE
  for (model in models) {
    scenarios <- config$primary_scenarios
    if (identical(model, "Slaughterhouse")) scenarios <- c(scenarios, "liver_sp90_independent")
    if (include_dependence) {
      scenarios <- c(scenarios, config$dependence_scenarios)
      if (identical(model, "Slaughterhouse")) scenarios <- c(scenarios, "liver_sp90")
    }
    scenarios <- unique(scenarios[scenarios %in% valid_scenarios(model, config)])
    for (scenario in scenarios) {
      tryCatch(
        run_one(model, scenario, file.path(project_dir, "outputs"), file.path(project_dir, "models")),
        error = function(e) {
          failed <<- TRUE
          message("FAILED ", model, " ", scenario, ": ", conditionMessage(e))
          writeLines(conditionMessage(e), file.path(project_dir, "outputs", paste0(model, "_", scenario, "_ERROR.txt")))
        }
      )
    }
  }
  writeLines(capture.output(sessionInfo()), file.path(project_dir, "outputs", "sessionInfo.txt"))
  if (failed) quit(status = 1)
}

