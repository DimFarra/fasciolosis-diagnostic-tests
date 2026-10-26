#!/usr/bin/env Rscript
# GitHub entry point. Examples:
# Rscript scripts/run_models.R
# Rscript scripts/run_models.R --model Herd --scenario baseline_independent
# Rscript scripts/run_models.R --include-dependence  # optional structural sensitivity

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_argument[1]), mustWork = TRUE)
project_dir <- dirname(dirname(script_file))
source(file.path(project_dir, "R", "load_modules.R")); load_blcm_modules(project_dir)
args <- commandArgs(TRUE)
value_after <- function(flag) {
  index <- match(flag, args)
  if (is.na(index)) return(NULL)
  if (index == length(args)) stop("Missing value after ", flag)
  args[index + 1]
}
config <- blcm_config(); requested_model <- value_after("--model"); requested_scenario <- value_after("--scenario")
models <- if (is.null(requested_model)) config$models else requested_model
stopifnot(all(models %in% config$models)); if (!is.null(requested_scenario)) stopifnot(requested_scenario %in% config$scenarios)
include_dependence <- "--include-dependence" %in% args
failed <- FALSE
for (model in models) {
  if (is.null(requested_scenario)) {
    scenarios <- config$primary_scenarios
    if (identical(model, "Slaughterhouse")) scenarios <- c(scenarios, "liver_sp90_independent")
    if (include_dependence) {
      scenarios <- c(scenarios, config$dependence_scenarios)
      if (identical(model, "Slaughterhouse")) scenarios <- c(scenarios, "liver_sp90")
    }
  } else {
    scenarios <- requested_scenario
  }
  scenarios <- unique(scenarios[scenarios %in% valid_scenarios(model, config)])
  if (!length(scenarios)) next
  for (scenario in scenarios) tryCatch(
    run_one(model, scenario, file.path(project_dir, "outputs"), file.path(project_dir, "models"), config),
    error = function(e) {
      failed <<- TRUE; message("FAILED ", model, " ", scenario, ": ", conditionMessage(e))
      writeLines(conditionMessage(e), file.path(project_dir, "outputs", paste0(model, "_", scenario, "_ERROR.txt")))
    }
  )
}
writeLines(capture.output(sessionInfo()), file.path(project_dir, "outputs", "sessionInfo.txt"))
if (failed) quit(status = 1)

