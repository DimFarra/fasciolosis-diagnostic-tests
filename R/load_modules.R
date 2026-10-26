# Load the analysis modules into the current R session.

load_blcm_modules <- function(project_dir = ".") {
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  suppressPackageStartupMessages({ library(rjags); library(coda) })
  source(file.path(project_dir, "data_counts.R"), local = .GlobalEnv)
  source(file.path(project_dir, "model_helpers.R"), local = .GlobalEnv)
  source(file.path(project_dir, "R", "config.R"), local = .GlobalEnv)
  source(file.path(project_dir, "R", "priors.R"), local = .GlobalEnv)
  source(file.path(project_dir, "R", "mcmc.R"), local = .GlobalEnv)
  source(file.path(project_dir, "R", "fit.R"), local = .GlobalEnv)
  # Shared names used by the reporting scripts.
  assign("scenario_names", blcm_config()$scenarios, envir = .GlobalEnv)
  assign("init_for", initial_values, envir = .GlobalEnv)
  assign("summarize", summarize_chains, envir = .GlobalEnv)
  invisible(blcm_config())
}

