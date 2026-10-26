#!/usr/bin/env Rscript
# Fast, R-only structural and saved-output checks for the final repository.

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_argument[1]), mustWork = TRUE)
project_dir <- dirname(dirname(script_file))

all_files <- list.files(project_dir, recursive = TRUE, full.names = TRUE, all.files = TRUE)
all_files <- all_files[!grepl("/(\\.git)(/|$)", all_files)]
forbidden <- all_files[grepl("\\.(py|ipynb)$", all_files, ignore.case = TRUE)]
stopifnot(length(forbidden) == 0L)

r_files <- all_files[grepl("\\.R$", all_files)]
invisible(lapply(r_files, parse))

source(file.path(project_dir, "R", "load_modules.R"))
load_blcm_modules(project_dir)
config <- blcm_config()

comparison_file <- file.path(project_dir, "outputs", "comparison_all.csv")
if (file.exists(comparison_file)) {
  comparison <- read.csv(comparison_file)
} else {
  summary_files <- list.files(file.path(project_dir, "outputs"),
    pattern = "^(Herd|Slaughterhouse)_.*_summary\\.csv$",
    full.names = TRUE)
  comparison <- do.call(rbind, lapply(summary_files, read.csv, stringsAsFactors = FALSE))
}

required <- rbind(
  data.frame(model = "Herd", scenario = config$primary_scenarios),
  data.frame(model = "Slaughterhouse",
             scenario = c(config$primary_scenarios, "liver_sp90_independent"))
)
observed <- unique(comparison[c("model", "scenario")])
reported <- do.call(rbind, lapply(config$models, function(model)
  data.frame(model = model, scenario = reported_scenarios(model, config))))
stopifnot(nrow(reported) == 24L,
  all(apply(reported, 1, function(z)
    any(observed$model == z[["model"]] & observed$scenario == z[["scenario"]]))))
stopifnot(
  all(apply(required, 1, function(z)
    any(observed$model == z[["model"]] & observed$scenario == z[["scenario"]]))),
  all(comparison$passes)
)

laboratory_file <- file.path(project_dir, "laboratory_comparison", "posterior_comparison.csv")
if (file.exists(laboratory_file)) {
  laboratory <- read.csv(laboratory_file)
  stopifnot(
    nrow(unique(laboratory[c("model", "scenario", "lab", "scope")])) == 16L,
    all(laboratory$passes)
  )
}

transfer_file <- file.path(project_dir, "outputs", "transferred_test_priors.csv")
if (file.exists(transfer_file)) {
  transfer <- read.csv(transfer_file)
  stopifnot(
    identical(transfer$herd_parameter, c("SeS", "SpS", "SeC", "SpC")),
    all(transfer$beta_alpha > 0),
    all(transfer$beta_beta > 0)
  )
} else {
  transfer <- herd_posterior_test_prior(file.path(project_dir, "outputs"), config)
  stopifnot(all(vapply(transfer, function(x) all(is.finite(x) & x > 0), logical(1))))
}

cat("Repository check passed: all required conditional-independence primary/prior-sensitivity fits are present; covariance fits, if present, are separate structural sensitivity analyses.\n")


