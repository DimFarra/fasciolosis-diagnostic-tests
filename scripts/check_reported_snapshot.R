#!/usr/bin/env Rscript
# Validate the 24-fit reporting snapshot without fitting models or requiring JAGS.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", arg[1]))))
source(file.path(root, "R", "config.R"))
source(file.path(root, "data_counts.R"))
config <- blcm_config()
manifest <- read.csv(file.path(root, "outputs", "reporting_manifest.csv"))
comparison <- read.csv(file.path(root, "outputs", "comparison_all.csv"))
comparison$passes <- as.logical(comparison$passes)
required <- do.call(rbind, lapply(config$models, function(model)
  data.frame(model = model, scenario = reported_scenarios(model, config))))
key <- function(x) paste(x$model, x$scenario, sep = "/")
stopifnot(nrow(required) == 24L, nrow(manifest) == 24L,
          !anyDuplicated(key(manifest)),
          setequal(key(required), key(manifest)),
          setequal(key(required), key(comparison)),
          all(comparison$passes),
          all(comparison$rhat <= config$max_rhat),
          all(comparison$ess >= config$min_ess),
          sum(ferik_counts) == 217L, sum(slaughter_counts) == 109L)
for (i in seq_len(nrow(manifest))) {
  record <- manifest[i, ]
  summary <- read.csv(file.path(root, record$summary_file))
  summary$passes <- as.logical(summary$passes)
  saved <- comparison[comparison$model == record$model &
                      comparison$scenario == record$scenario, ]
  summary <- summary[order(summary$parameter), ]
  saved <- saved[order(saved$parameter), ]
  rownames(summary) <- rownames(saved) <- NULL
  stopifnot(nrow(summary) == record$rows,
            isTRUE(all.equal(summary, saved, tolerance = 1e-12)),
            record$n == if (record$model == "Herd") 217L else 109L)
}
stopifnot(!"uniform_prevalence_independent" %in% config$primary_scenarios,
          "uniform_prevalence_independent" %in% config$optional_unreported_scenarios,
          "uniform_prevalence_independent" %in% config$scenarios)
cat("Validated 24 reported fits, 222 parameter summaries, dataset sizes and convergence thresholds.\n")

