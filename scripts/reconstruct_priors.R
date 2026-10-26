#!/usr/bin/env Rscript
# Optional PriorGen reconstruction of the reported elicitation targets.

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_argument[1]), mustWork = TRUE)
project_dir <- dirname(dirname(script_file))
old_directory <- setwd(project_dir)
on.exit(setwd(old_directory), add = TRUE)
source("reconstruct_priors.R", local = .GlobalEnv)

