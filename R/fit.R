# Construction, caching and fitting of one model/scenario combination.

model_seed <- function(model, scenario, config = blcm_config()) {
  offset <- if (identical(model, "Herd")) 1L else 2L
  config$base_seed + 1000L * offset + 10L * match(scenario, config$scenarios)
}

model_components <- function(model, scenario, output_dir = "outputs", config = blcm_config()) {
  dependent <- !grepl("independent", scenario)
  constrained_liver <- grepl("liver_sp90", scenario)
  prior <- make_prior(model, scenario, output_dir, config)
  text <- if (identical(model, "Herd")) ferik_model_string(dependent) else slaughter_model_string(dependent, constrained_liver)
  text <- sub("model {", "model {\n orientation <- step(SeS + SpS - 1 - 1.0E-8)\n oneOrient ~ dbern(orientation)", text, fixed = TRUE)
  data <- if (identical(model, "Herd")) make_ferik_data(prior) else make_slaughter_data(prior)
  if (!dependent) { data$oneD <- NULL; data$oneN <- NULL }
  data$oneOrient <- 1L
  monitor <- if (identical(model, "Herd")) ferik_monitor(dependent) else slaughter_monitor(dependent)
  list(prior = prior, text = text, data = data, monitor = monitor,
       dependent = dependent, constrained_liver = constrained_liver)
}

run_one <- function(model, scenario, output_dir = "outputs", model_dir = "models", config = blcm_config()) {
  stopifnot(model %in% config$models, scenario %in% valid_scenarios(model, config))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE); dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(output_dir, paste0(model, "_", scenario)); model_file <- file.path(model_dir, paste0(model, "_", scenario, ".jags"))
  cat(format(Sys.time()), model, scenario, "START\n"); flush.console()
  component <- model_components(model, scenario, output_dir, config)

  # Reuse a cached fit only when its exact data, priors and JAGS model match.
  cache_ok <- FALSE
  if (file.exists(paste0(prefix, "_summary.csv")) &&
      file.exists(paste0(prefix, "_inputs.rds")) &&
      file.exists(model_file)) {
    old <- readRDS(paste0(prefix, "_inputs.rds"))
    old_text <- paste(readLines(model_file), collapse = "\n")
    cache_ok <- identical(old$data, component$data) &&
                identical(old$prior, component$prior) &&
                identical(trimws(old_text), trimws(component$text)) &&
                file.exists(paste0(prefix, "_draws.rds"))
  }
  if (cache_ok) {
    cached <- read.csv(paste0(prefix, "_summary.csv"))
    stopifnot(all(cached$passes))
    return(cached)
  }

  writeLines(component$text, model_file); saveRDS(list(data = component$data, prior = component$prior), paste0(prefix, "_inputs.rds"))
  seed <- model_seed(model, scenario, config)
  jags <- jags.model(textConnection(component$text), data = component$data,
                     inits = lapply(seq_len(config$chains), function(chain) initial_values(model, component$dependent, seed + chain, chain, component$constrained_liver)),
                     n.chains = config$chains, n.adapt = config$adapt, quiet = TRUE)
  update(jags, config$burnin, progress.bar = "none")
  chains <- coda.samples(jags, component$monitor, n.iter = config$sample_block, thin = config$thin, progress.bar = "none")
  iterations <- config$sample_block
  repeat {
    summary <- summarize_chains(chains, model, scenario)
    if (all(summary$rhat <= config$max_rhat & summary$ess >= config$min_ess) || iterations >= config$max_sample) break
    additional <- coda.samples(jags, component$monitor, n.iter = config$sample_block, thin = config$thin, progress.bar = "none")
    chains <- mcmc.list(lapply(seq_len(config$chains), function(chain) mcmc(rbind(as.matrix(chains[[chain]]), as.matrix(additional[[chain]])), start = 1, thin = config$thin)))
    iterations <- iterations + config$sample_block
  }
  summary$iterations_per_chain <- iterations; summary$passes <- summary$rhat <= config$max_rhat & summary$ess >= config$min_ess
  saveRDS(chains, paste0(prefix, "_draws.rds")); write.csv(summary, paste0(prefix, "_summary.csv"), row.names = FALSE)
  stopifnot(all(summary$passes))
  cat(format(Sys.time()), model, scenario, "DONE", iterations, "maxRhat", max(summary$rhat), "minESS", min(summary$ess), "\n"); flush.console()
  summary
}

