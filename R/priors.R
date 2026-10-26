# Prior construction and transformations used by the sensitivity scenarios.

beta_from_moments <- function(mean, sd) {
  stopifnot(length(mean) == 1L, length(sd) == 1L,
            is.finite(mean), is.finite(sd), mean > 0, mean < 1, sd > 0)
  variance <- sd^2
  max_variance <- mean * (1 - mean)
  if (variance >= max_variance) {
    stop("Posterior moments cannot be represented by a Beta distribution: variance too large.")
  }
  concentration <- max_variance / variance - 1
  c(mean * concentration, (1 - mean) * concentration)
}

herd_posterior_test_prior <- function(output_dir = "outputs", config = blcm_config()) {
  scenario <- config$herd_transfer_scenario
  summary_file <- file.path(output_dir, paste0("Herd_", scenario, "_summary.csv"))
  if (!file.exists(summary_file)) {
    stop(
      "Sequential prior transfer requires the herd model to be fitted first. Missing: ",
      summary_file,
      ". Run the Herd ", scenario, " model before fitting slaughterhouse baseline scenarios."
    )
  }

  herd <- read.csv(summary_file, stringsAsFactors = FALSE)
  required <- c("SeS", "SpS", "SeC", "SpC")
  if (!all(required %in% herd$parameter)) {
    stop("Herd summary does not contain all four test-performance parameters required for transfer.")
  }

  get_beta <- function(parameter) {
    row <- herd[herd$parameter == parameter, , drop = FALSE]
    stopifnot(nrow(row) == 1L, isTRUE(row$passes))
    beta_from_moments(row$mean, row$sd)
  }

  list(
    se_sero = get_beta("SeS"),
    sp_sero = get_beta("SpS"),
    se_copro = get_beta("SeC"),
    sp_copro = get_beta("SpC")
  )
}

make_prior <- function(model, scenario, output_dir = "outputs", config = blcm_config()) {
  if (identical(model, "Herd")) {
    prior <- ferik_prior
  } else {
    prior <- slaughter_prior

    # For baseline-derived slaughterhouse scenarios, transfer the posterior
    # test-performance information from the primary herd model by fitting
    # Beta distributions to the herd posterior mean and variance.
    uses_transferred_test_priors <- !grepl("uniform_all|beta22_tests|uniform_tests", scenario)
    if (uses_transferred_test_priors) {
      transferred <- herd_posterior_test_prior(output_dir, config)
      prior[names(transferred)] <- transferred
    }
  }

  is_prevalence <- grepl("^p($|_)", names(prior))
  is_test <- !is_prevalence

  if (grepl("widened", scenario)) {
    prior[is_test] <- lapply(prior[is_test], function(x) if (identical(x, c(1, 1))) x else x / 5)
    prior[is_prevalence] <- lapply(prior[is_prevalence], function(x) if (identical(x, c(1, 1))) x else x / 5)
  }
  if (grepl("uniform_all", scenario)) prior <- lapply(prior, function(x) c(1, 1))
  if (grepl("beta22_tests", scenario)) prior[is_test] <- lapply(prior[is_test], function(x) c(2, 2))
  if (grepl("uniform_tests", scenario)) prior[is_test] <- lapply(prior[is_test], function(x) c(1, 1))
  if (grepl("^uniform_prevalence", scenario)) prior[is_prevalence] <- lapply(prior[is_prevalence], function(x) c(1, 1))
  prior
}

