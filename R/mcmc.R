# Chain initialization and posterior summaries.

zero_covariance_position <- function(a, b) {
  lower <- max(-a * b, -(1 - a) * (1 - b))
  upper <- min(a * (1 - b), (1 - a) * b)
  -lower / (upper - lower)
}

initial_values <- function(model, dependent, seed, chain, constrained_liver = FALSE) {
  set.seed(seed)
  se_s <- c(.91, .97, .75, .85)[chain]; sp_s <- c(.94, .86, .98, .90)[chain]
  se_c <- c(.15, .40, .65, .25)[chain]; sp_c <- c(.90, .85, .95, .88)[chain]
  z <- list(SeS = se_s, SpS = sp_s, SeC = se_c, SpC = sp_c,
            .RNG.name = "base::Mersenne-Twister", .RNG.seed = seed)
  if (identical(model, "Herd")) {
    z$p <- c(.08, .015, .65)
    if (dependent) {
      z$uSe <- zero_covariance_position(se_s, se_c)
      z$uSp <- zero_covariance_position(sp_s, sp_c)
    }
  } else {
    z$p <- c(.08, .20, .12, .35)[chain]
    z$SeL <- c(.30, .60, .20, .45)[chain]
    z$SpL <- if (constrained_liver) .96 else c(.93, .85, .97, .90)[chain]
    if (dependent) {
      z$uSeLC <- zero_covariance_position(z$SeL, se_c); z$uSeLS <- zero_covariance_position(z$SeL, se_s)
      z$uSeCS <- zero_covariance_position(se_c, se_s); z$uSpLC <- zero_covariance_position(z$SpL, sp_c)
      z$uSpLS <- zero_covariance_position(z$SpL, sp_s); z$uSpCS <- zero_covariance_position(sp_c, sp_s)
    }
  }
  z
}

summarize_chains <- function(chains, model, scenario) {
  draws <- as.matrix(chains)
  quantiles <- t(apply(draws, 2, quantile, c(.025, .5, .975)))
  ess <- effectiveSize(chains)
  rhat <- gelman.diag(chains, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
  posterior_sd <- apply(draws, 2, sd)
  data.frame(model = model, scenario = scenario, parameter = colnames(draws),
             lower = quantiles[, 1], median = quantiles[, 2], upper = quantiles[, 3],
             mean = colMeans(draws), sd = posterior_sd, ess = as.numeric(ess),
             rhat = as.numeric(rhat), mcse_mean = posterior_sd / sqrt(ess),
             mc_pct_sd = 100 / sqrt(ess), row.names = NULL)
}

