# Empirical Bayes estimates computed with deconvolveR::deconv() (version 1.2.1)
# for three simulated data sets with unit standard errors (K = 20, 50, 100).
parity_load_fixture <- function() {
  readRDS(test_path("_fixtures", "deconvolveR_baseline.rds"))
}

parity_panel <- function(fixture, K) {
  fixture$panels[[sprintf("K%d_homo", as.integer(K))]]
}

# Posterior means of the site effects given a prior g on the grid tau.
parity_deconvolveR_theta_mean <- function(theta_hat, tau, g, sigma = 1) {
  vapply(
    theta_hat,
    function(z) {
      weights <- g * stats::dnorm(z, mean = tau, sd = sigma)
      weights <- weights / sum(weights)
      sum(weights * tau)
    },
    numeric(1)
  )
}

# Agreement between two sets of effect estimates on the scale of the standard
# error: the largest |actual - expected| / sigma, together with the
# correlation. A relative error is not usable here because the reference
# estimates pass through zero. The plug-in and the fully Bayesian estimates are
# different estimators, so this is a check on the implementation and not a
# claim that the two should coincide.
parity_scaled_error <- function(actual, expected, sigma) {
  abs(as.numeric(actual) - as.numeric(expected)) / as.numeric(sigma)
}

parity_expect_scaled_agreement <- function(actual, expected, sigma,
                                           max_scaled = 1.0,
                                           min_cor = 0.99) {
  err <- parity_scaled_error(actual, expected, sigma)
  expect_true(all(is.finite(err)))
  expect_lte(max(err), max_scaled)
  expect_gte(stats::cor(as.numeric(actual), as.numeric(expected)), min_cor)
}
