# Structure and internal consistency of a fitted object, checked on the short
# five-site fit shipped in inst/examples (one chain, 100 draws; no CmdStan
# needed).

example_fit <- function() {
  readRDS(test_path("_fixtures", "short_fit.rds"))
}

site_draws <- function(fit, name) {
  getFromNamespace(".bef_site_draws", "bayesEfron")(fit, name)
}

test_that("a fit has draws, metadata and posterior with the documented fields", {
  fit <- example_fit()

  expect_s3_class(fit, c("bef_fit_re", "bef_fit"), exact = TRUE)
  expect_named(fit, c("draws", "metadata", "posterior"))
  expect_s3_class(fit$draws, "draws_array")

  expect_named(
    fit$metadata,
    c(
      "model_family", "grid_method", "seed", "cmdstan_version",
      "stan_file_sha256", "data_list", "runtime_seconds", "mean_g_summary",
      "var_g_summary", "theta_summary", "effective_params_summary",
      "log_marginal_likelihood_summary"
    )
  )
  expect_setequal(
    setdiff(names(attributes(fit$metadata)), "names"),
    c(
      "sd_g_summary", "diagnostics", "diagnostic_skipped",
      "sampler_diagnostics_failed", "sampler_diagnostics_warned"
    )
  )
  expect_named(
    fit$metadata$data_list,
    c("K", "theta_hat", "sigma", "L", "grid", "M", "B", "store_grid_quantities")
  )
  expect_match(fit$metadata$stan_file_sha256, "^[0-9a-f]{64}$")
  expect_match(fit$metadata$cmdstan_version, "^[0-9]+[.][0-9]+[.][0-9]+")
})

test_that("posterior holds the five scalar series, one value per draw", {
  fit <- example_fit()
  n_draws <- posterior::ndraws(fit$draws)

  expect_named(
    fit$posterior,
    c("mean_g", "var_g", "sd_g", "effective_params", "log_marginal_likelihood")
  )
  for (series in fit$posterior) {
    expect_type(series, "double")
    expect_length(series, n_draws)
    expect_true(all(is.finite(series)))
  }
  expect_true(all(fit$posterior$var_g >= 0))
  expect_equal(fit$posterior$sd_g, sqrt(fit$posterior$var_g), tolerance = 1e-8)
})

test_that("site-level draws have one column per site and finite values", {
  fit <- example_fit()
  K <- fit$metadata$data_list$K

  for (name in c("theta_mean", "theta_sd", "theta_map", "theta_rep")) {
    draws <- site_draws(fit, name)
    expect_equal(ncol(draws), K)
    expect_true(all(is.finite(draws)))
  }
  expect_true(all(site_draws(fit, "theta_sd") >= 0))

  # sum_i Var(theta_i | data) / sigma_i^2 is bounded below by zero. It is not
  # bounded above by K: under a multimodal g a site between two modes can have
  # a posterior variance larger than sigma_i^2.
  summary_mean <- fit$metadata$effective_params_summary$mean
  expect_true(is.finite(summary_mean))
  expect_gte(summary_mean, 0)
})

test_that("pointwise log-likelihoods sum to the log marginal likelihood", {
  fit <- example_fit()
  draws <- posterior::as_draws_matrix(fit$draws)
  K <- fit$metadata$data_list$K
  log_lik <- draws[, sprintf("log_lik[%d]", seq_len(K)), drop = FALSE]

  expect_equal(
    as.numeric(rowSums(log_lik)),
    as.numeric(draws[, "log_marginal_likelihood"]),
    tolerance = 1e-6
  )
  expect_equal(
    as.numeric(draws[, "log_marginal_likelihood"]),
    fit$posterior$log_marginal_likelihood,
    tolerance = 1e-12
  )
})

test_that("the prior g recovered from log_g is a probability vector in every draw", {
  fit <- example_fit()
  g <- getFromNamespace(".bef_grid_draws", "bayesEfron")(fit, "g")

  expect_equal(ncol(g), fit$metadata$data_list$L)
  expect_true(all(is.finite(g)))
  expect_true(all(g >= 0))
  expect_equal(as.numeric(rowSums(g)), rep(1, nrow(g)), tolerance = 1e-8)
})

test_that("convergence diagnostics are reported per parameter", {
  diagnostic <- diagnose(example_fit())

  expect_gt(length(diagnostic$rhat), 1L)
  expect_false(is.null(names(diagnostic$rhat)))
  expect_identical(names(diagnostic$rhat), names(diagnostic$ess_bulk))
  expect_identical(names(diagnostic$rhat), names(diagnostic$ess_tail))
})

test_that("prepared Stan data keep the number of sites and the basis has full rank", {
  theta_hat <- c(-1, 0, 2, 3)
  sigma <- c(0.2, 0.3, 0.4, 0.5)
  grid <- make_efron_grid(theta_hat, sigma, M = 4L)
  stan_data <- getFromNamespace("prepare_stan_data", "bayesEfron")(
    list(theta_hat = theta_hat, sigma = sigma),
    grid
  )

  expect_identical(stan_data$K, length(theta_hat))
  expect_length(stan_data$theta_hat, stan_data$K)
  expect_length(stan_data$sigma, stan_data$K)
  expect_true(all(is.finite(stan_data$B)))
  # The softmax ignores a constant shift, so the basis must not contain one:
  # with an intercept column added, the matrix has to keep full column rank.
  expect_identical(
    qr(cbind(1, stan_data$B), tol = sqrt(.Machine$double.eps))$rank,
    stan_data$M + 1L
  )
})
