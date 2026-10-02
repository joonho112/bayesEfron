# Comparison with deconvolveR (Narasimhan and Efron 2020), the reference
# implementation of the empirical Bayes estimator. Its estimates are stored in
# tests/testthat/_fixtures, so deconvolveR itself is not needed to run these
# tests. The fits use the default sampler settings.

parity_fit <- function(panel, seed) {
  # With K = 20 a few divergent transitions can occur for some seeds or
  # platforms. The agreement measured below does not change with them, so the
  # diagnostic warning is not part of what these tests check.
  suppressWarnings(
    bayes_efron_fit(
      theta_hat = panel$theta_hat,
      sigma = panel$sigma,
      L = panel$grid_count,
      M = panel$basis_df,
      expansion = 0.5,
      seed = seed
    ),
    classes = "bef_sampler_diagnostics_failed"
  )
}

test_that("the stored deconvolveR posterior means follow from the stored prior", {
  fixture <- parity_load_fixture()

  for (K in c(20L, 50L, 100L)) {
    panel <- parity_panel(fixture, K)
    expect_equal(
      parity_deconvolveR_theta_mean(panel$theta_hat, panel$tau, panel$g),
      panel$theta_mean_eb,
      tolerance = 1e-12
    )
  }
})

test_that("posterior means agree with deconvolveR within one standard error", {
  skip_if_no_cmdstan()

  fixture <- parity_load_fixture()
  for (K in c(20L, 50L, 100L)) {
    panel <- parity_panel(fixture, K)
    fit <- parity_fit(panel, seed = 26600L + K)

    parity_expect_scaled_agreement(
      coef(fit, type = "mean"),
      panel$theta_mean_eb,
      panel$sigma,
      max_scaled = 1,
      min_cor = 0.99
    )
    expect_equal(fit$metadata$data_list$L, panel$grid_count)
    expect_equal(fit$metadata$data_list$M, panel$basis_df)
  }
})

test_that("the estimated prior agrees with deconvolveR at K = 100", {
  skip_if_no_cmdstan()

  panel <- parity_panel(parity_load_fixture(), 100L)

  fit <- parity_fit(panel, seed = 26700L)
  g_draws <- getFromNamespace(".bef_grid_draws", "bayesEfron")(fit, "g")
  g_fb <- colMeans(as.matrix(g_draws))

  # Two densities on a grid are compared by their total variation distance. A
  # relative comparison is meaningless in the tails, where both are near zero.
  g_eb <- as.numeric(panel$g)
  tv <- 0.5 * sum(abs(g_fb / sum(g_fb) - g_eb / sum(g_eb)))
  expect_lte(tv, 0.25)
})
