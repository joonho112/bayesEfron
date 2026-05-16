parity_target_ids <- function() {
  c(
    "tier2_deconvolveR_theta_mean_K20_homo",
    "tier2_deconvolveR_theta_mean_K50_homo",
    "tier2_deconvolveR_theta_mean_K100_homo",
    "tier2_deconvolveR_g_mean_K100_homo",
    "tier2_deconvolveR_grid_count",
    "tier2_deconvolveR_basis_df"
  )
}

parity_run_requested <- function() {
  identical(Sys.getenv("BAYESEFRON_RUN_PARITY"), "1")
}

parity_fixture_path <- function(mustWork = TRUE) {
  .bef_target_fixture_path(
    "tier2_deconvolveR_theta_mean_K20_homo",
    mustWork = mustWork
  )
}

parity_load_fixture <- function(path = parity_fixture_path()) {
  readRDS(path)
}

parity_panel_name <- function(K) {
  sprintf("K%d_homo", as.integer(K))
}

parity_panel <- function(fixture, K) {
  fixture$panels[[parity_panel_name(K)]]
}

parity_relative_error <- function(actual, expected) {
  denominator <- pmax(abs(expected), sqrt(.Machine$double.eps))
  abs(actual - expected) / denominator
}

parity_expect_relative_lte <- function(actual, expected, tolerance, info = NULL) {
  error <- parity_relative_error(actual, expected)
  expect_true(
    all(is.finite(error)),
    info = sprintf("%s relative errors must be finite", info %||% "parity")
  )
  expect_lte(
    max(error),
    tolerance
  )
}

parity_validate_fixture <- function(fixture, targets = .bef_load_targets()) {
  expect_type(fixture, "list")
  expect_named(fixture, c("metadata", "targets", "panels"))
  expect_equal(fixture$metadata$fixture_format_version, "v1")
  expect_equal(fixture$metadata$generator, "deconvolveR::deconv")
  expect_equal(fixture$metadata$family, "Normal")
  expect_equal(fixture$metadata$sigma, 1)
  expect_equal(fixture$metadata$grid_count, 101L)
  expect_equal(fixture$metadata$basis_df, 6L)
  expect_equal(names(fixture$panels), parity_panel_name(c(20L, 50L, 100L)))

  target_rows <- targets[targets$target_id %in% parity_target_ids(), , drop = FALSE]
  expect_equal(target_rows$status, rep("deferred", length(parity_target_ids())))
  expect_false(any(target_rows$release_blocking))
  expect_equal(
    unique(target_rows$fixture_path),
    "_fixtures/deconvolveR_baseline.rds"
  )
  expect_setequal(fixture$targets$target_id, parity_target_ids())
  expect_equal(fixture$targets$status, target_rows$status)

  for (K in c(20L, 50L, 100L)) {
    panel <- parity_panel(fixture, K)
    expect_named(
      panel,
      c(
        "K", "seed", "theta_true", "theta_hat", "sigma", "tau",
        "deconv_stats", "g", "theta_mean_eb", "grid_count", "basis_df"
      )
    )
    expect_equal(panel$K, K)
    expect_length(panel$theta_hat, K)
    expect_equal(panel$sigma, rep(1, K))
    expect_length(panel$tau, 101L)
    expect_equal(panel$grid_count, 101L)
    expect_equal(panel$basis_df, 6L)
    expect_equal(nrow(panel$deconv_stats), 101L)
    expect_true(all(c("theta", "g") %in% colnames(panel$deconv_stats)))
    expect_equal(panel$tau, panel$deconv_stats[, "theta"], tolerance = 1e-12)
    expect_equal(panel$g, panel$deconv_stats[, "g"], tolerance = 1e-12)
    expect_equal(sum(panel$g), 1, tolerance = 1e-10)
    expect_length(panel$theta_mean_eb, K)
    expect_true(all(is.finite(panel$theta_mean_eb)))
  }

  invisible(TRUE)
}

parity_targets_active <- function(targets = .bef_load_targets()) {
  target_rows <- targets[targets$target_id %in% parity_target_ids(), , drop = FALSE]
  all(target_rows$status == "active")
}

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

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
