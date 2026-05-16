test_that("Tier 2 deconvolveR fixture policy is registered and default-light", {
  targets <- .bef_load_targets()
  fixture <- parity_load_fixture()

  parity_validate_fixture(fixture, targets = targets)

  for (id in parity_target_ids()) {
    target <- .bef_target(id, targets = targets, statuses = "deferred")
    expect_false(target$release_blocking)
  }
})

test_that("Tier 2 deconvolveR reference posterior means are reproducible from fixture g", {
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

test_that("Tier 2 bayesEfron live parity runs only under BAYESEFRON_RUN_PARITY", {
  skip_if_not(
    parity_run_requested(),
    "Tier 2 deconvolveR parity is gated by BAYESEFRON_RUN_PARITY=1."
  )
  skip_if_not(
    parity_targets_active(),
    "Tier 2 deconvolveR live parity targets remain deferred pending calibration."
  )
  live_cmdstan_skip_if_unavailable()

  fixture <- parity_load_fixture()
  with_live_cmdstan_cache_root({
    for (K in c(20L, 50L, 100L)) {
      panel <- parity_panel(fixture, K)
      target <- .bef_target(
        sprintf("tier2_deconvolveR_theta_mean_K%d_homo", K),
        statuses = "active"
      )

      fit <- bayes_efron_fit(
        theta_hat = panel$theta_hat,
        sigma = panel$sigma,
        L = panel$grid_count,
        M = panel$basis_df,
        expansion = 0.5,
        chains = 1L,
        iter_warmup = 150L,
        iter_sampling = 20L,
        seed = 26600L + K,
        keep_cmdstan_fit = FALSE
      )

      parity_expect_relative_lte(
        coef(fit, type = "mean"),
        panel$theta_mean_eb,
        as.numeric(target$tolerance_value),
        info = sprintf("K=%d theta_mean", K)
      )
      expect_equal(fit$metadata$data_list$L, panel$grid_count)
      expect_equal(fit$metadata$data_list$M, panel$basis_df)
    }
  })
})

test_that("Tier 2 K100 prior-g parity runs only under BAYESEFRON_RUN_PARITY", {
  skip_if_not(
    parity_run_requested(),
    "Tier 2 deconvolveR parity is gated by BAYESEFRON_RUN_PARITY=1."
  )
  skip_if_not(
    parity_targets_active(),
    "Tier 2 deconvolveR live parity targets remain deferred pending calibration."
  )
  live_cmdstan_skip_if_unavailable()

  fixture <- parity_load_fixture()
  panel <- parity_panel(fixture, 100L)
  target <- .bef_target("tier2_deconvolveR_g_mean_K100_homo", statuses = "active")

  with_live_cmdstan_cache_root({
    fit <- bayes_efron_fit(
      theta_hat = panel$theta_hat,
      sigma = panel$sigma,
      L = panel$grid_count,
      M = panel$basis_df,
      expansion = 0.5,
      chains = 1L,
      iter_warmup = 150L,
      iter_sampling = 20L,
      seed = 26700L,
      keep_cmdstan_fit = FALSE
    )

    draw_matrix <- posterior::as_draws_matrix(fit$draws)
    g_columns <- sprintf("g[%d]", seq_len(panel$grid_count))
    skip_if_not(
      all(g_columns %in% colnames(draw_matrix)),
      "CmdStan draws do not include transformed-parameter g columns."
    )

    parity_expect_relative_lte(
      colMeans(as.matrix(draw_matrix[, g_columns, drop = FALSE])),
      panel$g,
      as.numeric(target$tolerance_value),
      info = "K=100 prior g"
    )
  })
})
