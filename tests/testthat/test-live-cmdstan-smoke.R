test_that("live CmdStan compile and fit smoke test succeeds when configured", {
  skip_if_not(
    live_cmdstan_run_requested(),
    "Live CmdStan smoke is gated by BAYESEFRON_RUN_LIVE=1 or BAYESEFRON_RUN_FULL_LIVE=1."
  )
  live_cmdstan_skip_if_unavailable()

  with_live_cmdstan_cache_root({
    model <- bayes_efron_compile(
      quiet = TRUE,
      force_recompile = TRUE,
      seed_for_check = 47L
    )
    expect_s3_class(model, "CmdStanModel")

    fit <- bayes_efron_fit(
      theta_hat = c(-0.45, -0.1, 0, 0.25, 0.55),
      sigma = c(0.12, 0.18, 0.15, 0.22, 0.2),
      L = 51L,
      M = 3L,
      chains = 1L,
      iter_warmup = 150L,
      iter_sampling = 4L,
      seed = 4701L,
      keep_cmdstan_fit = FALSE
    )

    expect_s3_class(fit, "bef_fit_re")
    expect_s3_class(fit, "bef_fit")
    expect_false("cmdstan_fit" %in% names(fit))
    expect_equal(fit$metadata$model_family, "RE")
    expect_identical(fit$metadata$seed, 4701L)
    expect_equal(fit$metadata$data_list$K, 5L)
    expect_equal(fit$metadata$data_list$L, 51L)
    expect_equal(dim(fit$posterior$theta_rep), c(4L, 5L))
    expect_true(all(is.finite(fit$metadata$theta_summary$mean)))
  })
})
