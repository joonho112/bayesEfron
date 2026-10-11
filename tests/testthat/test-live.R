# Tests that compile the Stan model and sample with CmdStan. They run wherever
# CmdStan is installed and are skipped elsewhere. The fits are deliberately
# short: they check that the pieces fit together, not that the sampler has
# converged, so the convergence warnings of such runs are suppressed.

test_that("the model compiles and a short fit returns a complete object", {
  skip_if_no_cmdstan()

  model <- bayes_efron_compile(quiet = TRUE, force_recompile = TRUE)
  expect_s3_class(model, "CmdStanModel")
  expect_true(file.exists(model$exe_file()))
  expect_match(basename(model$exe_file()), "^efron_mix-[0-9a-f]{12}-cmdstan")
  # The same call compiles the program of version 0.3.0.
  old <- getFromNamespace(".bef_model", "bayesEfron")(program = "efron_re")
  expect_true(file.exists(old$exe_file()))
  expect_match(basename(old$exe_file()), "^efron_re-[0-9a-f]{12}-cmdstan")

  # A second call attaches to the same file without compiling again.
  compiled_at <- file.mtime(model$exe_file())
  again <- bayes_efron_compile()
  expect_identical(again$exe_file(), model$exe_file())
  expect_identical(file.mtime(again$exe_file()), compiled_at)

  fit <- suppressWarnings(bayes_efron_fit(
    theta_hat = c(-0.45, -0.1, 0, 0.25, 0.55),
    sigma = c(0.12, 0.18, 0.15, 0.22, 0.2),
    L = 51L,
    M = 3L,
    chains = 1L,
    iter_warmup = 150L,
    iter_sampling = 4L,
    seed = 4701L,
    keep_cmdstan_fit = FALSE,
    site_draws = TRUE
  ))

  expect_s3_class(fit, c("bef_fit_re", "bef_fit"), exact = TRUE)
  expect_false("cmdstan_fit" %in% names(fit))
  expect_equal(fit$metadata$model_family, "RE")
  expect_identical(fit$metadata$seed, 4701L)
  expect_equal(fit$metadata$data_list$K, 5L)
  expect_equal(fit$metadata$data_list$L, 51L)
  expect_identical(
    fit$metadata$cmdstan_version,
    as.character(cmdstanr::cmdstan_version())
  )
  expect_setequal(
    names(fit$posterior),
    c("mean_g", "var_g", "sd_g", "effective_params", "log_marginal_likelihood")
  )
  expect_equal(
    dim(getFromNamespace(".bef_site_draws", "bayesEfron")(fit, "theta_rep")),
    c(4L, 5L)
  )
  expect_true(all(is.finite(fit$metadata$theta_summary$mean)))
})

test_that("an escalc object goes through as_bef_data(), a short fit and the methods", {
  skip_if_no_cmdstan()
  skip_if_not_installed("metafor")

  es <- metafor::escalc(measure = "GEN", yi = yi, vi = vi, data = raudenbush1985)
  rownames(es) <- make.unique(raudenbush1985$author)
  dat <- as_bef_data(es)
  expect_s3_class(dat, "bef_data")
  expect_identical(dat$source, "metafor::escalc")
  expect_identical(dat$names, make.unique(raudenbush1985$author))
  expect_equal(dat$sigma, raudenbush1985$sei)

  fit <- suppressWarnings(bayes_efron_fit(
    dat$theta_hat, dat$sigma,
    L = 51L, M = 3L, chains = 1L, iter_warmup = 150L, iter_sampling = 50L,
    seed = 1234L
  ))
  K <- nrow(raudenbush1985)
  expect_s3_class(fit, "bef_fit_re")
  expect_false("cmdstan_fit" %in% names(fit))
  expect_equal(fit$metadata$data_list$K, K)
  expect_s3_class(summary(fit), "summary.bef_fit_re")
  expect_s3_class(diagnose(fit), "bef_diagnostic")

  intervals <- confint(fit, level = 0.9, type = "theta")
  expect_named(intervals, c("site", "mean", "lower", "upper", "point"))
  expect_equal(nrow(intervals), K)
  expect_true(all(is.finite(intervals$lower)))
  expect_true(all(intervals$lower <= intervals$upper))

  grDevices::pdf(tempfile(fileext = ".pdf"))
  withr::defer(grDevices::dev.off())
  drawn <- plot(fit, type = "caterpillar")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(drawn, "ggplot")
  } else {
    expect_null(drawn)
  }
})

test_that("interval_coverage() counts the sites whose true effect is covered", {
  draws <- matrix(
    c(
      -1.0, 0.1, 1.0,
      -0.8, 0.2, 1.2,
      -0.6, 0.3, 1.4,
      -0.4, 0.4, 1.6,
      -0.2, 0.5, 1.8
    ),
    nrow = 5L,
    ncol = 3L,
    byrow = TRUE
  )
  theta_true <- c(-0.7, 0.35, 3.0)

  expect_equal(interval_coverage(draws, theta_true, level = 0.8), 2 / 3)
})

test_that("a short fit to one K = 50 benchmark replication gives a coverage in [0, 1]", {
  skip_if_no_cmdstan()

  benchmark <- readRDS(test_path("_fixtures", "lee_sui_K50.rds"))
  replication <- benchmark$replications[[1L]]
  expect_length(replication$theta_hat, 50L)
  expect_length(replication$theta_true, 50L)

  fit <- suppressWarnings(bayes_efron_fit(
    theta_hat = replication$theta_hat,
    sigma = replication$sigma,
    L = replication$L,
    M = replication$M,
    grid_method = replication$grid_method,
    expansion = replication$expansion,
    site_draws = TRUE,
    chains = 1L,
    iter_warmup = 150L,
    iter_sampling = 20L,
    seed = 20260509L
  ))
  theta_rep <- getFromNamespace(".bef_site_draws", "bayesEfron")(fit, "theta_rep")

  expect_equal(dim(theta_rep), c(20L, 50L))
  coverage <- interval_coverage(theta_rep, replication$theta_true)
  expect_true(is.finite(coverage))
  expect_gte(coverage, 0)
  expect_lte(coverage, 1)
})

test_that("a fit keeps the draws of g by default, and the same ones as a fit with site-level draws", {
  skip_if_no_cmdstan()
  call <- function(site_draws) {
    suppressWarnings(bayes_efron_fit(
      theta_hat = c(north = -0.45, south = -0.1, east = 0, west = 0.25, centre = 0.55),
      sigma = c(0.12, 0.18, 0.15, 0.22, 0.2),
      L = 51L, M = 3L, chains = 2L, iter_warmup = 150L, iter_sampling = 40L,
      seed = 4702L, site_draws = site_draws
    ))
  }
  compact <- call(FALSE)
  full <- call(TRUE)
  blocks <- function(fit) unique(sub("\\[.*$", "", dimnames(fit$draws)[[3L]]))
  expect_setequal(
    blocks(compact),
    c("lp__", "alpha", "lambda", "log_g", "mean_g", "var_g", "sd_g",
      "effective_params", "log_marginal_likelihood")
  )
  expect_true(all(c("theta_mean", "theta_sd", "theta_map", "theta_rep", "log_lik") %in%
                    blocks(full)))
  shared <- dimnames(compact$draws)[[3L]]
  expect_identical(compact$draws[, , shared], full$draws[, , shared])
  expect_equal(site_summary(compact), site_summary(full))
  expect_lt(as.numeric(object.size(compact)), as.numeric(object.size(full)))
  # The labels of the sites follow the estimates into the fit.
  expect_named(coef(compact), c("north", "south", "east", "west", "centre"))
  expect_identical(rownames(confint(compact, parm = "east")), "east")
  # Every method works on a fit without site-level draws.
  expect_s3_class(summary(compact), "summary.bef_fit_re")
  expect_s3_class(diagnose(compact), "bef_diagnostic")
  expect_length(predict(compact, n = 5L, seed = 1L), 5L)
  expect_equal(dim(vcov(compact)), c(5L, 5L))
  expect_output(print(compact), "site-level draws not kept")
  if (requireNamespace("loo", quietly = TRUE)) {
    expect_equal(suppressWarnings(loo::loo(compact))$estimates,
                 suppressWarnings(loo::loo(full))$estimates)
    expect_equal(suppressWarnings(loo::waic(compact))$estimates,
                 suppressWarnings(loo::waic(full))$estimates)
    # The log likelihood computed from log_g is the one that Stan wrote.
    expect_equal(
      .bef_log_lik(.bef_engine(full)),
      unname(unclass(.bef_log_lik_draws(full))),
      tolerance = 1e-6, ignore_attr = TRUE
    )
  }
  grDevices::pdf(tempfile(fileext = ".pdf"))
  withr::defer(grDevices::dev.off())
  for (type in c("caterpillar", "g", "diagnostic")) {
    expect_no_error(plot(compact, type = type, backend = "base"))
  }
  # The sampler's quantities of every iteration are kept for bayesplot.
  nuts <- attr(compact$metadata, "nuts")
  expect_equal(dim(nuts)[1:2], c(40L, 2L))
  expect_true(all(c("divergent__", "treedepth__", "energy__", "accept_stat__",
                    "stepsize__", "n_leapfrog__") %in% dimnames(nuts)[[3L]]))
  if (requireNamespace("bayesplot", quietly = TRUE)) {
    np <- bayesplot::nuts_params(compact)
    expect_named(np, c("Chain", "Iteration", "Parameter", "Value"))
    expect_identical(nrow(np), 40L * 2L * dim(nuts)[[3L]])
    expect_s3_class(
      bayesplot::mcmc_nuts_energy(np),
      "ggplot"
    )
  }
})
