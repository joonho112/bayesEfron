post_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

post_theta_hat <- function() {
  c(site_a = -0.45, site_b = -0.1, site_c = 0, site_d = 0.25, site_e = 0.55)
}

post_sigma <- function() {
  c(0.12, 0.18, 0.15, 0.22, 0.2)
}

post_sha <- function(letter = "a") {
  paste(rep(letter, 64L), collapse = "")
}

post_context <- function(keep_cmdstan_fit = FALSE,
                         draws = post_draws_array(),
                         runtime_seconds = 3.25) {
  ctx <- post_ns(".bef_fit_prepare_context")(
    call = quote(bayes_efron_fit(theta_hat, sigma)),
    theta_hat = post_theta_hat(),
    sigma = post_sigma(),
    seed = 321L,
    keep_cmdstan_fit = keep_cmdstan_fit,
    check_installed = FALSE
  )
  ctx$cmdstan_fit <- post_fake_cmdstan_fit(draws)
  ctx$runtime_seconds <- runtime_seconds
  ctx$sampler_settings <- post_ns(".bef_fit_sampler_settings")(
    ctx, interactive_fun = function() FALSE
  )
  ctx
}

post_draws_array <- function(K = 5L, iterations = 2L, chains = 2L, L = 101L) {
  vector_fields <- unlist(
    lapply(
      c("theta_map", "theta_mean", "theta_sd", "theta_rep"),
      function(field) sprintf("%s[%d]", field, seq_len(K))
    ),
    use.names = FALSE
  )
  variables <- c(
    sprintf("log_g[%d]", seq_len(L)),
    "mean_g", "var_g", "sd_g",
    vector_fields,
    "effective_params", "log_marginal_likelihood"
  )
  draws <- array(
    NA_real_,
    dim = c(iterations, chains, length(variables)),
    dimnames = list(NULL, NULL, variables)
  )
  n <- iterations * chains

  fill <- function(variable, values) {
    draws[, , variable] <<- matrix(values, nrow = iterations, ncol = chains)
  }

  # log g for each draw: a normal shape on the grid whose center moves with
  # the draw
  for (point in seq_len(L)) {
    position <- (point - 1) / (L - 1)
    centers <- seq(0.4, 0.6, length.out = n)
    weights <- vapply(
      centers,
      function(center) {
        all_points <- stats::dnorm((seq_len(L) - 1) / (L - 1), center, 0.15)
        stats::dnorm(position, center, 0.15) / sum(all_points)
      },
      numeric(1L)
    )
    fill(sprintf("log_g[%d]", point), log(weights))
  }

  fill("mean_g", seq(0.1, 0.4, length.out = n))
  fill("var_g", seq(0.8, 1.1, length.out = n))
  fill("sd_g", seq(0.9, 1.2, length.out = n))
  fill("effective_params", seq(2.5, 3.5, length.out = n))
  fill("log_marginal_likelihood", seq(-12, -9, length.out = n))

  for (site in seq_len(K)) {
    fill(sprintf("theta_map[%d]", site), rep(-0.3 + site / 10, n))
    fill(sprintf("theta_mean[%d]", site), seq(-0.4, 0.4, length.out = n) + site / 20)
    fill(sprintf("theta_sd[%d]", site), rep(0.08 + site / 100, n))
    fill(sprintf("theta_rep[%d]", site), seq(-0.5, 0.5, length.out = n) + site / 15)
  }

  draws
}

post_fake_cmdstan_fit <- function(draws = post_draws_array(),
                                  record = new.env(parent = emptyenv()),
                                  error = NULL,
                                  diagnostic_summary = list(
                                    num_divergent = 0,
                                    num_max_treedepth = 0,
    ebfmi = 0.9
                                  )) {
  structure(
    list(
      draws = function(format = "draws_array") {
        record$format <- format
        if (!is.null(error)) {
          stop(error, call. = FALSE)
        }
        draws
      },
      diagnostic_summary = function(diagnostics = c("divergences", "treedepth"),
                                    quiet = TRUE) {
        if (is.function(diagnostic_summary)) {
          return(diagnostic_summary(diagnostics = diagnostics, quiet = quiet))
        }
        diagnostic_summary
      }
    ),
    class = "fake_cmdstan_mcmc"
  )
}

post_assemble <- function(context = post_context(),
                          cmdstan_version = "2.38.0",
                          sha = post_sha("a"),
                          postprocess_fun = post_ns("postprocess_stan_draws")) {
  post_ns(".bef_fit_assemble")(
    context,
    postprocess_fun = postprocess_fun,
    cmdstan_version_fun = function() cmdstan_version,
    stan_sha_fun = function(model_family) sha,
    package_version_fun = function(package) {
      switch(package, bayesEfron = "0.3.0", cmdstanr = "0.8.0")
    }
  )
}

expect_post_extraction_failed <- function(err) {
  expect_s3_class(err, "bef_extraction_failed")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
}

test_that("draws are turned into the posterior series and summaries", {
  ctx <- post_context()
  processed <- post_ns("postprocess_stan_draws")(
    cmdstan_fit = ctx$cmdstan_fit,
    stan_data = ctx$stan_data,
    model_family = "RE"
  )

  expect_named(
    processed,
    c(
      "draws", "posterior", "diagnostics", "diagnostic_skipped",
      "sampler_diagnostics_failed",
      "sampler_diagnostics_warned", "mean_g_summary", "var_g_summary",
      "sd_g_summary", "effective_params_summary",
      "log_marginal_likelihood_summary", "theta_summary"
    )
  )
  expect_named(processed$posterior, post_ns(".bef_generated_quantity_fields")())
  expect_named(
    processed$diagnostics,
    c("rhat", "ess_bulk", "ess_tail", "divergences", "max_treedepth", "ebfmi")
  )
  expect_equal(processed$diagnostics$divergences, 0)
  expect_equal(processed$diagnostics$max_treedepth, 0)
  expect_type(processed$diagnostic_skipped, "character")
  expect_type(processed$sampler_diagnostics_failed, "character")
  expect_equal(length(processed$posterior$mean_g), 4L)
  expect_false("theta_mean" %in% names(processed$posterior))
  expect_equal(nrow(processed$theta_summary), 5L)
  expect_named(
    processed$theta_summary,
    c("site", "mean", "sd", "hpdi_lower", "hpdi_upper", "map")
  )
  expect_equal(processed$theta_summary$site, seq_len(5L))
  expect_true(is.finite(processed$mean_g_summary$mean))
  expect_true(is.finite(processed$sd_g_summary$mean))
})

test_that("per-chain sampler counts are added up", {
  ctx <- post_context()
  ctx$cmdstan_fit <- post_fake_cmdstan_fit(
    diagnostic_summary = list(
      num_divergent = c(0, 1, 0, 2),
      num_max_treedepth = c(1, 0, 1, 0)
    )
  )

  expect_warning(
    processed <- post_ns("postprocess_stan_draws")(
      cmdstan_fit = ctx$cmdstan_fit,
      stan_data = ctx$stan_data,
      model_family = "RE"
    ),
    class = "bef_sampler_diagnostics_failed"
  )

  expect_equal(processed$diagnostics$divergences, 3)
  expect_equal(processed$diagnostics$max_treedepth, 2)
  expect_false("sampler_diagnostics" %in% processed$diagnostic_skipped)
  expect_identical(processed$sampler_diagnostics_failed, "divergences")
  expect_true("max_treedepth" %in% processed$sampler_diagnostics_warned)
})

test_that("malformed per-chain sampler counts are recorded as skipped", {
  ctx <- post_context()
  ctx$cmdstan_fit <- post_fake_cmdstan_fit(
    diagnostic_summary = list(
      num_divergent = c(0, NA),
      num_max_treedepth = c(0, 0)
    )
  )

  expect_warning(
    processed <- post_ns("postprocess_stan_draws")(
      cmdstan_fit = ctx$cmdstan_fit,
      stan_data = ctx$stan_data,
      model_family = "RE"
    ),
    class = "bef_diagnostic_skipped"
  )

  expect_true("sampler_diagnostics" %in% processed$diagnostic_skipped)
  expect_true(is.na(processed$diagnostics$divergences))
  expect_true(is.na(processed$diagnostics$max_treedepth))
})

test_that("a failure to read sampler diagnostics is a warning, not an error", {
  ctx <- post_context()
  ctx$cmdstan_fit <- post_fake_cmdstan_fit(
    diagnostic_summary = function(...) {
      stop("diagnostic csv unavailable", call. = FALSE)
    }
  )

  expect_warning(
    processed <- post_ns("postprocess_stan_draws")(
      cmdstan_fit = ctx$cmdstan_fit,
      stan_data = ctx$stan_data,
      model_family = "RE"
    ),
    class = "bef_diagnostic_skipped"
  )

  expect_true("sampler_diagnostics" %in% processed$diagnostic_skipped)
  expect_true(is.na(processed$diagnostics$divergences))
  expect_true(is.na(processed$diagnostics$max_treedepth))
})

test_that("a failure to compute draw diagnostics is a warning, not an error", {
  draws <- post_draws_array()

  expect_warning(
    rhat <- post_ns(".bef_draw_diagnostic")(
      draws,
      "rhat",
      function(...) {
        stop("rhat unavailable", call. = FALSE)
      }
    ),
    class = "bef_diagnostic_skipped"
  )

  expect_true(is.na(rhat))
  expect_equal(
    attr(rhat, "bef_skipped_diagnostic", exact = TRUE),
    "rhat"
  )
})

test_that("diagnostics beyond the failure thresholds give a warning", {
  draws <- post_draws_array(iterations = 250L, chains = 2L)
  ctx <- post_context(draws = draws)
  ctx$cmdstan_fit <- post_fake_cmdstan_fit(
    draws = draws,
    diagnostic_summary = list(num_divergent = 20, num_max_treedepth = 0)
  )

  expect_warning(
    processed <- post_ns("postprocess_stan_draws")(
      cmdstan_fit = ctx$cmdstan_fit,
      stan_data = ctx$stan_data,
      model_family = "RE"
    ),
    class = "bef_sampler_diagnostics_failed"
  )

  expect_true("divergences" %in% processed$sampler_diagnostics_failed)
  expect_equal(processed$diagnostics$divergences, 20)
})

test_that("draw arrays are converted with the posterior package", {
  draws <- post_draws_array()
  draw_matrix <- post_ns(".bef_draws_matrix")(draws)

  expect_s3_class(draw_matrix, "draws_matrix")
  expect_s3_class(draw_matrix, "draws")
  expect_true(is.matrix(draw_matrix))
  expect_equal(dim(draw_matrix), c(4L, dim(draws)[[3L]]))
  expect_equal(colnames(draw_matrix), dimnames(draws)[[3L]])
  expect_equal(as.numeric(draw_matrix[, "mean_g"]), as.numeric(draws[, , "mean_g"]))
})

test_that("the assembled fit is valid and carries no CmdStan handle by default", {
  ctx <- post_context(keep_cmdstan_fit = FALSE)
  fit <- post_assemble(ctx)

  expect_s3_class(fit, "bef_fit_re")
  expect_s3_class(fit, "bef_fit")
  expect_false("cmdstan_fit" %in% names(fit))
  expect_named(fit, c("draws", "metadata", "posterior"))
  expect_equal(names(fit$metadata), post_ns(".bef_fit_re_metadata_fields")())
  expect_equal(fit$metadata$model_family, "RE")
  expect_equal(fit$metadata$grid_method, "paper_realdata")
  expect_identical(fit$metadata$seed, 321L)
  expect_equal(fit$metadata$cmdstan_version, "2.38.0")
  expect_equal(fit$metadata$stan_file_sha256, post_sha("a"))
  expect_identical(fit$metadata$data_list, ctx$stan_data)
  expect_equal(fit$metadata$runtime_seconds, 3.25)
  expect_identical(attr(fit$metadata, "sampler_settings", exact = TRUE),
                   ctx$sampler_settings)
  expect_identical(attr(fit$metadata, "package_version", exact = TRUE), "0.3.0")
  expect_identical(attr(fit$metadata, "cmdstanr_version", exact = TRUE), "0.8.0")
  # The site-level quantities live only in `draws`.
  theta_rep <- post_ns(".bef_site_draws")(fit, "theta_rep")
  expect_equal(dim(theta_rep), c(4L, 5L))
  expect_false(inherits(theta_rep, "draws"))
  expect_false("theta_rep_draws" %in% names(fit$metadata))
  expect_false("diagnostics" %in% names(fit$metadata))
  expect_false("sd_g_summary" %in% names(fit$metadata))
  expect_named(
    attr(fit$metadata, "diagnostics", exact = TRUE),
    c("rhat", "ess_bulk", "ess_tail", "divergences", "max_treedepth", "ebfmi")
  )
  expect_named(attr(fit$metadata, "sd_g_summary", exact = TRUE), c("mean", "sd", "q5", "q50", "q95"))
})

test_that("the CmdStan fit is kept only when requested", {
  ctx <- post_context(keep_cmdstan_fit = TRUE)
  fit <- post_assemble(ctx)

  expect_true("cmdstan_fit" %in% names(fit))
  expect_identical(fit$cmdstan_fit, ctx$cmdstan_fit)
})

test_that("draw extraction errors become bef_extraction_failed", {
  ctx <- post_context()
  ctx$cmdstan_fit <- post_fake_cmdstan_fit(error = "draw extraction failed")

  err <- tryCatch(post_assemble(ctx), error = identity)

  expect_post_extraction_failed(err)
  expect_s3_class(err$parent, "simpleError")
})

test_that("malformed draw arrays and missing generated quantities are rejected", {
  ctx <- post_context(draws = array(1, dim = c(2L, 2L, 2L)))
  err <- tryCatch(post_assemble(ctx), error = identity)
  expect_post_extraction_failed(err)

  draws <- post_draws_array()
  draws <- draws[, , dimnames(draws)[[3L]] != "theta_rep[5]", drop = FALSE]
  ctx <- post_context(draws = draws)
  err <- tryCatch(post_assemble(ctx), error = identity)
  expect_post_extraction_failed(err)
  expect_equal(err$field, "theta_rep")
})

test_that("a missing draws method, missing names, non-finite draws and missing scalars are rejected", {
  ctx <- post_context()
  ctx$cmdstan_fit <- structure(list(), class = "fake_cmdstan_mcmc")
  err <- tryCatch(post_assemble(ctx), error = identity)
  expect_post_extraction_failed(err)

  draws <- post_draws_array()
  dimnames(draws)[[3L]] <- NULL
  ctx <- post_context(draws = draws)
  err <- tryCatch(post_assemble(ctx), error = identity)
  expect_post_extraction_failed(err)

  draws <- post_draws_array()
  draws[1, 1, "mean_g"] <- Inf
  ctx <- post_context(draws = draws)
  err <- tryCatch(post_assemble(ctx), error = identity)
  expect_post_extraction_failed(err)

  draws <- post_draws_array()
  draws <- draws[, , dimnames(draws)[[3L]] != "mean_g", drop = FALSE]
  ctx <- post_context(draws = draws)
  err <- tryCatch(post_assemble(ctx), error = identity)
  expect_post_extraction_failed(err)
  expect_equal(err$field, "mean_g")
})

test_that("other postprocessing errors become bef_extraction_failed", {
  ctx <- post_context()
  err <- tryCatch(
    post_assemble(
      ctx,
      postprocess_fun = function(...) {
        stop("postprocess bug", call. = FALSE)
      }
    ),
    error = identity
  )

  expect_post_extraction_failed(err)
  expect_s3_class(err$parent, "simpleError")
})

test_that("validation failures are bef_invalid_fit errors", {
  ctx <- post_context()
  err <- tryCatch(post_assemble(ctx, sha = "not-a-sha"), error = identity)

  expect_s3_class(err, "bef_invalid_fit")
  expect_s3_class(err, "bef_validate_error")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
})

test_that("validation failures in posterior and metadata are bef_invalid_fit errors", {
  ctx <- post_context()

  negative_sd <- function(...) {
    processed <- post_ns("postprocess_stan_draws")(...)
    # `posterior` keeps only the scalar series, so the malformed case is a
    # non-finite scalar.
    processed$posterior$sd_g[1] <- NA_real_
    processed
  }
  err <- tryCatch(post_assemble(ctx, postprocess_fun = negative_sd), error = identity)
  expect_s3_class(err, "bef_invalid_fit")
  expect_s3_class(err, "bef_validate_error")

  bad_theta_summary <- function(...) {
    processed <- post_ns("postprocess_stan_draws")(...)
    processed$theta_summary$sd[1] <- -0.1
    processed
  }
  err <- tryCatch(
    post_assemble(ctx, postprocess_fun = bad_theta_summary),
    error = identity
  )
  expect_s3_class(err, "bef_invalid_fit")
  expect_s3_class(err, "bef_validate_error")

  bad_metadata_attr <- function(...) {
    processed <- post_ns("postprocess_stan_draws")(...)
    processed$diagnostics <- NULL
    processed
  }
  err <- tryCatch(
    post_assemble(ctx, postprocess_fun = bad_metadata_attr),
    error = identity
  )
  expect_s3_class(err, "bef_invalid_fit")
  expect_s3_class(err, "bef_validate_error")

  # A theta_summary that disagrees with the draws it claims to summarize.
  bad_theta_summary_mean <- function(...) {
    processed <- post_ns("postprocess_stan_draws")(...)
    processed$theta_summary$mean[1] <- processed$theta_summary$mean[1] + 1
    processed
  }
  err <- tryCatch(
    post_assemble(ctx, postprocess_fun = bad_theta_summary_mean),
    error = identity
  )
  expect_s3_class(err, "bef_invalid_fit")
  expect_s3_class(err, "bef_validate_error")
})

test_that("the Stan file hash is the SHA-256 of the installed Stan file", {
  expected <- digest::digest(
    post_ns(".bef_stan_file")("RE"),
    algo = "sha256",
    file = TRUE
  )

  expect_equal(post_ns(".bef_fit_stan_sha256")("RE"), expected)
})
