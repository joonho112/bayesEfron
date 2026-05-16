fit8_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

fit8_theta_hat <- function() {
  c(site_a = -0.45, site_b = -0.1, site_c = 0, site_d = 0.25, site_e = 0.55)
}

fit8_sigma <- function() {
  c(0.12, 0.18, 0.15, 0.22, 0.2)
}

fit8_sha <- function(letter = "a") {
  paste(rep(letter, 64L), collapse = "")
}

fit8_context <- function(keep_cmdstan_fit = FALSE,
                         draws = fit8_draws_array(),
                         runtime_seconds = 3.25) {
  ctx <- fit8_ns(".bef_fit_prepare_context")(
    call = quote(bayes_efron_fit(theta_hat, sigma)),
    theta_hat = fit8_theta_hat(),
    sigma = fit8_sigma(),
    seed = 321L,
    keep_cmdstan_fit = keep_cmdstan_fit,
    check_installed = FALSE
  )
  ctx$cmdstan_fit <- fit8_fake_cmdstan_fit(draws)
  ctx$runtime_seconds <- runtime_seconds
  ctx
}

fit8_draws_array <- function(K = 5L, iterations = 2L, chains = 2L) {
  vector_fields <- unlist(
    lapply(
      c("theta_map", "theta_mean", "theta_sd", "theta_rep"),
      function(field) sprintf("%s[%d]", field, seq_len(K))
    ),
    use.names = FALSE
  )
  variables <- c(
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

fit8_fake_cmdstan_fit <- function(draws = fit8_draws_array(),
                                  record = new.env(parent = emptyenv()),
                                  error = NULL,
                                  diagnostic_summary = list(
                                    num_divergent = 0,
                                    num_max_treedepth = 0
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

fit8_assemble <- function(context = fit8_context(),
                          cmdstan_version = "2.38.0",
                          sha = fit8_sha("a"),
                          postprocess_fun = fit8_ns("postprocess_stan_draws")) {
  fit8_ns(".bef_fit_assemble")(
    context,
    postprocess_fun = postprocess_fun,
    cmdstan_version_fun = function() cmdstan_version,
    stan_sha_fun = function(model_family) sha
  )
}

expect_fit8_extraction_failed <- function(err) {
  expect_s3_class(err, "bef_extraction_failed")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
  expect_equal(err$stage, 8L)
}

test_that("Stage 8 postprocesses draws into the RE posterior and summaries", {
  ctx <- fit8_context()
  processed <- fit8_ns("postprocess_stan_draws")(
    cmdstan_fit = ctx$cmdstan_fit,
    stan_data = ctx$stan_data,
    model_family = "RE"
  )

  expect_named(
    processed,
    c(
      "draws", "posterior", "diagnostics", "diagnostic_skipped",
      "sampler_diagnostics_failed", "mean_g_summary", "var_g_summary",
      "sd_g_summary", "effective_params_summary",
      "log_marginal_likelihood_summary", "theta_summary", "theta_rep_draws"
    )
  )
  expect_named(processed$posterior, fit8_ns(".bef_generated_quantity_fields")())
  expect_named(
    processed$diagnostics,
    c("rhat", "ess_bulk", "ess_tail", "divergences", "max_treedepth")
  )
  expect_equal(processed$diagnostics$divergences, 0)
  expect_equal(processed$diagnostics$max_treedepth, 0)
  expect_type(processed$diagnostic_skipped, "character")
  expect_type(processed$sampler_diagnostics_failed, "character")
  expect_equal(dim(processed$posterior$theta_mean), c(4L, 5L))
  expect_equal(dim(processed$theta_rep_draws), c(4L, 5L))
  expect_equal(nrow(processed$theta_summary), 5L)
  expect_named(
    processed$theta_summary,
    c("site", "mean", "sd", "hpdi_lower", "hpdi_upper", "map")
  )
  expect_equal(processed$theta_summary$site, seq_len(5L))
  expect_true(is.finite(processed$mean_g_summary$mean))
  expect_true(is.finite(processed$sd_g_summary$mean))
})

test_that("Stage 8 aggregates per-chain sampler diagnostic counts", {
  ctx <- fit8_context()
  ctx$cmdstan_fit <- fit8_fake_cmdstan_fit(
    diagnostic_summary = list(
      num_divergent = c(0, 1, 0, 2),
      num_max_treedepth = c(1, 0, 1, 0)
    )
  )

  processed <- fit8_ns("postprocess_stan_draws")(
    cmdstan_fit = ctx$cmdstan_fit,
    stan_data = ctx$stan_data,
    model_family = "RE"
  )

  expect_equal(processed$diagnostics$divergences, 3)
  expect_equal(processed$diagnostics$max_treedepth, 2)
  expect_false("sampler_diagnostics" %in% processed$diagnostic_skipped)
})

test_that("Stage 8 records malformed per-chain sampler diagnostic counts as skipped", {
  ctx <- fit8_context()
  ctx$cmdstan_fit <- fit8_fake_cmdstan_fit(
    diagnostic_summary = list(
      num_divergent = c(0, NA),
      num_max_treedepth = c(0, 0)
    )
  )

  expect_warning(
    processed <- fit8_ns("postprocess_stan_draws")(
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

test_that("Stage 8 diagnostics downgrade sampler diagnostic extraction failures", {
  ctx <- fit8_context()
  ctx$cmdstan_fit <- fit8_fake_cmdstan_fit(
    diagnostic_summary = function(...) {
      stop("diagnostic csv unavailable", call. = FALSE)
    }
  )

  expect_warning(
    processed <- fit8_ns("postprocess_stan_draws")(
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

test_that("Stage 8 diagnostics downgrade draw diagnostic failures", {
  draws <- fit8_draws_array()

  expect_warning(
    rhat <- fit8_ns(".bef_draw_diagnostic")(
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

test_that("Stage 8 diagnostics warn when sampler health thresholds fail", {
  draws <- fit8_draws_array(iterations = 250L, chains = 2L)
  ctx <- fit8_context(draws = draws)
  ctx$cmdstan_fit <- fit8_fake_cmdstan_fit(
    draws = draws,
    diagnostic_summary = list(num_divergent = 20, num_max_treedepth = 0)
  )

  expect_warning(
    processed <- fit8_ns("postprocess_stan_draws")(
      cmdstan_fit = ctx$cmdstan_fit,
      stan_data = ctx$stan_data,
      model_family = "RE"
    ),
    class = "bef_sampler_diagnostics_failed"
  )

  expect_true("divergences" %in% processed$sampler_diagnostics_failed)
  expect_equal(processed$diagnostics$divergences, 20)
})

test_that("Stage 8 normalizes draw arrays through posterior", {
  draws <- fit8_draws_array()
  draw_matrix <- fit8_ns(".bef_draws_matrix")(draws)

  expect_s3_class(draw_matrix, "draws_matrix")
  expect_s3_class(draw_matrix, "draws")
  expect_true(is.matrix(draw_matrix))
  expect_equal(dim(draw_matrix), c(4L, dim(draws)[[3L]]))
  expect_equal(colnames(draw_matrix), dimnames(draws)[[3L]])
  expect_equal(as.numeric(draw_matrix[, "mean_g"]), as.numeric(draws[, , "mean_g"]))
})

test_that("Stage 8 assembles and validates a portable bef_fit_re by default", {
  ctx <- fit8_context(keep_cmdstan_fit = FALSE)
  fit <- fit8_assemble(ctx)

  expect_s3_class(fit, "bef_fit_re")
  expect_s3_class(fit, "bef_fit")
  expect_false("cmdstan_fit" %in% names(fit))
  expect_named(fit, c("draws", "metadata", "posterior"))
  expect_equal(names(fit$metadata), fit8_ns(".bef_fit_re_metadata_fields")())
  expect_equal(fit$metadata$model_family, "RE")
  expect_equal(fit$metadata$grid_method, "paper_realdata")
  expect_identical(fit$metadata$seed, 321L)
  expect_equal(fit$metadata$cmdstan_version, "2.38.0")
  expect_equal(fit$metadata$stan_file_sha256, fit8_sha("a"))
  expect_identical(fit$metadata$data_list, ctx$stan_data)
  expect_equal(fit$metadata$runtime_seconds, 3.25)
  expect_equal(dim(fit$posterior$theta_rep), c(4L, 5L))
  expect_false(inherits(fit$posterior$theta_rep, "draws"))
  expect_false(inherits(fit$metadata$theta_rep_draws, "draws"))
  expect_false("diagnostics" %in% names(fit$metadata))
  expect_false("sd_g_summary" %in% names(fit$metadata))
  expect_named(
    attr(fit$metadata, "diagnostics", exact = TRUE),
    c("rhat", "ess_bulk", "ess_tail", "divergences", "max_treedepth")
  )
  expect_named(attr(fit$metadata, "sd_g_summary", exact = TRUE), c("mean", "sd", "q5", "q50", "q95"))
})

test_that("Stage 8 retains raw CmdStan fit only when requested", {
  ctx <- fit8_context(keep_cmdstan_fit = TRUE)
  fit <- fit8_assemble(ctx)

  expect_true("cmdstan_fit" %in% names(fit))
  expect_identical(fit$cmdstan_fit, ctx$cmdstan_fit)
})

test_that("Stage 8 wraps draw extraction failures as bef_extraction_failed", {
  ctx <- fit8_context()
  ctx$cmdstan_fit <- fit8_fake_cmdstan_fit(error = "draw extraction failed")

  err <- tryCatch(fit8_assemble(ctx), error = identity)

  expect_fit8_extraction_failed(err)
  expect_s3_class(err$parent, "simpleError")
})

test_that("Stage 8 rejects malformed draw arrays and missing generated quantities", {
  ctx <- fit8_context(draws = array(1, dim = c(2L, 2L, 2L)))
  err <- tryCatch(fit8_assemble(ctx), error = identity)
  expect_fit8_extraction_failed(err)

  draws <- fit8_draws_array()
  draws <- draws[, , dimnames(draws)[[3L]] != "theta_rep[5]", drop = FALSE]
  ctx <- fit8_context(draws = draws)
  err <- tryCatch(fit8_assemble(ctx), error = identity)
  expect_fit8_extraction_failed(err)
  expect_equal(err$field, "theta_rep")
})

test_that("Stage 8 rejects missing draw method, missing names, nonfinite draws, and scalar omissions", {
  ctx <- fit8_context()
  ctx$cmdstan_fit <- structure(list(), class = "fake_cmdstan_mcmc")
  err <- tryCatch(fit8_assemble(ctx), error = identity)
  expect_fit8_extraction_failed(err)

  draws <- fit8_draws_array()
  dimnames(draws)[[3L]] <- NULL
  ctx <- fit8_context(draws = draws)
  err <- tryCatch(fit8_assemble(ctx), error = identity)
  expect_fit8_extraction_failed(err)

  draws <- fit8_draws_array()
  draws[1, 1, "mean_g"] <- Inf
  ctx <- fit8_context(draws = draws)
  err <- tryCatch(fit8_assemble(ctx), error = identity)
  expect_fit8_extraction_failed(err)

  draws <- fit8_draws_array()
  draws <- draws[, , dimnames(draws)[[3L]] != "mean_g", drop = FALSE]
  ctx <- fit8_context(draws = draws)
  err <- tryCatch(fit8_assemble(ctx), error = identity)
  expect_fit8_extraction_failed(err)
  expect_equal(err$field, "mean_g")
})

test_that("Stage 8 wraps untyped postprocess failures as extraction failures", {
  ctx <- fit8_context()
  err <- tryCatch(
    fit8_assemble(
      ctx,
      postprocess_fun = function(...) {
        stop("postprocess bug", call. = FALSE)
      }
    ),
    error = identity
  )

  expect_fit8_extraction_failed(err)
  expect_s3_class(err$parent, "simpleError")
})

test_that("Stage 8 validation failures remain typed bef_invalid_fit errors", {
  ctx <- fit8_context()
  err <- tryCatch(fit8_assemble(ctx, sha = "not-a-sha"), error = identity)

  expect_s3_class(err, "bef_invalid_fit")
  expect_s3_class(err, "bayesEfron_validate_error")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
})

test_that("Stage 8 validation failures for posterior and metadata stay typed", {
  ctx <- fit8_context()

  negative_sd <- function(...) {
    processed <- fit8_ns("postprocess_stan_draws")(...)
    processed$posterior$theta_sd[1, 1] <- -0.01
    processed
  }
  err <- tryCatch(fit8_assemble(ctx, postprocess_fun = negative_sd), error = identity)
  expect_s3_class(err, "bef_invalid_fit")
  expect_s3_class(err, "bayesEfron_validate_error")

  bad_theta_summary <- function(...) {
    processed <- fit8_ns("postprocess_stan_draws")(...)
    processed$theta_summary$sd[1] <- -0.1
    processed
  }
  err <- tryCatch(
    fit8_assemble(ctx, postprocess_fun = bad_theta_summary),
    error = identity
  )
  expect_s3_class(err, "bef_invalid_fit")
  expect_s3_class(err, "bayesEfron_validate_error")

  bad_metadata_attr <- function(...) {
    processed <- fit8_ns("postprocess_stan_draws")(...)
    processed$diagnostics <- NULL
    processed
  }
  err <- tryCatch(
    fit8_assemble(ctx, postprocess_fun = bad_metadata_attr),
    error = identity
  )
  expect_s3_class(err, "bef_invalid_fit")
  expect_s3_class(err, "bayesEfron_validate_error")

  bad_theta_rep_draws <- function(...) {
    processed <- fit8_ns("postprocess_stan_draws")(...)
    processed$theta_rep_draws[1, 1] <- processed$theta_rep_draws[1, 1] + 1
    processed
  }
  err <- tryCatch(
    fit8_assemble(ctx, postprocess_fun = bad_theta_rep_draws),
    error = identity
  )
  expect_s3_class(err, "bef_invalid_fit")
  expect_s3_class(err, "bayesEfron_validate_error")
})

test_that("Stage 8 Stan SHA helper hashes the package Stan source", {
  expected <- digest::digest(
    fit8_ns(".bef_stan_file")("RE"),
    algo = "sha256",
    file = TRUE
  )

  expect_equal(fit8_ns(".bef_fit_stan_sha256")("RE"), expected)
})
