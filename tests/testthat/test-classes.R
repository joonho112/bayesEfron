bef_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

synthetic_summary <- function(mean = 0) {
  list(
    mean = mean,
    sd = 0.1,
    q5 = mean - 0.1,
    q50 = mean,
    q95 = mean + 0.1
  )
}

synthetic_data_list <- function(K = 5L, L = 51L, M = 3L) {
  list(
    K = as.integer(K),
    theta_hat = seq(-0.4, 0.4, length.out = K),
    sigma = rep(0.2, K),
    L = as.integer(L),
    grid = seq(-1, 1, length.out = L),
    M = as.integer(M),
    B = matrix(seq_len(L * M) / 100, nrow = L, ncol = M),
    store_grid_quantities = 0L
  )
}

synthetic_metadata <- function(K = 5L, L = 51L, M = 3L, S = 4L) {
  metadata <- list(
    model_family = "RE",
    grid_method = "paper_realdata",
    seed = 123L,
    cmdstan_version = "2.34.0",
    stan_file_sha256 = paste(rep("a", 64L), collapse = ""),
    data_list = synthetic_data_list(K = K, L = L, M = M),
    runtime_seconds = 1.25,
    mean_g_summary = synthetic_summary(0),
    var_g_summary = synthetic_summary(1),
    theta_summary = data.frame(
      site = seq_len(K),
      mean = rep(0, K),
      sd = rep(0.1, K),
      hpdi_lower = rep(0, K),
      hpdi_upper = rep(0, K),
      map = rep(0, K)
    ),
    effective_params_summary = synthetic_summary(3),
    log_marginal_likelihood_summary = synthetic_summary(-10)
  )
  attr(metadata, "sd_g_summary") <- synthetic_summary(1)
  attr(metadata, "diagnostics") <- list(
    rhat = 1,
    ess_bulk = 100,
    ess_tail = 100,
    divergences = 0,
    max_treedepth = 0,
    ebfmi = 0.9
  )
  attr(metadata, "diagnostic_skipped") <- character()
  attr(metadata, "sampler_diagnostics_failed") <- character()
attr(metadata, "sampler_diagnostics_warned") <- character()
  metadata
}

synthetic_posterior <- function(K = 5L, S = 4L) {
  list(
    mean_g = rep(0, S),
    var_g = rep(1, S),
    sd_g = rep(1, S),
    effective_params = rep(3, S),
    log_marginal_likelihood = rep(-10, S)
  )
}

synthetic_fit_re <- function(K = 5L, L = 51L, M = 3L, S = 4L) {
  bef_ns("new_bef_fit_re")(
    draws = bef_fixture_draws(
      theta_map = matrix(0, nrow = S, ncol = K),
      theta_mean = matrix(0, nrow = S, ncol = K),
      theta_sd = matrix(0.1, nrow = S, ncol = K),
      theta_rep = matrix(0, nrow = S, ncol = K)
    ),
    metadata = synthetic_metadata(K = K, L = L, M = M, S = S),
    posterior = synthetic_posterior(K = K, S = S)
  )
}

expect_bef_invalid_args <- function(err) {
  expect_s3_class(err, "bef_invalid_args")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
}

expect_bef_invalid_fit <- function(err) {
  expect_s3_class(err, "bef_invalid_fit")
  expect_s3_class(err, "bef_validate_error")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
}

test_that("constructors assemble objects without validating them", {
  data <- bef_ns("new_bef_data")(
    theta_hat = 1,
    sigma = -1,
    names = NULL,
    source = ""
  )
  expect_s3_class(data, "bef_data")
  expect_named(data, c("theta_hat", "sigma", "names", "source"))
  expect_equal(data$theta_hat, 1)
  expect_equal(data$sigma, -1)

  fit <- bef_ns("new_bef_fit")(draws = "bad", metadata = NULL)
  expect_s3_class(fit, "bef_fit")
  expect_named(fit, c("draws", "metadata", "posterior"))
  expect_false("cmdstan_fit" %in% names(fit))

  fit_re <- bef_ns("new_bef_fit_re")(
    draws = "bad",
    metadata = NULL,
    cmdstan_fit = list(raw = TRUE)
  )
  expect_equal(class(fit_re)[1:2], c("bef_fit_re", "bef_fit"))
  expect_true("cmdstan_fit" %in% names(fit_re))

  diag <- bef_ns("new_bef_diagnostic")(
    rhat = NA_real_,
    ess_bulk = NA_real_,
    ess_tail = NA_real_,
    divergences = 0,
    max_treedepth = 0,
    model_family = "BAD",
    stan_file_sha256 = "not-a-sha"
  )
  expect_s3_class(diag, "bef_diagnostic")
})

test_that("validators accept minimal valid objects", {
  data <- bef_ns("new_bef_data")(
    theta_hat = seq(-0.4, 0.4, length.out = 5L),
    sigma = rep(0.2, 5L),
    names = paste0("study", 1:5),
    source = "test"
  )
  expect_identical(bef_ns("validate_bef_data")(data), data)

  fit <- synthetic_fit_re()
  expect_identical(bef_ns("validate_bef_fit")(fit), fit)
  expect_identical(bef_ns("validate_bef_fit_re")(fit), fit)

  diag <- bef_ns("new_bef_diagnostic")(
    rhat = c(1.01, NA_real_),
    ess_bulk = c(100, NA_real_),
    ess_tail = c(90, NA_real_),
    divergences = c(0, NA_real_),
    max_treedepth = c(5, NA_real_),
    model_family = "RE",
    stan_file_sha256 = paste(rep("b", 64L), collapse = ""),
    effective_params_summary = synthetic_summary(3),
    runtime_seconds = 2,
    diagnostic_skipped = "rhat",
    sampler_diagnostics_failed = character()
  )
  expect_identical(bef_ns("validate_bef_diagnostic")(diag), diag)
})

test_that("validator failures carry the package condition classes", {
  bad_data <- bef_ns("new_bef_data")(
    theta_hat = 1:4,
    sigma = rep(0.2, 4L)
  )
  err <- tryCatch(bef_ns("validate_bef_data")(bad_data), error = identity)
  expect_bef_invalid_args(err)
  expect_s3_class(err, "bef_validate_error")

  bad_fit <- synthetic_fit_re()
  bad_fit$metadata$model_family <- "HE"
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- synthetic_fit_re()
  bad_fit$posterior$sd_g[1] <- NA_real_
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_diag <- bef_ns("new_bef_diagnostic")(
    rhat = Inf,
    ess_bulk = 100,
    ess_tail = 90,
    divergences = 0,
    max_treedepth = 5,
    model_family = "RE",
    stan_file_sha256 = paste(rep("b", 64L), collapse = "")
  )
  err <- tryCatch(bef_ns("validate_bef_diagnostic")(bad_diag), error = identity)
  expect_bef_invalid_fit(err)
})

test_that("the fit validator rejects missing, extra and reordered metadata fields", {
  expected <- bef_ns(".bef_fit_re_metadata_fields")()
  fit <- synthetic_fit_re()
  expect_equal(names(fit$metadata), expected)

  bad_fit <- synthetic_fit_re()
  bad_fit$metadata$extra <- TRUE
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)
  expect_equal(err$extra_fields, "extra")

  bad_fit <- synthetic_fit_re()
  bad_fit$metadata$theta_summary <- NULL
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)
  expect_equal(err$missing_fields, "theta_summary")

  bad_fit <- synthetic_fit_re()
  bad_fit$metadata <- rev(bad_fit$metadata)
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- synthetic_fit_re()
  class(bad_fit) <- c(class(bad_fit), "extra")
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- synthetic_fit_re()
  bad_fit$draws <- array(1, dim = c(4L, 2L))
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)
})

test_that("the fit validator rejects malformed summaries and site-level output", {
  bad_fit <- synthetic_fit_re()
  bad_fit$metadata$mean_g_summary$q50 <- NULL
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- synthetic_fit_re()
  bad_fit$metadata$mean_g_summary$extra <- 1
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- synthetic_fit_re()
  bad_fit$metadata$theta_summary$site <- rev(bad_fit$metadata$theta_summary$site)
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- synthetic_fit_re()
  bad_fit$metadata$theta_summary$mean[1] <- 1
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  # `sd` must be the posterior standard deviation that the draws imply.
  bad_fit <- synthetic_fit_re()
  bad_fit$metadata$theta_summary$sd[1] <- 0.05
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)
  expect_match(conditionMessage(err), "posterior standard deviation")

  # `map` must be a grid point.
  bad_fit <- synthetic_fit_re()
  grid <- bad_fit$metadata$data_list$grid
  bad_fit$metadata$theta_summary$map[1] <- mean(grid[1:2])
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)
  expect_match(conditionMessage(err), "points of the grid")

  bad_fit <- synthetic_fit_re()
  bad_fit$metadata$theta_summary$hpdi_lower[1] <-
    bad_fit$metadata$theta_summary$hpdi_upper[1] + 1
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  # The site-level quantities are stored only in `draws`, so the malformed case
  # is a draws array that is missing one of them.
  bad_fit <- synthetic_fit_re()
  keep <- dimnames(bad_fit$draws)[[3L]] != "theta_rep[1]"
  bad_fit$draws <- bad_fit$draws[, , keep, drop = FALSE]
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- synthetic_fit_re()
  bad_fit$posterior$extra <- 1
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)
  expect_equal(err$extra_fields, "extra")

  bad_fit <- synthetic_fit_re()
  attr(bad_fit$metadata, "diagnostic_skipped") <- "rhat"
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- synthetic_fit_re()
  diagnostics <- attr(bad_fit$metadata, "diagnostics", exact = TRUE)
  diagnostics$rhat <- NA_real_
  attr(bad_fit$metadata, "diagnostics") <- diagnostics
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- synthetic_fit_re()
  attr(bad_fit$metadata, "sampler_diagnostics_failed") <- "not_real"
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)
})

test_that("fit argument validation normalizes valid inputs", {
  out <- bef_ns(".bef_validate_fit_args")(
    theta_hat = c(a = 1L, b = 2L, c = 3L, d = 4L, e = 5L),
    sigma = c(0.2, 0.3, 0.4, 0.3, 0.2),
    seed = 7L
  )

  expect_named(
    out,
    c(
      "theta_hat", "sigma", "grid_method", "L", "expansion", "M",
      "theta_true", "bound_expansion", "model_family", "chains",
      "parallel_chains", "iter_warmup", "iter_sampling", "adapt_delta",
      "max_treedepth", "seed",
      "keep_cmdstan_fit", "store_grid_quantities"
    )
  )
  expect_type(out$theta_hat, "double")
  expect_equal(names(out$theta_hat), c("a", "b", "c", "d", "e"))
  expect_equal(out$grid_method, "paper_realdata")
  expect_type(out$L, "integer")
  expect_type(out$M, "integer")
  expect_equal(out$model_family, "RE")
  expect_equal(out$seed, 7L)
  expect_false(out$keep_cmdstan_fit)
})

test_that("fit argument validation failures are typed and payloaded", {
  err <- tryCatch(
    bef_ns(".bef_validate_fit_args")(1:5 / 10, c(0.1, 0.2, 0.1, 0.2, 0)),
    error = identity
  )
  expect_bef_invalid_args(err)
  expect_equal(err$arg, "sigma")
  expect_s3_class(err$parent, "error")

  err <- tryCatch(
    bef_ns(".bef_validate_fit_args")(1:5 / 10, rep(0.2, 5), foo = 1),
    error = identity
  )
  expect_bef_invalid_args(err)
  expect_equal(err$arg, "...")

  bad_cases <- list(
    function() bef_ns(".bef_validate_fit_args")(1:5 / 10, rep(0.2, 5), L = NULL),
    function() bef_ns(".bef_validate_fit_args")(1:5 / 10, rep(0.2, 5), L = 50L),
    function() bef_ns(".bef_validate_fit_args")(1:5 / 10, rep(0.2, 5), M = 11L),
    function() bef_ns(".bef_validate_fit_args")(1:5 / 10, rep(0.2, 5), chains = 17L),
    function() bef_ns(".bef_validate_fit_args")(1:5 / 10, rep(0.2, 5), iter_warmup = -1L),
    function() bef_ns(".bef_validate_fit_args")(1:5 / 10, rep(0.2, 5), adapt_delta = 1),
    function() bef_ns(".bef_validate_fit_args")(1:5 / 10, rep(0.2, 5), seed = -1L),
    function() bef_ns(".bef_validate_fit_args")(1:5 / 10, rep(0.2, 5), keep_cmdstan_fit = NA),
    function() bef_ns(".bef_validate_fit_args")(1:5 / 10, rep(0.2, 5), bound_expansion = 0)
  )
  errs <- lapply(
    bad_cases,
    function(case) tryCatch(case(), error = identity)
  )
  expect_true(all(vapply(errs, inherits, logical(1), "bef_invalid_args")))
})

test_that("as_bef_data converts lists and preserves labels", {
  x <- as_bef_data(list(
    theta_hat = c(a = 1, b = 2, c = 3, d = 4, e = 5),
    sigma = rep(0.2, 5L)
  ))

  expect_s3_class(x, "bef_data")
  expect_equal(x$source, "list")
  expect_equal(x$names, letters[1:5])
  expect_equal(x$sigma, rep(0.2, 5L))

  explicit <- as_bef_data(list(
    theta_hat = c(a = 1, b = 2, c = 3, d = 4, e = 5),
    sigma = rep(0.2, 5L),
    names = paste0("site", 1:5),
    ignored = "extra"
  ))
  expect_equal(explicit$names, paste0("site", 1:5))

  unnamed <- as_bef_data(list(
    theta_hat = c(a = 1, b = 2, 3, d = 4, e = 5),
    sigma = rep(0.2, 5L)
  ))
  expect_null(unnamed$names)
})

test_that("as_bef_data converts escalc-class data without metafor", {
  es <- data.frame(
    yi = seq(0.1, 0.5, length.out = 5L),
    vi = rep(0.04, 5L),
    slab = paste0("ignored", 1:5),
    row.names = paste0("study", 1:5)
  )
  class(es) <- c("escalc", "data.frame")

  x <- as_bef_data(es)

  expect_s3_class(x, "bef_data")
  expect_equal(x$source, "metafor::escalc")
  expect_equal(x$theta_hat, es$yi)
  expect_equal(x$sigma, rep(0.2, 5L))
  expect_equal(x$names, paste0("study", 1:5))

  row.names(es) <- as.character(seq_len(5L))
  x <- as_bef_data(es)
  expect_null(x$names)
})

test_that("as_bef_data conversion failures are typed", {
  err <- tryCatch(as_bef_data(1:5), error = identity)
  expect_bef_invalid_args(err)

  err <- tryCatch(
    as_bef_data(list(theta_hat = 1:4, sigma = rep(0.2, 4L))),
    error = identity
  )
  expect_bef_invalid_args(err)
  expect_s3_class(err, "bef_validate_error")

  err <- tryCatch(
    as_bef_data(list(theta_hat = 1:5, sigma = rep(0.2, 5L)), extra = 1),
    error = identity
  )
  expect_bef_invalid_args(err)
  expect_equal(err$arg, "...")

  es <- data.frame(yi = 1:5, vi = c(0.04, 0.04, -0.1, 0.04, 0.04))
  class(es) <- c("escalc", "data.frame")
  err <- tryCatch(as_bef_data(es), error = identity)
  expect_bef_invalid_args(err)
  expect_equal(err$arg, "x$vi")

  err <- tryCatch(as_bef_data(data.frame(yi = 1:5, vi = rep(0.04, 5L))), error = identity)
  expect_bef_invalid_args(err)
})

test_that("as_bef_data validates existing bef_data input", {
  data <- bef_ns("new_bef_data")(
    theta_hat = seq(-0.4, 0.4, length.out = 5L),
    sigma = rep(0.2, 5L),
    names = NULL,
    source = "test"
  )
  expect_identical(as_bef_data(data), data)

  bad_data <- data
  bad_data$sigma[1] <- 0
  err <- tryCatch(as_bef_data(bad_data), error = identity)
  expect_bef_invalid_args(err)
  expect_s3_class(err, "bef_validate_error")
})
