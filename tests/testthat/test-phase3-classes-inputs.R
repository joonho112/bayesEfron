bef_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

phase3_summary <- function(mean = 0) {
  list(
    mean = mean,
    sd = 0.1,
    q5 = mean - 0.1,
    q50 = mean,
    q95 = mean + 0.1
  )
}

phase3_data_list <- function(K = 5L, L = 51L, M = 3L) {
  list(
    K = as.integer(K),
    theta_hat = seq(-0.4, 0.4, length.out = K),
    sigma = rep(0.2, K),
    L = as.integer(L),
    grid = seq(-1, 1, length.out = L),
    M = as.integer(M),
    B = matrix(seq_len(L * M) / 100, nrow = L, ncol = M)
  )
}

phase3_metadata <- function(K = 5L, L = 51L, M = 3L, S = 4L) {
  metadata <- list(
    model_family = "RE",
    grid_method = "paper_realdata",
    seed = 123L,
    cmdstan_version = "2.34.0",
    stan_file_sha256 = paste(rep("a", 64L), collapse = ""),
    data_list = phase3_data_list(K = K, L = L, M = M),
    runtime_seconds = 1.25,
    mean_g_summary = phase3_summary(0),
    var_g_summary = phase3_summary(1),
    theta_summary = data.frame(
      site = seq_len(K),
      mean = rep(0, K),
      sd = rep(0.1, K),
      hpdi_lower = rep(0, K),
      hpdi_upper = rep(0, K),
      map = rep(0, K)
    ),
    theta_rep_draws = matrix(0, nrow = S, ncol = K),
    effective_params_summary = phase3_summary(3),
    log_marginal_likelihood_summary = phase3_summary(-10)
  )
  attr(metadata, "sd_g_summary") <- phase3_summary(1)
  attr(metadata, "diagnostics") <- list(
    rhat = 1,
    ess_bulk = 100,
    ess_tail = 100,
    divergences = 0,
    max_treedepth = 0
  )
  attr(metadata, "diagnostic_skipped") <- character()
  attr(metadata, "sampler_diagnostics_failed") <- character()
  metadata
}

phase3_posterior <- function(K = 5L, S = 4L) {
  list(
    mean_g = rep(0, S),
    var_g = rep(1, S),
    sd_g = rep(1, S),
    theta_map = matrix(0, nrow = S, ncol = K),
    theta_mean = matrix(0, nrow = S, ncol = K),
    theta_sd = matrix(0.1, nrow = S, ncol = K),
    theta_rep = matrix(0, nrow = S, ncol = K),
    effective_params = rep(3, S),
    log_marginal_likelihood = rep(-10, S)
  )
}

phase3_fit_re <- function(K = 5L, L = 51L, M = 3L, S = 4L) {
  bef_ns("new_bef_fit_re")(
    draws = array(seq_len(S * 1L * 2L), dim = c(S, 1L, 2L)),
    metadata = phase3_metadata(K = K, L = L, M = M, S = S),
    posterior = phase3_posterior(K = K, S = S)
  )
}

expect_bef_invalid_args <- function(err) {
  expect_s3_class(err, "bef_invalid_args")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
}

expect_bef_invalid_fit <- function(err) {
  expect_s3_class(err, "bef_invalid_fit")
  expect_s3_class(err, "bayesEfron_validate_error")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
}

test_that("Phase 3 constructors are O(1) structural assemblers", {
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

test_that("Phase 3 validators accept minimal valid synthetic objects", {
  data <- bef_ns("new_bef_data")(
    theta_hat = seq(-0.4, 0.4, length.out = 5L),
    sigma = rep(0.2, 5L),
    names = paste0("study", 1:5),
    source = "test"
  )
  expect_identical(bef_ns("validate_bef_data")(data), data)

  fit <- phase3_fit_re()
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
    effective_params_summary = phase3_summary(3),
    runtime_seconds = 2,
    diagnostic_skipped = "rhat",
    sampler_diagnostics_failed = character()
  )
  expect_identical(bef_ns("validate_bef_diagnostic")(diag), diag)
})

test_that("Phase 3 validator failures carry typed class chains", {
  bad_data <- bef_ns("new_bef_data")(
    theta_hat = 1:4,
    sigma = rep(0.2, 4L)
  )
  err <- tryCatch(bef_ns("validate_bef_data")(bad_data), error = identity)
  expect_bef_invalid_args(err)
  expect_s3_class(err, "bayesEfron_validate_error")

  bad_fit <- phase3_fit_re()
  bad_fit$metadata$model_family <- "HE"
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- phase3_fit_re()
  bad_fit$posterior$theta_sd[1, 1] <- -0.1
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

test_that("Phase 3 fit validator closes the Step 5.2 metadata schema", {
  expected <- bef_ns(".bef_fit_re_metadata_fields")()
  fit <- phase3_fit_re()
  expect_equal(names(fit$metadata), expected)

  bad_fit <- phase3_fit_re()
  bad_fit$metadata$extra <- TRUE
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)
  expect_equal(err$extra_fields, "extra")

  bad_fit <- phase3_fit_re()
  bad_fit$metadata$theta_summary <- NULL
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)
  expect_equal(err$missing_fields, "theta_summary")

  bad_fit <- phase3_fit_re()
  bad_fit$metadata <- rev(bad_fit$metadata)
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- phase3_fit_re()
  class(bad_fit) <- c(class(bad_fit), "extra")
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- phase3_fit_re()
  bad_fit$draws <- array(1, dim = c(4L, 2L))
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)
})

test_that("Phase 3 fit validator rejects malformed summaries and theta outputs", {
  bad_fit <- phase3_fit_re()
  bad_fit$metadata$mean_g_summary$q50 <- NULL
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- phase3_fit_re()
  bad_fit$metadata$mean_g_summary$extra <- 1
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- phase3_fit_re()
  bad_fit$metadata$theta_summary$site <- rev(bad_fit$metadata$theta_summary$site)
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- phase3_fit_re()
  bad_fit$metadata$theta_summary$mean[1] <- 1
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- phase3_fit_re()
  bad_fit$metadata$theta_summary$hpdi_lower[1] <-
    bad_fit$metadata$theta_summary$hpdi_upper[1] + 1
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- phase3_fit_re()
  bad_fit$metadata$theta_rep_draws <- bad_fit$metadata$theta_rep_draws[-1, , drop = FALSE]
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- phase3_fit_re()
  bad_fit$posterior$extra <- 1
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)
  expect_equal(err$extra_fields, "extra")

  bad_fit <- phase3_fit_re()
  attr(bad_fit$metadata, "diagnostic_skipped") <- "rhat"
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- phase3_fit_re()
  diagnostics <- attr(bad_fit$metadata, "diagnostics", exact = TRUE)
  diagnostics$rhat <- NA_real_
  attr(bad_fit$metadata, "diagnostics") <- diagnostics
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)

  bad_fit <- phase3_fit_re()
  attr(bad_fit$metadata, "sampler_diagnostics_failed") <- "not_real"
  err <- tryCatch(bef_ns("validate_bef_fit_re")(bad_fit), error = identity)
  expect_bef_invalid_fit(err)
})

test_that("central typed condition catalog and representative helpers are stable", {
  expect_equal(length(bef_ns(".bef_condition_catalog")()), 19L)

  expect_equal(
    bef_ns(".bef_condition_classes")("bef_invalid_args"),
    c("bef_invalid_args", "bef_pipeline_error", "bef_error")
  )
  expect_equal(
    bef_ns(".bef_condition_classes")("bef_grid_rank_deficient"),
    c("bef_grid_rank_deficient", "bef_grid_error", "bef_error")
  )
  expect_equal(
    bef_ns(".bef_condition_classes")("bef_cache_corruption"),
    c("bef_cache_corruption", "bef_cache_warning", "bef_warning")
  )
  expect_equal(
    bef_ns(".bef_condition_classes")(
      "bef_invalid_fit",
      extra_class = "bayesEfron_validate_error"
    ),
    c(
      "bef_invalid_fit",
      "bayesEfron_validate_error",
      "bef_pipeline_error",
      "bef_error"
    )
  )

  err <- tryCatch(
    bef_ns(".bef_condition_classes")("bef_not_real"),
    error = identity
  )
  expect_s3_class(err, "rlang_error")

  parent <- simpleError("parent")
  err <- tryCatch(
    bef_ns(".bef_abort_invalid_args")(
      "bad args",
      arg = "theta_hat",
      parent = parent
    ),
    error = identity
  )
  expect_bef_invalid_args(err)
  expect_equal(err$arg, "theta_hat")
  expect_identical(err$parent, parent)

  warning_classes <- NULL
  withCallingHandlers(
    bef_ns(".bef_warn_cache_corruption")("cache issue"),
    bef_cache_warning = function(w) {
      warning_classes <<- class(w)
      invokeRestart("muffleWarning")
    }
  )
  expect_true("bef_cache_corruption" %in% warning_classes)
  expect_true("bef_warning" %in% warning_classes)
})

test_that("fit argument validation normalizes valid inputs", {
  out <- bef_ns(".bef_validate_fit_args")(
    theta_hat = c(a = 1L, b = 2L, c = 3L),
    sigma = c(0.2, 0.3, 0.4),
    seed = 7L
  )

  expect_named(
    out,
    c(
      "theta_hat", "sigma", "grid_method", "L", "expansion", "M",
      "theta_true", "bound_expansion", "model_family", "chains",
      "iter_warmup", "iter_sampling", "adapt_delta", "seed",
      "keep_cmdstan_fit"
    )
  )
  expect_type(out$theta_hat, "double")
  expect_equal(names(out$theta_hat), c("a", "b", "c"))
  expect_equal(out$grid_method, "paper_realdata")
  expect_type(out$L, "integer")
  expect_type(out$M, "integer")
  expect_equal(out$model_family, "RE")
  expect_equal(out$seed, 7L)
  expect_false(out$keep_cmdstan_fit)
})

test_that("fit argument validation failures are typed and payloaded", {
  err <- tryCatch(
    bef_ns(".bef_validate_fit_args")(c(1, 2), c(0.1, 0)),
    error = identity
  )
  expect_bef_invalid_args(err)
  expect_equal(err$arg, "sigma")
  expect_equal(err$module, "input-validate")
  expect_equal(err$stage, 2L)
  expect_s3_class(err$parent, "error")

  err <- tryCatch(
    bef_ns(".bef_validate_fit_args")(c(1, 2), c(0.1, 0.2), foo = 1),
    error = identity
  )
  expect_bef_invalid_args(err)
  expect_equal(err$arg, "...")
  expect_equal(err$predicate, "empty dots")

  bad_cases <- list(
    function() bef_ns(".bef_validate_fit_args")(c(1, 2), c(0.1, 0.2), L = NULL),
    function() bef_ns(".bef_validate_fit_args")(c(1, 2), c(0.1, 0.2), L = 50L),
    function() bef_ns(".bef_validate_fit_args")(c(1, 2), c(0.1, 0.2), M = 11L),
    function() bef_ns(".bef_validate_fit_args")(c(1, 2), c(0.1, 0.2), chains = 17L),
    function() bef_ns(".bef_validate_fit_args")(c(1, 2), c(0.1, 0.2), iter_warmup = -1L),
    function() bef_ns(".bef_validate_fit_args")(c(1, 2), c(0.1, 0.2), adapt_delta = 1),
    function() bef_ns(".bef_validate_fit_args")(c(1, 2), c(0.1, 0.2), seed = -1L),
    function() bef_ns(".bef_validate_fit_args")(c(1, 2), c(0.1, 0.2), keep_cmdstan_fit = NA),
    function() bef_ns(".bef_validate_fit_args")(c(1, 2), c(0.1, 0.2), bound_expansion = 0)
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
  expect_equal(err$module, "as-bef-data")
  expect_equal(err$stage, 4L)

  err <- tryCatch(
    as_bef_data(list(theta_hat = 1:4, sigma = rep(0.2, 4L))),
    error = identity
  )
  expect_bef_invalid_args(err)
  expect_s3_class(err, "bayesEfron_validate_error")

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
  expect_s3_class(err, "bayesEfron_validate_error")
})
