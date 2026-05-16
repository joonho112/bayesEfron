methods_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

methods_summary <- function(mean = 0) {
  list(
    mean = mean,
    sd = 0.1,
    q5 = mean - 0.1,
    q50 = mean,
    q95 = mean + 0.1
  )
}

methods_metadata <- function(K = 5L, S = 4L) {
  metadata <- list(
    model_family = "RE",
    grid_method = "paper_realdata",
    seed = 123L,
    cmdstan_version = "2.38.0",
    stan_file_sha256 = paste(rep("c", 64L), collapse = ""),
    data_list = list(
      K = as.integer(K),
      theta_hat = seq(-0.4, 0.4, length.out = K),
      sigma = rep(0.2, K),
      L = 51L,
      grid = seq(-1, 1, length.out = 51L),
      M = 3L,
      B = matrix(seq_len(51L * 3L) / 100, nrow = 51L, ncol = 3L)
    ),
    runtime_seconds = 2.5,
    mean_g_summary = methods_summary(0),
    var_g_summary = methods_summary(1),
    theta_summary = data.frame(
      site = seq_len(K),
      mean = rep(0, K),
      sd = rep(0.1, K),
      hpdi_lower = rep(0, K),
      hpdi_upper = rep(0, K),
      map = rep(0, K)
    ),
    theta_rep_draws = matrix(0, nrow = S, ncol = K),
    effective_params_summary = methods_summary(3),
    log_marginal_likelihood_summary = methods_summary(-10)
  )
  attr(metadata, "sd_g_summary") <- methods_summary(1)
  attr(metadata, "diagnostics") <- list(
    rhat = 1.01,
    ess_bulk = 500,
    ess_tail = 450,
    divergences = 0,
    max_treedepth = 0
  )
  attr(metadata, "diagnostic_skipped") <- character()
  attr(metadata, "sampler_diagnostics_failed") <- character()
  metadata
}

methods_fit <- function(K = 5L, S = 4L) {
  draws <- array(
    seq_len(S * 1L * 2L),
    dim = c(S, 1L, 2L),
    dimnames = list(NULL, NULL, c("mean_g", "var_g"))
  )
  methods_ns("validate_bef_fit_re")(
    methods_ns("new_bef_fit_re")(
      draws = draws,
      metadata = methods_metadata(K = K, S = S),
      posterior = list(
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
    )
  )
}

test_that("parent bef_fit methods expose portable universal summaries", {
  fit <- methods_fit()
  class(fit) <- "bef_fit"
  fit_summary <- summary(fit)

  expect_s3_class(fit_summary, "summary.bef_fit")
  expect_named(fit_summary, c("prior_summary", "diagnostics"))
  expect_equal(attr(fit_summary, "level", exact = TRUE), 0.9)
  expect_named(fit_summary$prior_summary, c("mean", "var", "sd"))
  expect_equal(fit_summary$prior_summary$mean, 0)
  expect_equal(fit_summary$prior_summary$var, 1)
  expect_equal(fit_summary$prior_summary$sd, 1)
  expect_false("theta_summary" %in% names(fit_summary))

  expect_named(
    fit_summary$diagnostics,
    c(
      "rhat", "ess_bulk", "ess_tail", "divergences", "max_treedepth",
      "effective_params", "log_marginal_likelihood", "model_family",
      "stan_file_sha256", "runtime_seconds", "diagnostic_skipped",
      "sampler_diagnostics_failed"
    )
  )
  expect_equal(fit_summary$diagnostics$effective_params$mean, 3)
  expect_equal(fit_summary$diagnostics$log_marginal_likelihood$mean, -10)
})

test_that("parent bef_fit print and format return stable character surfaces", {
  fit <- methods_fit()
  class(fit) <- "bef_fit"
  lines <- format(fit)

  expect_type(lines, "character")
  expect_true(any(grepl("<bayesEfron fit>", lines, fixed = TRUE)))
  expect_true(any(grepl("Sites: 5", lines, fixed = TRUE)))
  expect_true(any(grepl("Diagnostics:", lines, fixed = TRUE)))
  fit_print <- capture.output(out <- print(fit))
  expect_identical(fit_print, lines)
  expect_identical(out, fit)

  fit_summary <- summary(fit)
  summary_lines <- format(fit_summary)
  expect_type(summary_lines, "character")
  expect_true(any(grepl("<summary.bef_fit>", summary_lines, fixed = TRUE)))
  expect_true(any(grepl("Prior g:", summary_lines, fixed = TRUE)))
  expect_true(any(grepl("Log marginal lik.", summary_lines, fixed = TRUE)))
  summary_print <- capture.output(out <- print(fit_summary))
  expect_identical(summary_print, summary_lines)
  expect_identical(out, fit_summary)
})

test_that("parent bef_fit nobs, logLik, and as_draws use universal fields", {
  fit <- methods_fit()
  class(fit) <- "bef_fit"
  ll <- logLik(fit)

  expect_identical(nobs(fit), 5L)
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), -10)
  expect_equal(attr(ll, "df", exact = TRUE), 3)
  expect_equal(attr(ll, "nobs", exact = TRUE), 5L)
  draws <- posterior::as_draws(fit)
  expect_s3_class(draws, "draws_array")
  expect_equal(dim(draws), dim(fit$draws))
  expect_equal(dimnames(draws)[[3L]], dimnames(fit$draws)[[3L]])
})

test_that("summary.bef_fit validates level", {
  fit <- methods_fit()
  class(fit) <- "bef_fit"

  err <- tryCatch(summary(fit, level = 1), error = identity)
  expect_s3_class(err, "bef_invalid_args")
  expect_equal(err$arg, "level")

  err <- tryCatch(summary(fit, level = NA_real_), error = identity)
  expect_s3_class(err, "bef_invalid_args")
  expect_equal(err$arg, "level")
})

test_that("parent bef_fit format surfaces diagnostic skipped and failed flags", {
  fit <- methods_fit()
  class(fit) <- "bef_fit"
  diagnostics <- attr(fit$metadata, "diagnostics", exact = TRUE)
  diagnostics$rhat <- NA_real_
  attr(fit$metadata, "diagnostics") <- diagnostics
  attr(fit$metadata, "diagnostic_skipped") <- "rhat"
  attr(fit$metadata, "sampler_diagnostics_failed") <- "divergences"

  lines <- format(fit)
  expect_true(any(grepl("Skipped diagnostics: rhat", lines, fixed = TRUE)))
  expect_true(any(grepl("Diagnostics over warning thresholds: divergences", lines, fixed = TRUE)))

  summary_lines <- format(summary(fit))
  expect_true(any(grepl("Skipped:", summary_lines, fixed = TRUE)))
  expect_true(any(grepl("Warning flags:", summary_lines, fixed = TRUE)))
})

test_that("parent bef_fit format honors base formatting escape hatch", {
  fit <- methods_fit()
  class(fit) <- "bef_fit"
  withr::local_envvar(c(BAYESEFRON_NO_CLI = "1"))

  expect_identical(format(fit, use_cli = NULL), format(fit, use_cli = FALSE))
  expect_identical(
    format(summary(fit), use_cli = NULL),
    format(summary(fit), use_cli = FALSE)
  )
})

test_that("parent bef_fit format exposes optional cli path with stable shape", {
  fit <- methods_fit()
  class(fit) <- "bef_fit"
  fit_summary <- summary(fit)

  expect_error(format(fit, use_cli = NA), class = "bef_invalid_args")
  expect_error(format(fit_summary, use_cli = c(TRUE, FALSE)), class = "bef_invalid_args")

  if (requireNamespace("cli", quietly = TRUE)) {
    expect_length(
      format(fit, use_cli = TRUE),
      length(format(fit, use_cli = FALSE))
    )
    expect_length(
      format(fit_summary, use_cli = TRUE),
      length(format(fit_summary, use_cli = FALSE))
    )
  }

  withr::local_envvar(c(BAYESEFRON_NO_CLI = "1"))
  expect_identical(format(fit, use_cli = TRUE), format(fit, use_cli = FALSE))
  expect_identical(
    format(fit_summary, use_cli = TRUE),
    format(fit_summary, use_cli = FALSE)
  )
})
