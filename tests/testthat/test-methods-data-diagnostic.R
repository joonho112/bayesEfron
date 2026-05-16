diag_methods_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

diag_methods_summary <- function(mean = 0) {
  list(
    mean = mean,
    sd = 0.1,
    q5 = mean - 0.1,
    q50 = mean,
    q95 = mean + 0.1
  )
}

diag_methods_metadata <- function(K = 5L, S = 4L) {
  metadata <- list(
    model_family = "RE",
    grid_method = "paper_realdata",
    seed = 123L,
    cmdstan_version = "2.38.0",
    stan_file_sha256 = paste(rep("d", 64L), collapse = ""),
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
    mean_g_summary = diag_methods_summary(0),
    var_g_summary = diag_methods_summary(1),
    theta_summary = data.frame(
      site = seq_len(K),
      mean = rep(0, K),
      sd = rep(0.1, K),
      hpdi_lower = rep(-0.2, K),
      hpdi_upper = rep(0.2, K),
      map = rep(0, K)
    ),
    theta_rep_draws = matrix(0, nrow = S, ncol = K),
    effective_params_summary = diag_methods_summary(3),
    log_marginal_likelihood_summary = diag_methods_summary(-10)
  )
  attr(metadata, "sd_g_summary") <- diag_methods_summary(1)
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

diag_methods_fit <- function(K = 5L, S = 4L) {
  draws <- array(
    seq_len(S * 1L * 2L),
    dim = c(S, 1L, 2L),
    dimnames = list(NULL, NULL, c("mean_g", "var_g"))
  )
  diag_methods_ns("validate_bef_fit_re")(
    diag_methods_ns("new_bef_fit_re")(
      draws = draws,
      metadata = diag_methods_metadata(K = K, S = S),
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

test_that("bef_data summary, format, and print expose compact input diagnostics", {
  data <- as_bef_data(list(
    theta_hat = stats::setNames(seq(-0.4, 0.4, length.out = 5L), LETTERS[1:5]),
    sigma = rep(0.2, 5L)
  ))
  data_summary <- summary(data)

  expect_named(data_summary, c("K", "theta_hat", "sigma", "source", "names"))
  expect_equal(data_summary$K, 5L)
  expect_equal(data_summary$theta_hat$min, -0.4)
  expect_equal(data_summary$sigma$median, 0.2)
  expect_equal(data_summary$source, "list")
  expect_equal(data_summary$names, LETTERS[1:5])

  lines <- format(data)
  expect_type(lines, "character")
  expect_true(any(grepl("<bef_data>", lines, fixed = TRUE)))
  expect_true(any(grepl("Sites: 5", lines, fixed = TRUE)))
  expect_true(any(grepl("theta_hat:", lines, fixed = TRUE)))
  expect_true(any(grepl("sigma:", lines, fixed = TRUE)))

  data_print <- capture.output(out <- print(data))
  expect_identical(data_print, lines)
  expect_identical(out, data)
})

test_that("diagnose returns a validated bef_diagnostic from fit metadata attributes", {
  fit <- diag_methods_fit()
  diagnostic <- diagnose(fit)

  expect_s3_class(diagnostic, "bef_diagnostic")
  expect_named(diagnostic, diag_methods_ns(".bef_diagnostic_fields")())
  expect_equal(diagnostic$rhat, 1.01)
  expect_equal(diagnostic$ess_bulk, 500)
  expect_equal(diagnostic$ess_tail, 450)
  expect_equal(diagnostic$divergences, 0)
  expect_equal(diagnostic$model_family, "RE")
  expect_equal(diagnostic$effective_params_summary$mean, 3)
  expect_equal(diagnostic$runtime_seconds, 2.5)
  expect_identical(diagnostic$diagnostic_skipped, character())
  expect_identical(diagnostic$sampler_diagnostics_failed, character())
})

test_that("diagnose preserves skipped and failed diagnostic flags", {
  fit <- diag_methods_fit()
  diagnostics <- attr(fit$metadata, "diagnostics", exact = TRUE)
  diagnostics$rhat <- NA_real_
  attr(fit$metadata, "diagnostics") <- diagnostics
  attr(fit$metadata, "diagnostic_skipped") <- "rhat"
  attr(fit$metadata, "sampler_diagnostics_failed") <- "divergences"

  diagnostic <- diagnose(fit)
  diagnostic_summary <- summary(diagnostic)

  expect_true(is.na(diagnostic$rhat))
  expect_identical(diagnostic$diagnostic_skipped, "rhat")
  expect_identical(diagnostic$sampler_diagnostics_failed, "divergences")
  expect_true(is.na(diagnostic_summary$rhat$value))
  expect_true(is.na(diagnostic_summary$rhat$index))

  lines <- format(diagnostic)
  expect_true(any(grepl("Skipped diagnostics: rhat", lines, fixed = TRUE)))
  expect_true(any(grepl(
    "Diagnostics over warning thresholds: divergences",
    lines,
    fixed = TRUE
  )))
})

test_that("bef_diagnostic summary, format, and print expose sampler health", {
  diagnostic <- diagnose(diag_methods_fit())
  diagnostic_summary <- summary(diagnostic)

  expect_named(
    diagnostic_summary,
    c(
      "rhat", "ess_bulk", "ess_tail", "divergences", "max_treedepth",
      "effective_params", "model_family", "stan_file_sha256",
      "runtime_seconds", "diagnostic_skipped", "sampler_diagnostics_failed"
    )
  )
  expect_equal(diagnostic_summary$rhat$value, 1.01)
  expect_equal(diagnostic_summary$rhat$index, 1L)
  expect_equal(diagnostic_summary$ess_bulk$value, 500)
  expect_equal(diagnostic_summary$ess_tail$value, 450)

  lines <- format(diagnostic)
  expect_type(lines, "character")
  expect_true(any(grepl("<bef_diagnostic>", lines, fixed = TRUE)))
  expect_true(any(grepl("Rhat max:", lines, fixed = TRUE)))
  expect_true(any(grepl("ESS bulk min:", lines, fixed = TRUE)))
  expect_true(any(grepl("Divergences:", lines, fixed = TRUE)))

  diagnostic_print <- capture.output(out <- print(diagnostic))
  expect_identical(diagnostic_print, lines)
  expect_identical(out, diagnostic)
})

test_that("diagnose argument failures are typed", {
  err <- tryCatch(diagnose(list()), error = identity)
  expect_s3_class(err, "bef_invalid_fit")

  err <- tryCatch(diagnose(diag_methods_fit(), extra = 1), error = identity)
  expect_s3_class(err, "bef_invalid_args")
  expect_equal(err$arg, "...")
})

test_that("bef_data and bef_diagnostic format expose optional cli path with escape hatch", {
  data <- as_bef_data(list(
    theta_hat = stats::setNames(seq(-0.4, 0.4, length.out = 5L), LETTERS[1:5]),
    sigma = rep(0.2, 5L)
  ))
  diagnostic <- diagnose(diag_methods_fit())

  expect_error(format(data, use_cli = NA), class = "bef_invalid_args")
  expect_error(format(diagnostic, use_cli = c(TRUE, FALSE)), class = "bef_invalid_args")

  if (requireNamespace("cli", quietly = TRUE)) {
    expect_length(
      format(data, use_cli = TRUE),
      length(format(data, use_cli = FALSE))
    )
    expect_length(
      format(diagnostic, use_cli = TRUE),
      length(format(diagnostic, use_cli = FALSE))
    )
  }

  withr::local_envvar(c(BAYESEFRON_NO_CLI = "1"))
  expect_identical(format(data, use_cli = TRUE), format(data, use_cli = FALSE))
  expect_identical(
    format(diagnostic, use_cli = TRUE),
    format(diagnostic, use_cli = FALSE)
  )
})
