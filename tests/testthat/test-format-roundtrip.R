format_rt_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

format_rt_summary <- function(mean = 0) {
  list(
    mean = mean,
    sd = 0.1,
    q5 = mean - 0.1,
    q50 = mean,
    q95 = mean + 0.1
  )
}

format_rt_metadata <- function(K = 5L, S = 4L) {
  metadata <- list(
    model_family = "RE",
    grid_method = "paper_realdata",
    seed = 123L,
    cmdstan_version = "2.38.0",
    stan_file_sha256 = paste(rep("f", 64L), collapse = ""),
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
    mean_g_summary = format_rt_summary(0),
    var_g_summary = format_rt_summary(1),
    theta_summary = data.frame(
      site = seq_len(K),
      mean = rep(0, K),
      sd = rep(0.1, K),
      hpdi_lower = rep(-0.2, K),
      hpdi_upper = rep(0.2, K),
      map = rep(0, K)
    ),
    theta_rep_draws = matrix(0, nrow = S, ncol = K),
    effective_params_summary = format_rt_summary(3),
    log_marginal_likelihood_summary = format_rt_summary(-10)
  )
  attr(metadata, "sd_g_summary") <- format_rt_summary(1)
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

format_rt_fit <- function(K = 5L, S = 4L) {
  draws <- array(
    seq_len(S * 1L * 2L),
    dim = c(S, 1L, 2L),
    dimnames = list(NULL, NULL, c("mean_g", "var_g"))
  )
  format_rt_ns("validate_bef_fit_re")(
    format_rt_ns("new_bef_fit_re")(
      draws = draws,
      metadata = format_rt_metadata(K = K, S = S),
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

expect_format_print_roundtrip <- function(x, use_cli) {
  expect_identical(
    capture.output(print(x, use_cli = use_cli)),
    format(x, use_cli = use_cli)
  )
}

test_that("format and print round-trip for all format-print pairs", {
  fit <- format_rt_fit()
  parent_fit <- fit
  class(parent_fit) <- "bef_fit"
  objects <- list(
    parent_fit = parent_fit,
    summary_bef_fit = summary(parent_fit),
    bef_data = as_bef_data(list(
      theta_hat = stats::setNames(seq(-0.4, 0.4, length.out = 5L), LETTERS[1:5]),
      sigma = rep(0.2, 5L)
    )),
    bef_diagnostic = diagnose(fit)
  )

  lapply(objects, expect_format_print_roundtrip, use_cli = FALSE)
  if (requireNamespace("cli", quietly = TRUE)) {
    lapply(objects, expect_format_print_roundtrip, use_cli = TRUE)
    lapply(
      objects,
      function(x) {
        expect_identical(
          as.character(cli::ansi_strip(format(x, use_cli = TRUE))),
          format(x, use_cli = FALSE)
        )
      }
    )
  }

  withr::local_envvar(c(BAYESEFRON_NO_CLI = "1"))
  lapply(objects, expect_format_print_roundtrip, use_cli = NULL)
  lapply(
    objects,
    function(x) expect_identical(format(x, use_cli = NULL), format(x, use_cli = FALSE))
  )
})
