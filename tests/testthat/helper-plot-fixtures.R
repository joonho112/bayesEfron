plot_test_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

plot_test_summary <- function(mean = 0) {
  list(
    mean = mean,
    sd = 0.1,
    q5 = mean - 0.1,
    q50 = mean,
    q95 = mean + 0.1
  )
}

plot_test_fit <- function(K = 5L, S = 6L, L = 11L) {
  theta_rep <- matrix(
    seq(-0.5, 0.5, length.out = S * K),
    nrow = S,
    ncol = K
  )
  theta_mean <- matrix(rep(seq(-0.2, 0.2, length.out = K), each = S), nrow = S)
  theta_sd <- matrix(rep(seq(0.1, 0.2, length.out = K), each = S), nrow = S)
  theta_map <- matrix(rep(seq(-0.25, 0.15, length.out = K), each = S), nrow = S)
  grid <- seq(-1, 1, length.out = L)
  # theta_mean and theta_sd are the same in every draw, so the posterior
  # standard deviation is the standard deviation given g. The mode in the
  # table has to be a grid point.
  theta_summary <- data.frame(
    site = seq_len(K),
    mean = colMeans(theta_mean),
    sd = colMeans(theta_sd),
    hpdi_lower = apply(theta_rep, 2L, stats::quantile, probs = 0.05, names = FALSE),
    hpdi_upper = apply(theta_rep, 2L, stats::quantile, probs = 0.95, names = FALSE),
    map = grid[vapply(
      colMeans(theta_map),
      function(value) which.min(abs(value - grid)),
      integer(1L)
    )]
  )
  metadata <- list(
    model_family = "RE",
    grid_method = "paper_realdata",
    seed = 123L,
    cmdstan_version = "2.38.0",
    stan_file_sha256 = paste(rep("9", 64L), collapse = ""),
    data_list = list(
      K = as.integer(K),
      theta_hat = seq(-0.4, 0.4, length.out = K),
      sigma = seq(0.1, 0.3, length.out = K),
      L = as.integer(L),
      grid = grid,
      M = 3L,
      B = matrix(seq_len(L * 3L) / 100, nrow = L, ncol = 3L),
      store_grid_quantities = 0L
    ),
    runtime_seconds = 2.5,
    mean_g_summary = plot_test_summary(0),
    var_g_summary = plot_test_summary(1),
    theta_summary = theta_summary,
    effective_params_summary = plot_test_summary(3),
    log_marginal_likelihood_summary = plot_test_summary(-10)
  )
  attr(metadata, "sd_g_summary") <- plot_test_summary(1)
  attr(metadata, "diagnostics") <- list(
    rhat = 1.01,
    ess_bulk = 500,
    ess_tail = 450,
    divergences = 0,
    max_treedepth = 0,
    ebfmi = 0.9
  )
  attr(metadata, "diagnostic_skipped") <- character()
  attr(metadata, "sampler_diagnostics_failed") <- character()
attr(metadata, "sampler_diagnostics_warned") <- character()

  plot_test_ns("validate_bef_fit_re")(
    plot_test_ns("new_bef_fit_re")(
      draws = bef_fixture_draws(
        theta_map = theta_map, theta_mean = theta_mean,
        theta_sd = theta_sd, theta_rep = theta_rep,
        scalars = list(
          mean_g = seq(-0.1, 0.1, length.out = S),
          var_g = seq(0.8, 1.2, length.out = S)
        ),
        grid_points = L
      ),
      metadata = metadata,
      posterior = list(
        mean_g = seq(-0.1, 0.1, length.out = S),
        var_g = seq(0.8, 1.2, length.out = S),
        sd_g = seq(0.9, 1.1, length.out = S),
        effective_params = rep(3, S),
        log_marginal_likelihood = rep(-10, S)
      )
    )
  )
}

plot_test_with_pdf_device <- function(code) {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}
