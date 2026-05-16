re_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

re_summary <- function(mean = 0) {
  list(
    mean = mean,
    sd = 0.1,
    q5 = mean - 0.1,
    q50 = mean,
    q95 = mean + 0.1
  )
}

re_fit <- function(K = 5L, S = 6L) {
  theta_rep <- matrix(
    seq(-0.5, 0.5, length.out = S * K),
    nrow = S,
    ncol = K
  )
  theta_mean <- matrix(rep(seq(-0.2, 0.2, length.out = K), each = S), nrow = S)
  theta_sd <- matrix(rep(seq(0.1, 0.2, length.out = K), each = S), nrow = S)
  theta_map <- matrix(rep(seq(-0.25, 0.15, length.out = K), each = S), nrow = S)
  theta_summary <- data.frame(
    site = seq_len(K),
    mean = colMeans(theta_mean),
    sd = colMeans(theta_sd),
    hpdi_lower = vapply(
      seq_len(K),
      function(i) posterior::quantile2(theta_rep[, i], probs = 0.05, names = FALSE),
      numeric(1L)
    ),
    hpdi_upper = vapply(
      seq_len(K),
      function(i) posterior::quantile2(theta_rep[, i], probs = 0.95, names = FALSE),
      numeric(1L)
    ),
    map = colMeans(theta_map)
  )
  metadata <- list(
    model_family = "RE",
    grid_method = "paper_realdata",
    seed = 123L,
    cmdstan_version = "2.38.0",
    stan_file_sha256 = paste(rep("d", 64L), collapse = ""),
    data_list = list(
      K = as.integer(K),
      theta_hat = seq(-0.4, 0.4, length.out = K),
      sigma = seq(0.1, 0.3, length.out = K),
      L = 51L,
      grid = seq(-1, 1, length.out = 51L),
      M = 3L,
      B = matrix(seq_len(51L * 3L) / 100, nrow = 51L, ncol = 3L)
    ),
    runtime_seconds = 2.5,
    mean_g_summary = re_summary(0),
    var_g_summary = re_summary(1),
    theta_summary = theta_summary,
    theta_rep_draws = theta_rep,
    effective_params_summary = re_summary(3),
    log_marginal_likelihood_summary = re_summary(-10)
  )
  attr(metadata, "sd_g_summary") <- re_summary(1)
  attr(metadata, "diagnostics") <- list(
    rhat = 1.01,
    ess_bulk = 500,
    ess_tail = 450,
    divergences = 0,
    max_treedepth = 0
  )
  attr(metadata, "diagnostic_skipped") <- character()
  attr(metadata, "sampler_diagnostics_failed") <- character()
  draw_variables <- c("mean_g", "var_g", paste0("g[", seq_len(metadata$data_list$L), "]"))

  re_ns("validate_bef_fit_re")(
    re_ns("new_bef_fit_re")(
      draws = array(
        seq_len(S * 1L * length(draw_variables)) / 100,
        dim = c(S, 1L, length(draw_variables)),
        dimnames = list(NULL, NULL, draw_variables)
      ),
      metadata = metadata,
      posterior = list(
        mean_g = seq(-0.1, 0.1, length.out = S),
        var_g = seq(0.8, 1.2, length.out = S),
        sd_g = seq(0.9, 1.1, length.out = S),
        theta_map = theta_map,
        theta_mean = theta_mean,
        theta_sd = theta_sd,
        theta_rep = theta_rep,
        effective_params = rep(3, S),
        log_marginal_likelihood = rep(-10, S)
      )
    )
  )
}

test_that("RE summary augments the parent shell with theta_summary", {
  fit <- re_fit()
  fit_summary <- summary(fit)

  expect_s3_class(fit_summary, "summary.bef_fit_re")
  expect_s3_class(fit_summary, "summary.bef_fit")
  expect_named(fit_summary, c("prior_summary", "diagnostics", "theta_summary"))
  expect_identical(fit_summary$theta_summary, fit$metadata$theta_summary)
  expect_true(any(grepl("Theta summary:", format(fit_summary), fixed = TRUE)))

  fit_summary_80 <- summary(fit, level = 0.8)
  expect_equal(attr(fit_summary_80, "level", exact = TRUE), 0.8)
  expect_false(identical(
    fit_summary_80$theta_summary$hpdi_lower,
    fit$metadata$theta_summary$hpdi_lower
  ))
})

test_that("RE coef and vcov consume theta summaries only", {
  fit <- re_fit()
  theta <- fit$metadata$theta_summary

  expect_equal(coef(fit), stats::setNames(theta$mean, as.character(theta$site)))
  expect_equal(coef(fit, type = "map"), stats::setNames(theta$map, as.character(theta$site)))
  err <- tryCatch(coef(fit, type = "rep"), error = identity)
  expect_s3_class(err, "bef_invalid_args")
  expect_equal(err$arg, "type")

  V <- vcov(fit)
  expect_equal(dim(V), c(nrow(theta), nrow(theta)))
  expect_equal(unname(diag(V)), theta$sd^2)
  expect_true(all(V[row(V) != col(V)] == 0))
})

test_that("RE confint returns theta and g interval data frames", {
  fit <- re_fit()
  theta_ci <- confint(fit, level = 0.8, type = "theta")

  expect_named(theta_ci, c("site", "lower", "upper", "point"))
  expect_equal(nrow(theta_ci), fit$metadata$data_list$K)
  expect_equal(theta_ci$point, fit$metadata$theta_summary$mean)
  expect_true(all(theta_ci$lower <= theta_ci$upper))

  subset_ci <- confint(fit, parm = c(2L, 4L), level = 0.8)
  expect_equal(subset_ci$site, c(2L, 4L))

  named_ci <- confint(fit, parm = c("1", "3"), level = 0.8)
  expect_equal(named_ci$site, c(1L, 3L))

  g_ci <- confint(fit, level = 0.8, type = "g")
  expect_named(g_ci, c("site", "lower", "upper", "point"))
  expect_equal(g_ci$site, c("mean_g", "var_g", "sd_g"))
  expect_equal(g_ci$point, c(
    fit$metadata$mean_g_summary$mean,
    fit$metadata$var_g_summary$mean,
    attr(fit$metadata, "sd_g_summary", exact = TRUE)$mean
  ))
  expect_true(all(g_ci$lower <= g_ci$upper))

  err <- tryCatch(confint(fit, parm = 99L), error = identity)
  expect_s3_class(err, "bef_invalid_args")
  expect_equal(err$arg, "parm")

  err <- tryCatch(confint(fit, type = "t"), error = identity)
  expect_s3_class(err, "bef_invalid_args")
  expect_equal(err$arg, "type")
})

test_that("RE as.data.frame returns theta_summary with row-name control", {
  fit <- re_fit()

  out <- as.data.frame(fit)
  expect_identical(out, fit$metadata$theta_summary)

  out <- as.data.frame(fit, row.names = paste0("row", seq_len(nobs(fit))))
  expect_equal(row.names(out), paste0("row", seq_len(nobs(fit))))

  err <- tryCatch(as.data.frame(fit, row.names = "too-short"), error = identity)
  expect_s3_class(err, "bef_invalid_args")
})

test_that("RE plot shell validates arguments and returns numeric payloads", {
  fit <- re_fit()
  payload_fun <- re_ns(".bef_plot_payload_bef_fit_re")

  caterpillar_payload <- payload_fun(fit, type = "caterpillar", level = 0.8, sort_by = "mean")
  expect_type(caterpillar_payload, "list")
  expect_equal(caterpillar_payload$type, "caterpillar")
  expect_named(
    caterpillar_payload$data,
    c(
      "site", "lower", "upper", "point", "inner_lower", "inner_upper",
      "sd", "sigma", "position"
    )
  )
  expect_equal(caterpillar_payload$data$point, sort(caterpillar_payload$data$point))
  expect_equal(caterpillar_payload$reference, fit$metadata$mean_g_summary$mean)

  payload <- payload_fun(fit, type = "g", level = 0.8, sort_by = "mean")
  expect_equal(payload$type, "g")
  expect_named(payload$data, c("kind", "grid", "lower", "upper", "point"))
  expect_equal(payload$data$kind, rep("density", fit$metadata$data_list$L))
  expect_equal(payload$data$grid, fit$metadata$data_list$grid)

  plot_out <- plot(fit, type = "caterpillar", level = 0.8)
  if (requireNamespace("ggplot2", quietly = TRUE) &&
      !identical(Sys.getenv("BAYESEFRON_NO_GGPLOT2"), "1")) {
    expect_s3_class(plot_out, "ggplot")
    expect_equal(
      attr(plot_out, "bef_plot_payload", exact = TRUE)$data,
      caterpillar_payload$data
    )
  }

  err <- tryCatch(plot(fit, level = 1), error = identity)
  expect_s3_class(err, "bef_invalid_args")

  err <- tryCatch(plot(fit, type = "c"), error = identity)
  expect_s3_class(err, "bef_invalid_args")
  expect_equal(err$arg, "type")

  err <- tryCatch(plot(fit, sort_by = "m"), error = identity)
  expect_s3_class(err, "bef_invalid_args")
  expect_equal(err$arg, "sort_by")
})
