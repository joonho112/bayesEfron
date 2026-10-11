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
  grid <- seq(-1, 1, length.out = 51L)
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
    # the mode in the table has to be a grid point
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
    stan_file_sha256 = paste(rep("d", 64L), collapse = ""),
    data_list = list(
      K = as.integer(K),
      theta_hat = seq(-0.4, 0.4, length.out = K),
      sigma = seq(0.1, 0.3, length.out = K),
      L = 51L,
      grid = grid,
      M = 3L,
      B = matrix(seq_len(51L * 3L) / 100, nrow = 51L, ncol = 3L),
      store_grid_quantities = 0L
    ),
    runtime_seconds = 2.5,
    mean_g_summary = re_summary(0),
    var_g_summary = re_summary(1),
    theta_summary = theta_summary,
    effective_params_summary = re_summary(3),
    log_marginal_likelihood_summary = re_summary(-10)
  )
  attr(metadata, "sd_g_summary") <- re_summary(1)
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

  re_ns("validate_bef_fit_re")(fixture_consistent(
    re_ns("new_bef_fit_re")(
      draws = bef_fixture_draws(
        theta_map = theta_map, theta_mean = theta_mean,
        theta_sd = theta_sd, theta_rep = theta_rep,
        scalars = list(
          mean_g = seq(-0.1, 0.1, length.out = S),
          var_g = seq(0.8, 1.2, length.out = S)
        ),
        grid_points = metadata$data_list$L
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
  ))
}

test_that("summary() of a fit adds the site table, with intervals at the level asked for", {
  fit <- re_fit()
  fit_summary <- summary(fit)

  expect_s3_class(fit_summary, "summary.bef_fit_re")
  expect_s3_class(fit_summary, "summary.bef_fit")
  expect_named(fit_summary, c("prior_summary", "diagnostics", "theta_summary"))
  table <- fit_summary$theta_summary
  expect_named(table, c("site", "label", "estimate", "se", "mean", "sd", "lower", "upper",
                        "mode", "hpdi_lower", "hpdi_upper", "map"))
  expect_identical(table, re_ns(".bef_site_table")(fit))
  expect_identical(attr(table, "level"), 0.9)
  expect_identical(table$estimate, fit$metadata$data_list$theta_hat)
  expect_identical(table$se, fit$metadata$data_list$sigma)
  expect_identical(table$label, as.character(seq_len(nrow(table))))
  # The names of version 0.3.0 hold the same values.
  expect_identical(table$hpdi_lower, table$lower)
  expect_identical(table$hpdi_upper, table$upper)
  expect_identical(table$map, table$mode)
  expect_true(any(grepl("Site effects (posterior summaries, 90% intervals):",
                        format(fit_summary), fixed = TRUE)))

  fit_summary_80 <- summary(fit, level = 0.8)
  expect_equal(attr(fit_summary_80, "level", exact = TRUE), 0.8)
  narrower <- fit_summary_80$theta_summary
  expect_identical(attr(narrower, "level"), 0.8)
  expect_true(all(narrower$lower >= table$lower & narrower$upper <= table$upper))
  expect_false(identical(narrower$lower, table$lower))
  expect_identical(narrower$hpdi_lower, narrower$lower)
  expect_equal(narrower$lower, confint(fit, level = 0.8)$lower)
  expect_true(any(grepl("80% intervals", format(fit_summary_80), fixed = TRUE)))
})

test_that("coef() returns the means or the modes of the site table, and vcov() its variances on a diagonal", {
  fit <- re_fit()
  theta <- fit$metadata$theta_summary

  expect_equal(coef(fit), stats::setNames(theta$mean, as.character(theta$site)))
  expect_equal(coef(fit, type = "mode"), stats::setNames(theta$mode, as.character(theta$site)))
  expect_identical(coef(fit, type = "map"), coef(fit, type = "mode"))
  err <- tryCatch(coef(fit, type = "rep"), error = identity)
  expect_s3_class(err, "bef_invalid_args")
  expect_equal(err$arg, "type")

  V <- vcov(fit)
  expect_equal(dim(V), c(nrow(theta), nrow(theta)))
  expect_equal(unname(diag(V)), theta$sd^2)
  expect_true(all(V[row(V) != col(V)] == 0))
})

test_that("confint() gives equal-tailed or shortest intervals for the site effects and for the moments of g", {
  fit <- re_fit()
  theta_ci <- confint(fit, level = 0.8, type = "theta")

  expect_named(theta_ci, c("site", "mean", "lower", "upper", "point"))
  expect_s3_class(theta_ci, "data.frame")
  expect_identical(attr(theta_ci, "level"), 0.8)
  expect_equal(nrow(theta_ci), fit$metadata$data_list$K)
  expect_equal(theta_ci$mean, fit$metadata$theta_summary$mean)
  expect_identical(theta_ci$point, theta_ci$mean)
  expect_true(all(theta_ci$lower <= theta_ci$upper))
  expect_output(print(theta_ci), "80% posterior intervals (equal-tailed)", fixed = TRUE)

  # The limits come from the posterior distribution averaged over the draws
  # of g; method = "draws" gives the quantiles of the theta_rep draws.
  engine <- .bef_engine(fit)
  expect_equal(unname(as.matrix(theta_ci[, c("lower", "upper")])),
               .bef_site_quantile(engine, c(0.1, 0.9)))
  from_draws <- confint(fit, level = 0.8, method = "draws")
  theta_rep <- .bef_site_draws(fit, "theta_rep")
  expect_equal(from_draws$lower,
               unname(apply(theta_rep, 2L, posterior::quantile2, probs = 0.1)))
  wider <- confint(fit, level = 0.95)
  expect_true(all(wider$lower <= theta_ci$lower & wider$upper >= theta_ci$upper))

  # The shortest interval is never longer than the equal-tailed one, and
  # holds the same probability.
  shortest <- confint(fit, level = 0.8, interval = "hpd")
  expect_true(all(shortest$upper - shortest$lower <=
                    (theta_ci$upper - theta_ci$lower) * (1 + 1e-9)))
  held <- vapply(seq_len(nrow(shortest)), function(k) {
    diff(.bef_site_cdf(engine, c(shortest$lower[k], shortest$upper[k]), sites = k)[1L, ])
  }, numeric(1))
  expect_equal(held, rep(0.8, nrow(shortest)), tolerance = 1e-10)
  expect_output(print(shortest), "shortest", fixed = TRUE)
  expect_error(confint(fit, method = "sample"), class = "bef_invalid_args")
  expect_error(confint(fit, interval = "central"), class = "bef_invalid_args")

  subset_ci <- confint(fit, parm = c(2L, 4L), level = 0.8)
  expect_equal(subset_ci$site, c(2L, 4L))

  named_ci <- confint(fit, parm = c("1", "3"), level = 0.8)
  expect_equal(named_ci$site, c(1L, 3L))

  g_ci <- confint(fit, level = 0.8, type = "g")
  expect_named(g_ci, c("parameter", "mean", "lower", "upper", "site", "point"))
  expect_equal(g_ci$parameter, c("mean_g", "var_g", "sd_g"))
  expect_identical(g_ci$site, g_ci$parameter)
  expect_identical(g_ci$point, g_ci$mean)
  g_hpd <- confint(fit, level = 0.8, type = "g", interval = "hpd")
  expect_true(all(g_hpd$upper - g_hpd$lower <= g_ci$upper - g_ci$lower + 1e-12))
  expect_equal(g_ci$mean, c(
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

test_that("as.data.frame() of a fit is the site table, and row.names needs one name for each site", {
  fit <- re_fit()

  out <- as.data.frame(fit)
  expect_identical(out, fit$metadata$theta_summary)

  out <- as.data.frame(fit, row.names = paste0("row", seq_len(nobs(fit))))
  expect_equal(row.names(out), paste0("row", seq_len(nobs(fit))))

  err <- tryCatch(as.data.frame(fit, row.names = "too-short"), error = identity)
  expect_s3_class(err, "bef_invalid_args")
})

test_that("the data of the base plots have their columns, plot() draws the same intervals, and wrong arguments are errors", {
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
  expect_s3_class(plot_out, "ggplot")
  # The numbers of the plot are those of the data of the base plot.
  expect_equal(plot_out$data$lower, caterpillar_payload$data$lower)
  expect_equal(plot_out$data$upper, caterpillar_payload$data$upper)
  expect_equal(plot_out$data$mean, caterpillar_payload$data$point)

  err <- tryCatch(plot(fit, level = 1), error = identity)
  expect_s3_class(err, "bef_invalid_args")

  err <- tryCatch(plot(fit, type = "c"), error = identity)
  expect_s3_class(err, "bef_invalid_args")
  expect_equal(err$arg, "type")

  err <- tryCatch(plot(fit, sort_by = "m"), error = identity)
  expect_s3_class(err, "bef_invalid_args")
  expect_equal(err$arg, "sort_by")
})
