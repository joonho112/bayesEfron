# plot_ppc(): the probability integral transforms are those of the marginal
# model. The bands use an independent-uniform reference.

ppc_ns <- function(name) getFromNamespace(name, "bayesEfron")

test_that("the PIT is the distribution function of the marginal model at the estimate", {
  fit <- plot_fit(12L)
  engine <- ppc_ns(".bef_engine")(fit)
  pit <- ppc_ns(".bef_pit")(engine)
  g <- exp(engine$log_g)
  direct <- vapply(seq_len(engine$K), function(i) {
    mean(g %*% stats::pnorm((engine$theta_hat[[i]] - engine$grid) / engine$sigma[[i]]))
  }, numeric(1L))
  expect_equal(pit, direct)
  expect_true(all(pit > 0 & pit < 1))

  # With normal components their variance is added to that of the estimate.
  normal <- ppc_ns(".bef_engine")(plot_fit(12L, component = "normal"))
  expect_gt(normal$tau_c, 0)
  g <- exp(normal$log_g)
  direct <- vapply(seq_len(normal$K), function(i) {
    mean(g %*% stats::pnorm((normal$theta_hat[[i]] - normal$grid) /
                              sqrt(normal$sigma[[i]]^2 + normal$tau_c^2)))
  }, numeric(1L))
  expect_equal(ppc_ns(".bef_pit")(normal), direct)
})

test_that("the PIT values of data drawn from a known g are uniform", {
  # A known g on a grid, and estimates drawn from the model that it defines.
  grid <- seq(-3, 4, length.out = 101L)
  g <- 0.4 * stats::dnorm(grid, -1, 0.4) + 0.6 * stats::dnorm(grid, 1.5, 0.7)
  g <- g / sum(g)
  K <- 2000L
  data <- withr::with_seed(31L, {
    sigma <- stats::runif(K, 0.3, 1.2)
    theta <- sample(grid, K, replace = TRUE, prob = g)
    list(theta_hat = theta + stats::rnorm(K, 0, sigma), sigma = sigma)
  })
  engine <- ppc_ns(".bef_engine_from_parts")(
    log_g = matrix(log(g), nrow = 4L, ncol = 101L, byrow = TRUE), grid = grid,
    theta_hat = data$theta_hat, sigma = data$sigma
  )
  pit <- ppc_ns(".bef_pit")(engine)
  expect_gt(stats::ks.test(pit, "punif")$p.value, 0.01)
  expect_equal(mean(pit), 0.5, tolerance = 0.03)
  # Estimates with standard errors that are too small are found out.
  engine$sigma <- data$sigma / 2
  expect_lt(stats::ks.test(ppc_ns(".bef_pit")(engine), "punif")$p.value, 1e-3)
})

test_that("the difference of the distribution functions and its band are drawn", {
  fit <- plot_fit(150L)
  pit <- ppc_ns(".bef_pit")(ppc_ns(".bef_engine")(fit))
  plot <- plot_ppc(fit)
  expect_s3_class(plot, "ggplot")
  expect_identical(unname(plot_geoms(plot)), c("GeomRibbon", "GeomHline", "GeomStep"))
  expect_null(plot$labels$title)
  expect_identical(plot$labels$x, "PIT of the estimate")
  line <- plot_layer_data(plot, 3L)
  expect_equal(line$y, stats::ecdf(pit)(line$x) - line$x)
  expect_equal(range(line$x), c(0, 1))
  band <- plot_layer_data(plot, 1L)
  expect_equal(band$ymax, stats::qbinom(0.975, 150L, band$x) / 150 - band$x)
  expect_equal(band$ymin, stats::qbinom(0.025, 150L, band$x) / 150 - band$x)
  expect_match(plot_caption(plot), "pointwise 95% reference intervals for independent uniform values")
  expect_match(plot_caption(plot), "does not account for fitting g or dependence among PIT values")
  expect_no_match(plot_caption(plot), "if the model holds")

  narrow_plot <- plot_ppc(fit, level = 0.5)
  narrow <- plot_layer_data(narrow_plot, 1L)
  expect_true(all(narrow$ymax <= band$ymax))
  expect_equal(narrow$ymax, stats::qbinom(0.75, 150L, narrow$x) / 150 - narrow$x)
  expect_equal(narrow$ymin, stats::qbinom(0.25, 150L, narrow$x) / 150 - narrow$x)
  expect_match(plot_caption(narrow_plot), "pointwise 50% reference intervals")
})

test_that("the histogram counts every site, and few sites get the other plot", {
  fit <- plot_fit(150L)
  plot <- plot_ppc(fit, type = "pit_hist")
  expect_identical(unname(plot_geoms(plot)), c("GeomRect", "GeomRect", "GeomHline"))
  bars <- plot_layer_data(plot, 2L)
  expect_identical(sum(bars$ymax), 150)
  expect_equal(plot_layer_data(plot, 3L)$yintercept, 150 / nrow(bars))
  band <- plot_layer_data(plot, 1L)
  expect_equal(band$ymax, stats::qbinom(0.975, 150L, 1 / nrow(bars)))
  expect_equal(band$ymin, stats::qbinom(0.025, 150L, 1 / nrow(bars)))
  expect_match(plot_caption(plot), "pointwise 95% reference intervals for bin counts from independent uniform values")
  expect_match(plot_caption(plot), "does not account for fitting g or dependence among PIT values")
  narrow_plot <- plot_ppc(fit, type = "pit_hist", level = 0.8)
  narrow <- plot_layer_data(narrow_plot, 1L)
  expect_equal(narrow$ymax, stats::qbinom(0.9, 150L, 1 / nrow(bars)))
  expect_equal(narrow$ymin, stats::qbinom(0.1, 150L, 1 / nrow(bars)))
  expect_match(plot_caption(narrow_plot), "pointwise 80% reference intervals")

  few <- plot_fit(12L)
  expect_message(fallback <- plot_ppc(few, type = "pit_hist"), "With 12 sites",
                 class = "bef_plot_few_sites")
  expect_identical(unname(plot_geoms(fallback)), c("GeomRibbon", "GeomHline", "GeomStep"))
  expect_no_message(plot_ppc(few))
})

test_that("the marginal density of an estimate is drawn over the estimates", {
  fit <- plot_fit(150L)
  plot <- plot_ppc(fit, type = "marginal")
  expect_identical(unname(plot_geoms(plot)), c("GeomRect", "GeomLine"))
  expect_identical(plot$labels$x, "Estimate")
  curve <- plot_layer_data(plot, 2L)
  expect_equal(sum(curve$y) * diff(curve$x[1:2]), 1, tolerance = 0.02)
  bars <- plot_layer_data(plot, 1L)
  expect_equal(sum(bars$ymax * (bars$xmax - bars$xmin)), 1)

  scaled <- plot_layer_data(plot_ppc(fit, type = "marginal", scale = 10), 2L)
  expect_equal(scaled$x, 10 * curve$x)
  expect_equal(scaled$y, curve$y / 10)
})

test_that("the PIT can be computed without the site", {
  skip_if_not_installed("loo")
  fit <- plot_fit(150L)
  engine <- ppc_ns(".bef_engine")(fit)
  all_sites <- ppc_ns(".bef_pit")(engine)
  left_out <- ppc_ns(".bef_pit")(engine, loo = TRUE)
  expect_length(left_out, 150L)
  expect_true(all(left_out >= 0 & left_out <= 1))
  expect_false(isTRUE(all.equal(all_sites, left_out)))
  # One of 150 sites moves the model little.
  expect_lt(max(abs(all_sites - left_out)), 0.05)
  # The blocks of sites give what one block gives.
  expect_equal(ppc_ns(".bef_pit")(engine, loo = TRUE, max_cells = 4000), left_out)

  plot <- plot_ppc(fit, loo = TRUE)
  expect_identical(plot$labels$x, "PIT of the estimate, computed without its site")
  expect_match(plot_caption(plot), "posterior distribution without its site")
  expect_match(plot_caption(plot), "independent uniform values")
  expect_match(plot_caption(plot), "dependence among PIT values")
})

test_that("wrong arguments are errors", {
  fit <- plot_fit(12L)
  expect_error(plot_ppc(fit, type = "ecdf"), class = "bef_invalid_args")
  expect_error(plot_ppc(fit, loo = "yes"), class = "bef_invalid_args")
  expect_error(plot_ppc(fit, level = 95), class = "bef_invalid_args")
  expect_error(plot_ppc(fit, bins = 10), "no argument", class = "bef_invalid_args")
  expect_error(plot_ppc(NULL), class = "bef_invalid_fit")
})
