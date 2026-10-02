# The table of site effects: the posterior standard deviation and the
# posterior mode, which are computed in R from the draws (Lee and Sui 2025,
# Section 4.3 and Appendix A).

site_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

test_that("the posterior standard deviation adds the two parts of the variance", {
  posterior_sd <- site_ns(".bef_posterior_sd")
  theta_mean <- cbind(c(0, 0.2, 0.4), c(1, 1, 1))
  theta_sd <- cbind(c(0.3, 0.4, 0.5), c(0.2, 0.2, 0.2))

  expect_equal(
    posterior_sd(theta_mean, theta_sd),
    sqrt(colMeans(theta_sd^2) + apply(theta_mean, 2L, stats::var))
  )
  # When the mean given g varies over the draws, the posterior standard
  # deviation exceeds the average of the standard deviations given g.
  expect_gt(posterior_sd(theta_mean, theta_sd)[1L], mean(theta_sd[, 1L]))
  # When neither varies, it is the standard deviation given g.
  expect_equal(posterior_sd(theta_mean, theta_sd)[2L], 0.2)
  # One draw: the standard deviation given that g.
  expect_equal(
    posterior_sd(theta_mean[1L, , drop = FALSE], theta_sd[1L, , drop = FALSE]),
    c(0.3, 0.2)
  )
})

test_that("the posterior probabilities average Bayes' rule over the draws of g", {
  probabilities <- site_ns(".bef_posterior_probabilities")
  withr::local_seed(42)
  grid <- seq(-2, 2, length.out = 21L)
  n_draws <- 7L
  theta_hat <- c(-1.1, 0.2, 1.4)
  sigma <- c(0.3, 0.6, 0.25)
  log_g <- matrix(stats::rnorm(n_draws * length(grid)), nrow = n_draws)
  log_g <- log_g - log(rowSums(exp(log_g)))

  fast <- probabilities(log_g, theta_hat, sigma, grid)
  one_by_one <- t(vapply(
    seq_along(theta_hat),
    function(site) {
      total <- numeric(length(grid))
      for (draw in seq_len(n_draws)) {
        weights <- exp(log_g[draw, ]) *
          stats::dnorm(theta_hat[site], mean = grid, sd = sigma[site])
        total <- total + weights / sum(weights)
      }
      total / n_draws
    },
    numeric(length(grid))
  ))

  expect_equal(fast, one_by_one, tolerance = 1e-12)
  expect_equal(rowSums(fast), rep(1, 3L))
  # The sites are handled in blocks; one site per block gives the same.
  expect_equal(probabilities(log_g, theta_hat, sigma, grid, max_cells = 1), fast)
})

test_that("the posterior mode is the mode of the averaged probabilities", {
  posterior_mode <- site_ns(".bef_posterior_mode")
  grid <- c(-1, 0, 1)
  # A nearly flat likelihood, so that the posterior of the site is close to g.
  # The modes given g are -1, 1 and 1, whose average, 1/3, is not a grid
  # point; the averaged probabilities are largest at 1.
  g <- rbind(c(0.8, 0.1, 0.1), c(0.1, 0.1, 0.8), c(0.1, 0.1, 0.8))
  expect_identical(posterior_mode(log(g), theta_hat = 0, sigma = 50, grid = grid), 1)

  # Equal probabilities: the smallest grid point.
  expect_identical(
    posterior_mode(log(rbind(c(0.5, 0.5))), theta_hat = 0, sigma = 1, grid = c(-1, 1)),
    -1
  )
})

test_that("the posterior probabilities stay finite when a product underflows", {
  probabilities <- site_ns(".bef_posterior_probabilities")
  grid <- c(-5, -2.5, 0, 2.5, 5)
  # g is a point mass at 5 to machine precision. The first site has a precise
  # estimate at -5, where g underflows to zero.
  log_g <- matrix(-800, nrow = 3L, ncol = 5L)
  log_g[, 5L] <- 0
  out <- probabilities(log_g, theta_hat = c(-5, 5), sigma = c(0.01, 0.5), grid = grid)

  expect_true(all(is.finite(out)))
  expect_equal(rowSums(out), c(1, 1))
  expect_equal(out[1L, ], c(1, 0, 0, 0, 0))
  expect_equal(out[2L, ], c(0, 0, 0, 0, 1))
})

test_that("the stored example fit carries the summaries that the code computes", {
  fit <- bayesEfron::raudenbush_fit
  data <- fit$metadata$data_list
  table <- as.data.frame(fit)

  expect_equal(
    table$sd,
    site_ns(".bef_posterior_sd")(
      site_ns(".bef_site_draws")(fit, "theta_mean"),
      site_ns(".bef_site_draws")(fit, "theta_sd")
    )
  )
  log_g <- log(site_ns(".bef_plain_draw_matrix")(site_ns(".bef_grid_draws")(fit, "g")))
  expect_equal(
    table$map,
    site_ns(".bef_posterior_mode")(log_g, data$theta_hat, data$sigma, data$grid)
  )
  expect_true(all(table$map %in% data$grid))
  expect_equal(unname(diag(vcov(fit))), table$sd^2)
  expect_equal(unname(coef(fit, type = "map")), table$map)

  # The posterior standard deviation is at least the average standard
  # deviation given g, and each interval contains the posterior mode.
  expect_true(all(table$sd >= colMeans(site_ns(".bef_site_draws")(fit, "theta_sd"))))
  expect_true(all(table$map >= table$hpdi_lower & table$map <= table$hpdi_upper))
  expect_identical(site_ns("validate_bef_fit_re")(fit), fit)
})
