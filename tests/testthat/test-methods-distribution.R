# The methods that return the posterior distribution of the site effects:
# vcov(full = TRUE), quantile(), simulate() and residuals().

distribution_fit <- function() {
  readRDS(test_path("_fixtures", "short_fit.rds"))
}

test_that("the full covariance matrix has the posterior variances on its diagonal", {
  fit <- distribution_fit()
  table <- as.data.frame(fit)
  V <- vcov(fit, full = TRUE)
  expect_equal(dim(V), c(5L, 5L))
  expect_equal(unname(diag(V)), table$sd^2)
  expect_equal(V, t(V))
  expect_identical(rownames(V), table$label)
  expect_true(all(eigen(V, symmetric = TRUE, only.values = TRUE)$values > -1e-12))
  expect_equal(diag(vcov(fit)), diag(V))
  expect_true(all(vcov(fit)[upper.tri(V)] == 0))
  # The covariances are those of the means of the effects given g.
  conditional <- .bef_site_conditional(.bef_engine(fit))
  centred <- sweep(conditional$mean, 2L, colMeans(conditional$mean))
  expect_equal(unname(V[1L, 4L]), mean(centred[, 1L] * centred[, 4L]))
  expect_error(vcov(fit, full = NA), class = "bef_invalid_args")
})

test_that("quantile() returns a labelled matrix that agrees with confint()", {
  fit <- distribution_fit()
  q <- quantile(fit)
  expect_equal(dim(q), c(5L, 5L))
  expect_identical(colnames(q), c("5%", "25%", "50%", "75%", "95%"))
  expect_identical(rownames(q), as.character(1:5))
  expect_true(all(apply(q, 1L, diff) >= 0))
  intervals <- confint(fit)
  expect_equal(unname(q[, "5%"]), intervals$lower)
  expect_equal(unname(q[, "95%"]), intervals$upper)
  step <- quantile(fit, probs = 0.5, reading = "step")
  expect_true(all(step[, 1L] %in% fit$metadata$data_list$grid))
  expect_error(quantile(fit, probs = 1.2), class = "bef_invalid_args")
})

test_that("residuals() are the estimates less the posterior means", {
  fit <- distribution_fit()
  table <- as.data.frame(fit)
  expect_equal(residuals(fit), stats::setNames(table$estimate - table$mean, table$label))
})

test_that("simulate() draws the effects jointly, with the moments of the posterior", {
  fit <- distribution_fit()
  draws <- simulate(fit, nsim = 4000L, seed = 11L)
  expect_equal(dim(draws), c(4000L, 5L))
  expect_identical(colnames(draws), as.character(1:5))
  expect_true(all(draws %in% fit$metadata$data_list$grid))
  table <- as.data.frame(fit)
  # Means and variances agree with the engine within Monte Carlo error.
  se_mean <- table$sd / sqrt(4000)
  expect_true(all(abs(colMeans(draws) - table$mean) < 4 * se_mean))
  expect_equal(unname(apply(draws, 2L, stats::var)), table$sd^2, tolerance = 0.15)
  # The covariance of two sites: within four standard errors of a sample
  # covariance of 4,000 draws.
  V <- vcov(fit, full = TRUE)
  se_cov <- sqrt((V[1L, 1L] * V[5L, 5L] + V[1L, 5L]^2) / 4000)
  expect_lt(abs(stats::cov(draws[, 1L], draws[, 5L]) - V[1L, 5L]), 4 * se_cov)
  # The seed makes the draws repeatable and leaves the generator alone.
  set.seed(3)
  before <- stats::runif(1L)
  set.seed(3)
  again <- simulate(fit, nsim = 4000L, seed = 11L)
  expect_identical(stats::runif(1L), before)
  expect_identical(again, draws)
  expect_equal(dim(simulate(fit)), c(1L, 5L))
  expect_error(simulate(fit, nsim = 0), class = "bef_invalid_args")

  # With normal components the draws are not confined to the grid.
  engine <- .bef_engine(fit)
  engine$tau_c <- 2 / 3 * engine$spacing
  withr::local_seed(5)
  smooth <- .bef_simulate_sites(engine, 2000L)
  moments <- .bef_site_moments(engine)
  expect_false(all(smooth %in% engine$grid))
  expect_true(all(abs(colMeans(smooth) - moments$mean) < 4 * sqrt(moments$var / 2000)))
})
