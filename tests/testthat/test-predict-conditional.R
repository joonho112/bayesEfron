# Independent predictive moment and Bayes calculations on stored draws.

test_that("different future sites share g but have separate true effects", {
  for (component in c("point", "normal")) {
    fit <- plot_fit(12L, component = component)
    data <- fit$metadata$data_list
    dm <- posterior::as_draws_matrix(fit$draws)
    g <- exp(dm[, sprintf("log_g[%d]", seq_len(data$L)), drop = FALSE])
    g <- g / rowSums(g)
    means <- as.numeric(g %*% data$grid)
    within <- rowSums(g * outer(means, data$grid, function(m, x) (x - m)^2)) + data$tau_c^2
    expected_cov <- mean((means - mean(means))^2)
    for (jitter in c(FALSE, TRUE)) {
      se <- c(first = 0.001, second = 0.002)
      out <- predict(fit, newdata = se, type = "theta_hat", n = 30000L,
                     jitter = jitter, seed = 791L)
      # Conditional independence gives covariance Var(E[theta | g]), while
      # Var(theta1 - theta2) removes that shared uncertainty altogether.
      products <- (out[, 1L] - mean(means)) * (out[, 2L] - mean(means))
      expect_lt(abs(mean(products) - expected_cov), 6 * sd(products) / sqrt(nrow(out)))
      extra <- if (component == "point" && jitter) diff(data$grid)[1L]^2 / 12 else 0
      expected_difference_var <- 2 * (mean(within) + extra) + sum(se^2)
      squares <- (out[, 1L] - out[, 2L])^2
      expect_lt(abs(mean(squares) - expected_difference_var), 6 * sd(squares) / sqrt(nrow(out)))
    }
  }
})

test_that("several new sites give a column each, the first as for that site alone", {
  for (component in c("point", "normal")) {
    fit <- plot_fit(12L, component = component)
    se <- c(a = 0.1, b = 0.2, c = 0.4, d = 0.8)
    out <- predict(fit, newdata = se, type = "theta_hat", n = 5000L, seed = 12L)
    expect_identical(dim(out), c(5000L, 4L))
    expect_identical(colnames(out), names(se))
    expect_identical(out[, 1L], predict(fit, newdata = se[[1L]], type = "theta_hat",
                                        n = 5000L, seed = 12L))
    expect_null(colnames(predict(fit, newdata = unname(se), type = "theta_hat",
                                 n = 10L, seed = 12L)))
    # The effects have one distribution, so the variance of a column is
    # that of the effect and of its sampling error.
    theta <- predict(fit, n = 5000L, seed = 12L)
    spread <- apply(out, 2L, var)
    expect_true(all(diff(spread) > 0))
    expect_lt(max(abs(spread - (var(theta) + se^2)) / (var(theta) + se^2)), 0.15)
    # Different sites of a row do not share an effect: the columns are
    # nearly uncorrelated, far from the correlation that a shared effect gives.
    shared <- var(theta) / sqrt((var(theta) + se[[2L]]^2) * (var(theta) + se[[3L]]^2))
    expect_lt(abs(cor(out[, 2L], out[, 3L])), shared / 3)
  }
})

test_that("known new estimates condition the posterior mean g by Bayes' rule", {
  for (component in c("point", "normal")) {
    fit <- plot_fit(12L, component = component)
    data <- fit$metadata$data_list
    dm <- posterior::as_draws_matrix(fit$draws)
    g <- exp(dm[, sprintf("log_g[%d]", seq_len(data$L)), drop = FALSE])
    g <- g / rowSums(g)
    prior <- colMeans(g)
    y <- max(data$theta_hat) + mean(data$sigma)
    se <- median(data$sigma)
    level <- 0.8
    got <- predict(fit, newdata = data.frame(estimate = y, se = se), level = level)
    likelihood <- dnorm(y, data$grid, sqrt(se^2 + data$tau_c^2))
    mass <- prior * likelihood
    mass <- mass / sum(mass)
    shrink <- data$tau_c^2 / (data$tau_c^2 + se^2)
    centers <- data$grid + shrink * (y - data$grid)
    m <- sum(mass * centers)
    v <- sum(mass * (centers - m)^2) + shrink * se^2
    expect_equal(got$mean, m, tolerance = 1e-12)
    expect_equal(got$sd, sqrt(v), tolerance = 1e-12)
    if (component == "point") {
      probs <- c((1 - level) / 2, (1 + level) / 2)
      cumulative <- cumsum(mass)
      expected <- vapply(probs, function(p) {
        j <- which(cumulative >= p)[1L]
        before <- if (j > 1L) cumulative[j - 1L] else 0
        data$grid[j] + diff(data$grid)[1L] * ((p - before) / mass[j] - 0.5)
      }, numeric(1L))
      expect_equal(c(got$lower, got$upper), expected, tolerance = 1e-12)
      expect_equal(got$mode, data$grid[which.max(mass)])
    } else {
      cdf <- function(q) sum(mass * pnorm(q, centers, sqrt(shrink) * se))
      expect_equal(c(cdf(got$lower), cdf(got$upper)), c(0.1, 0.9), tolerance = 1e-10)
    }
    # A fixed posterior-mean mixing distribution supplies an independent
    # implementation of the same conditional predictive calculation.
    if (requireNamespace("ebnm", quietly = TRUE)) {
      fixed <- ebnm::ebnm_npmle(y, se, g_init = as_normalmix(fit), fix_g = TRUE)
      expect_equal(got$mean, fixed$posterior$mean, tolerance = 1e-10)
      expect_equal(got$sd, fixed$posterior$sd, tolerance = 1e-10)
    }
  }
})
