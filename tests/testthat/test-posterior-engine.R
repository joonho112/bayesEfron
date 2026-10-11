# The posterior engine: site summaries computed from the draws of log_g.

# Site weights, conditional moments and log likelihood computed draw by
# draw, without any of the engine's code.
brute_force <- function(log_g, grid, theta_hat, sigma, tau_c = 0) {
  S <- nrow(log_g)
  K <- length(theta_hat)
  L <- length(grid)
  weights <- matrix(0, K, L)
  cond_mean <- cond_var <- log_lik <- matrix(NA_real_, S, K)
  for (i in seq_len(K)) {
    ll <- dnorm(theta_hat[i], grid, sqrt(sigma[i]^2 + tau_c^2), log = TRUE)
    b <- tau_c^2 / (tau_c^2 + sigma[i]^2)
    centre <- grid + b * (theta_hat[i] - grid)
    for (s in seq_len(S)) {
      lp <- log_g[s, ] + ll
      top <- max(lp)
      w <- exp(lp - top)
      log_lik[s, i] <- top + log(sum(w))
      w <- w / sum(w)
      weights[i, ] <- weights[i, ] + w / S
      m <- sum(w * centre)
      cond_mean[s, i] <- m
      cond_var[s, i] <- sum(w * (centre - m)^2) + b * sigma[i]^2
    }
  }
  list(weights = weights, mean = cond_mean, var = cond_var, log_lik = log_lik)
}

random_case <- function(S = 40L, K = 7L, L = 31L, seed = 1L) {
  withr::local_seed(seed)
  grid <- seq(-2, 3, length.out = L)
  raw <- matrix(rnorm(S * L), S, L) + matrix(-0.5 * (grid - 0.4)^2, S, L, byrow = TRUE)
  log_g <- raw - log(rowSums(exp(raw)))
  list(log_g = log_g, grid = grid, theta_hat = runif(K, -1.5, 2.5),
       sigma = runif(K, 0.1, 0.9))
}

test_that("the averaged weights, the moments and the log likelihood match a draw-by-draw computation", {
  case <- random_case()
  for (tau_c in c(0, 0.11)) {
    engine <- .bef_engine_from_parts(case$log_g, case$grid, case$theta_hat, case$sigma, tau_c)
    truth <- brute_force(case$log_g, case$grid, case$theta_hat, case$sigma, tau_c)
    weights <- .bef_site_weights(engine)
    expect_equal(weights, truth$weights, tolerance = 1e-12)
    expect_equal(rowSums(weights), rep(1, engine$K), tolerance = 1e-12)
    conditional <- .bef_site_conditional(engine)
    expect_equal(conditional$mean, truth$mean, tolerance = 1e-12)
    expect_equal(conditional$var, truth$var, tolerance = 1e-10)
    expect_equal(.bef_log_lik(engine), truth$log_lik, tolerance = 1e-12)
  }
})

test_that("for point masses the weights are the posterior probabilities of version 0.3.0", {
  case <- random_case(seed = 2L)
  engine <- .bef_engine_from_parts(case$log_g, case$grid, case$theta_hat, case$sigma)
  expect_equal(
    .bef_site_weights(engine),
    .bef_posterior_probabilities(case$log_g, case$theta_hat, case$sigma, case$grid),
    tolerance = 1e-14
  )
})

test_that("the posterior variance is the within part plus the between part", {
  case <- random_case(seed = 3L)
  for (tau_c in c(0, 0.2)) {
    engine <- .bef_engine_from_parts(case$log_g, case$grid, case$theta_hat, case$sigma, tau_c)
    moments <- .bef_site_moments(engine)
    truth <- brute_force(case$log_g, case$grid, case$theta_hat, case$sigma, tau_c)
    expect_equal(moments$mean, colMeans(truth$mean), tolerance = 1e-12)
    expect_equal(moments$var_within, colMeans(truth$var), tolerance = 1e-10)
    expect_equal(moments$var_between,
                 colMeans(sweep(truth$mean, 2L, colMeans(truth$mean))^2),
                 tolerance = 1e-10)
    expect_equal(moments$var, moments$var_within + moments$var_between)
    # The same variance from the averaged mixture itself.
    weights <- .bef_site_weights(engine)
    b <- tau_c^2 / (tau_c^2 + case$sigma^2)
    centres <- outer(1 - b, case$grid) + b * case$theta_hat
    second <- rowSums(weights * (centres^2 + b * case$sigma^2))
    expect_equal(moments$var, second - moments$mean^2, tolerance = 1e-9)
  }
})

test_that("the results do not depend on how the sites are split into blocks", {
  case <- random_case(S = 30L, K = 11L, seed = 4L)
  engine <- .bef_engine_from_parts(case$log_g, case$grid, case$theta_hat, case$sigma, 0.05)
  expect_equal(.bef_site_weights(engine, max_cells = 60), .bef_site_weights(engine))
  expect_equal(.bef_site_moments(engine, max_cells = 60), .bef_site_moments(engine))
  expect_equal(.bef_log_lik(engine, max_cells = 60), .bef_log_lik(engine))
  expect_equal(.bef_site_weights(engine, sites = c(9L, 2L)),
               .bef_site_weights(engine)[c(9L, 2L), ])
  expect_equal(.bef_log_lik(engine, sites = 5L), .bef_log_lik(engine)[, 5L, drop = FALSE])
})

test_that("a site whose likelihood underflows where g has its mass is computed on the log scale", {
  grid <- seq(-50, 50, by = 1)
  log_g <- matrix(-800, 3L, length(grid))
  log_g[, 6L] <- 0
  log_g <- log_g - log(rowSums(exp(log_g)))
  engine <- .bef_engine_from_parts(log_g, grid, theta_hat = c(45, -45), sigma = c(1, 1))
  weights <- .bef_site_weights(engine)
  expect_true(all(is.finite(weights)))
  expect_equal(rowSums(weights), c(1, 1))
  moments <- .bef_site_moments(engine)
  expect_true(all(is.finite(as.matrix(moments))))
  # log_g is -800 away from the one point with mass, and the estimate of the
  # first site is 90 standard errors from that point: its posterior sits on
  # the grid point at its estimate.
  expect_equal(moments$mean, c(45, -45), tolerance = 1e-8)
  expect_true(all(is.finite(.bef_log_lik(engine))))
})

test_that("the distribution function is monotone and the quantiles invert it", {
  case <- random_case(seed = 5L)
  probs <- c(0.05, 0.25, 0.5, 0.9, 0.975)
  q <- seq(-3, 4, length.out = 57L)
  rows <- seq_along(case$theta_hat)
  for (tau_c in c(0, 0.11)) {
    engine <- .bef_engine_from_parts(case$log_g, case$grid, case$theta_hat, case$sigma, tau_c)
    weights <- .bef_site_weights(engine)
    cdf <- .bef_site_cdf(engine, q, weights = weights)
    expect_true(all(apply(cdf, 1L, diff) >= -1e-12))
    expect_equal(cdf[, 1L], rep(0, engine$K), tolerance = 1e-8)
    expect_equal(cdf[, length(q)], rep(1, engine$K), tolerance = 1e-8)
    quantiles <- .bef_site_quantile(engine, probs, weights = weights)
    expect_true(all(apply(quantiles, 1L, diff) >= 0))
    back <- vapply(seq_along(probs), function(k) {
      vapply(rows, function(i) {
        .bef_site_cdf(engine, quantiles[i, k], weights = weights[i, , drop = FALSE],
                      sites = i)[1L, 1L]
      }, numeric(1L))
    }, numeric(engine$K))
    expect_equal(back, matrix(probs, engine$K, length(probs), byrow = TRUE),
                 tolerance = 1e-8)
    expect_equal(.bef_site_exceedance(engine, 0.3, weights = weights),
                 1 - .bef_site_cdf(engine, 0.3, weights = weights)[, 1L])
  }
})

test_that("point masses are read as steps or as cells of one spacing", {
  grid <- c(0, 1, 2, 3)
  log_g <- matrix(log(c(0.1, 0.4, 0.3, 0.2)), 1L)
  # A standard error so large that the posterior of the site is g itself.
  engine <- .bef_engine_from_parts(log_g, grid, theta_hat = 1.5, sigma = 1e6)
  weights <- .bef_site_weights(engine)
  expect_equal(as.numeric(weights), c(0.1, 0.4, 0.3, 0.2), tolerance = 1e-8)

  step <- .bef_site_quantile(engine, c(0.05, 0.1, 0.5, 0.95), weights = weights, reading = "step")
  expect_equal(as.numeric(step), c(0, 0, 1, 3))
  # Interpolated: the cell of point 0 runs from -0.5 to 0.5 and holds the
  # cumulative probability up to 0.1; that of point 1 runs from 0.5 to 1.5
  # and holds it from 0.1 to 0.5.
  cells <- .bef_site_quantile(engine, c(0.025, 0.05, 0.3, 0.5, 0.9), weights = weights,
                              reading = "interpolated")
  expect_equal(as.numeric(cells), c(-0.25, 0, 1, 1.5, 3), tolerance = 1e-6)

  expect_equal(.bef_site_cdf(engine, c(1, 1.5), weights = weights, reading = "step")[1L, ],
               c(0.5, 0.5), tolerance = 1e-8)
  expect_equal(.bef_site_cdf(engine, c(1, 1.5), weights = weights, reading = "interpolated")[1L, ],
               c(0.3, 0.5), tolerance = 1e-8)
  # As steps a grid point equal to the threshold is not above it.
  expect_equal(.bef_site_exceedance(engine, 1, weights = weights, reading = "step"), 0.5,
               tolerance = 1e-8)
  expect_equal(.bef_site_exceedance(engine, 1, weights = weights, reading = "interpolated"),
               0.7, tolerance = 1e-8)
})

test_that("the shortest interval for point masses is the shortest of its reading", {
  withr::local_seed(8131)
  for (i in seq_len(30L)) {
    L <- sample(c(2L, 5L, 21L, 101L), 1L)
    grid <- seq(-1, 1, length.out = L) * sample(c(0.1, 1, 50), 1L)
    mass <- rexp(L)^sample(c(1, 5), 1L)
    if (i %% 3L == 0L) mass[sample(L, ceiling(L / 3))] <- 0
    if (sum(mass) == 0) mass[1L] <- 1
    mass <- mass / sum(mass)
    level <- sample(c(.2, .5, .9), 1L)
    spacing <- diff(grid)[1L]
    # On the points: every run of grid points is tried.
    cumulative <- c(0, cumsum(mass))
    runs <- expand.grid(first = seq_len(L), last = seq_len(L))
    runs <- runs[runs$last >= runs$first, ]
    runs$held <- cumulative[runs$last + 1L] - cumulative[runs$first]
    runs <- runs[runs$held >= level - 1e-12, ]
    fewest <- min(runs$last - runs$first)
    step <- .bef_point_mass_hpd_step(grid, mass, level)
    expect_equal(diff(step), fewest * spacing)
    expect_gte(sum(mass[grid >= step[1L] - 1e-9 & grid <= step[2L] + 1e-9]), level - 1e-12)
    # Over the cells: the interval holds the level, and no interval between
    # two values of a fine mesh that holds it is shorter by more than two
    # of its steps.
    cells <- .bef_point_mass_hpd_cells(grid, spacing, mass, level)
    x <- seq(grid[1L] - spacing / 2, grid[L] + spacing / 2, length.out = 20001L)
    cdf <- vapply(x, function(value) {
      sum(mass * pmin(pmax((value - grid) / spacing + 0.5, 0), 1))
    }, numeric(1))
    at <- function(value) sum(mass * pmin(pmax((value - grid) / spacing + 0.5, 0), 1))
    expect_equal(at(cells[2L]) - at(cells[1L]), level, tolerance = 1e-10)
    mesh <- .bef_shortest_interval(x, cdf, level)
    expect_lte(diff(cells), diff(mesh) + 1e-9 * spacing)
    expect_gte(diff(cells), diff(mesh) - 2 * (x[2L] - x[1L]) - 1e-9 * spacing)
  }
})

test_that("the shortest interval for point masses handles cells without probability", {
  # The three middle points hold 0.8, with an empty cell on each side. The
  # sums 0.1 + 0.8 and 0.9 - 0.8 are not the values of the distribution
  # function at the edges, by rounding.
  mass <- c(0.1, 0, 0.3, 0.2, 0.3, 0, 0.1)
  expect_equal(.bef_point_mass_hpd_cells(0:6, 1, mass, 0.8), c(1.5, 4.5))
  expect_equal(.bef_point_mass_hpd_step(0:6, mass, 0.8), c(2, 4))
  # The same anywhere and on any scale.
  expect_equal(.bef_point_mass_hpd_cells(1e6 + 0:6 * 1e-3, 1e-3, mass, 0.8),
               1e6 + c(1.5, 4.5) * 1e-3)
  # All the probability on one point.
  expect_equal(.bef_point_mass_hpd_cells(0:4, 1, c(0, 0, 1, 0, 0), 0.5), c(1.75, 2.25))
  expect_equal(.bef_point_mass_hpd_step(0:4, c(0, 0, 1, 0, 0), 0.5), c(2, 2))
})

test_that("among runs of grid points of one length the rules pick one", {
  grid <- 1:4
  # The run with the largest probability.
  expect_equal(.bef_point_mass_hpd_step(grid, c(.25, .3, .3, .15), .5), c(2, 3))
  # Three runs with the same probability: the middle one.
  expect_equal(.bef_point_mass_hpd_step(grid, c(.25, .25, .25, .25), .5), c(2, 3))
  # Two runs with the same probability: the lower one.
  expect_equal(.bef_point_mass_hpd_step(grid, c(.3, .2, .2, .3), .5), c(1, 2))
  # A level that the probabilities cannot reach by rounding: all the points.
  expect_equal(.bef_point_mass_hpd_step(grid, rep(.25, 4), 1 + 1e-9), c(1, 4))
})

test_that("a symmetric posterior has a symmetric shortest interval", {
  grid <- seq(-1, 1, length.out = 101L)
  engine <- .bef_engine_from_parts(matrix(stats::dnorm(grid, 0, 0.3, log = TRUE), nrow = 1L),
                                   grid, theta_hat = 0, sigma = 0.2)
  for (level in c(.5, .9)) {
    limits <- .bef_site_hpd(engine, level, reading = "interpolated")[1L, ]
    expect_equal(unname(limits[1L]), -unname(limits[2L]), tolerance = 1e-9)
    # On the points a run with an even number of points cannot be
    # symmetric about a grid point: it is one of two that mirror each other.
    points <- .bef_site_hpd(engine, level, reading = "step")[1L, ]
    expect_lte(abs(sum(points)), engine$spacing + 1e-9)
  }
})

test_that("the engine reproduces the site quantities that Stan wrote", {
  # A fit of version 0.2 with the site-level draws.
  fit <- readRDS(test_path("_fixtures", "short_fit.rds"))
  engine <- .bef_engine(fit)
  conditional <- .bef_site_conditional(engine)
  # CmdStan writes log_g with a fixed number of significant digits, so a mean
  # computed again from it is off by that rounding times the scale of the
  # effects, however small the mean itself is. The error of a mean is
  # therefore measured in posterior standard deviations and that of the log
  # likelihood in absolute terms; a relative error would blow up for values
  # near zero.
  relative <- function(a, b) max(abs(a - b) / abs(b))
  stan_sd <- unname(.bef_site_draws(fit, "theta_sd"))
  expect_lt(max(abs(conditional$mean - unname(.bef_site_draws(fit, "theta_mean"))) / stan_sd),
            1e-6)
  expect_lt(relative(sqrt(conditional$var), stan_sd), 1e-6)
  expect_lt(max(abs(.bef_log_lik(engine) - unname(unclass(.bef_log_lik_draws(fit))))), 1e-6)

  # The reference draws of version 0.3.0 keep two sites. log_g is formed
  # again from alpha and the basis, which the reference identifies.
  ref <- readRDS(test_path("_fixtures", "ref-v030.rds"))
  # The hash is of the bits of the basis, which can differ in the last place
  # on another platform than the one of the reference.
  same_platform <- identical(R.version$platform, ref$conditions$platform)
  for (case in ref$cases) {
    grid <- seq(case$grid_range[1L], case$grid_range[2L], length.out = case$L)
    B <- splines::ns(grid, df = case$M, intercept = FALSE)
    if (same_platform) {
      expect_identical(digest::digest(unname(matrix(as.numeric(B), case$L)), algo = "sha256"),
                       case$basis_sha256)
    }
    alpha <- case$draws[, sprintf("alpha[%d]", seq_len(case$M)), drop = FALSE]
    eta <- alpha %*% t(matrix(as.numeric(B), case$L))
    top <- apply(eta, 1L, max)
    log_g <- eta - (top + log(rowSums(exp(eta - top))))
    kept <- sprintf("log_g[%d]", case$grid_points)
    expect_lt(max(abs(log_g[, case$grid_points] - case$draws[, kept])), 1e-5)

    engine <- .bef_engine_from_parts(log_g, grid, case$args$theta_hat, case$args$sigma)
    conditional <- .bef_site_conditional(engine, sites = case$sites)
    pick <- function(name) unname(case$draws[, sprintf("%s[%d]", name, case$sites)])
    expect_lt(max(abs(conditional$mean - pick("theta_mean")) / pick("theta_sd")), 1e-5)
    expect_lt(relative(sqrt(conditional$var), pick("theta_sd")), 1e-5)
    expect_lt(max(abs(.bef_log_lik(engine, sites = case$sites) - pick("log_lik"))), 1e-5)
  }
})
