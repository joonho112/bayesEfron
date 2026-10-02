# Closed-form checks. The formulas of the model are written out in R below
# (softmax, discrete_posterior) and compared with independent references: a
# conjugate normal-normal posterior, Gauss-Hermite quadrature, optim. The last
# test then compares what the Stan program wrote for the stored fit with the
# same R formulas, so that the chain runs from the closed forms to the output
# of the package. Section numbers refer to Lee and Sui (2025),
# doi:10.3390/math13162639.

bef_internal <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

softmax <- function(log_w) {
  shifted <- log_w - max(log_w)
  exp_shifted <- exp(shifted)
  exp_shifted / sum(exp_shifted)
}

discrete_posterior <- function(theta_hat, sigma, grid, g) {
  K <- length(theta_hat)
  theta_map <- theta_mean <- theta_sd <- numeric(K)
  for (site in seq_len(K)) {
    log_post <- log(g) + stats::dnorm(
      theta_hat[site],
      mean = grid,
      sd = sigma[site],
      log = TRUE
    )
    weights <- softmax(log_post)
    theta_map[site] <- grid[which.max(log_post)]
    theta_mean[site] <- sum(weights * grid)
    theta_sd[site] <- sqrt(sum(weights * (grid - theta_mean[site])^2))
  }
  mean_g <- sum(g * grid)
  list(
    mean_g = mean_g,
    var_g = sum(g * (grid - mean_g)^2),
    theta_map = theta_map,
    theta_mean = theta_mean,
    theta_sd = theta_sd
  )
}

gauss_hermite_normal <- function(n, mu, tau) {
  J <- matrix(0, nrow = n, ncol = n)
  off_diagonal <- sqrt(seq_len(n - 1L))
  J[cbind(seq_len(n - 1L), 2:n)] <- off_diagonal
  J[cbind(2:n, seq_len(n - 1L))] <- off_diagonal

  eig <- eigen(J, symmetric = TRUE)
  order_nodes <- order(eig$values)
  weights <- eig$vectors[1L, order_nodes]^2
  weights <- weights / sum(weights)
  list(
    nodes = as.numeric(mu + tau * eig$values[order_nodes]),
    weights = as.numeric(weights)
  )
}

normal_normal_posterior <- function(theta_hat, sigma, mu, tau) {
  variance <- 1 / (1 / tau^2 + 1 / sigma^2)
  list(
    mean = variance * (mu / tau^2 + theta_hat / sigma^2),
    sd = sqrt(variance)
  )
}

# Section 2.4: the basis matrix Q is a natural cubic spline basis on the grid.
test_that("the grid basis agrees with splines::ns", {
  tol <- 1e-6

  theta_hat <- c(-0.85, -0.35, -0.05, 0.2, 0.65, 0.9)
  sigma <- c(0.14, 0.19, 0.16, 0.21, 0.18, 0.23)
  grid <- bayesEfron::make_efron_grid(
    theta_hat = theta_hat,
    sigma = sigma,
    L = 61L,
    expansion = 0.4,
    M = 5L,
    grid_method = "paper_realdata"
  )

  reference <- splines::ns(grid$grid, df = grid$M, intercept = FALSE)

  expect_equal(
    unname(as.matrix(grid$B)),
    unname(as.matrix(reference)),
    tolerance = tol
  )
  expect_equal(attr(grid$B, "knots"), attr(reference, "knots"))
  expect_equal(attr(grid$B, "Boundary.knots"), attr(reference, "Boundary.knots"))
  expect_identical(attr(grid$B, "intercept"), attr(reference, "intercept"))
})

# Section 2.4, Equation (4): g_j = exp(Q_j' alpha - phi(alpha)).
test_that("the discrete prior g agrees with prop.table", {
  tol <- 1e-6

  theta_hat <- c(-0.9, -0.2, 0.4, 1.1)
  sigma <- c(0.18, 0.27, 0.22, 0.31)
  alpha <- c(-0.35, 0.2, 0.65, -0.1)
  grid <- bayesEfron::make_efron_grid(
    theta_hat = theta_hat,
    sigma = sigma,
    L = 61L,
    M = length(alpha),
    expansion = 0.35,
    grid_method = "paper_realdata"
  )

  log_w <- drop(as.matrix(grid$B) %*% alpha)
  g_reference <- as.numeric(prop.table(exp(log_w)))
  g_stable <- softmax(log_w)

  expect_equal(g_stable, g_reference, tolerance = tol)
  expect_equal(sum(g_stable), 1, tolerance = tol)
  expect_true(all(g_stable > 0))

  perturbed_alpha <- alpha
  perturbed_alpha[[2L]] <- perturbed_alpha[[2L]] + 0.25
  g_perturbed <- softmax(drop(as.matrix(grid$B) %*% perturbed_alpha))
  expect_gt(max(abs(g_perturbed - g_reference)), tol)

  g_permuted <- softmax(drop(as.matrix(grid$B[nrow(grid$B):1L, ]) %*% alpha))
  expect_gt(max(abs(g_permuted - g_reference)), tol)
})

# Section 3.3: posterior weights, mean and variance of one site given g.
test_that("the two-site normal-normal posterior agrees with the closed form", {
  tol <- 1e-6

  mu <- 0.15
  tau <- 0.7
  theta_hat <- c(-0.35, 0.8)
  sigma <- c(0.2, 0.45)
  gh <- gauss_hermite_normal(n = 151L, mu = mu, tau = tau)

  posterior <- discrete_posterior(
    theta_hat = theta_hat,
    sigma = sigma,
    grid = gh$nodes,
    g = gh$weights
  )
  closed_form <- normal_normal_posterior(
    theta_hat = theta_hat,
    sigma = sigma,
    mu = mu,
    tau = tau
  )

  expect_equal(
    posterior$theta_mean,
    closed_form$mean,
    tolerance = tol
  )
  expect_equal(
    posterior$theta_sd,
    closed_form$sd,
    tolerance = tol
  )

  swapped <- discrete_posterior(
    theta_hat = theta_hat,
    sigma = rev(sigma),
    grid = gh$nodes,
    g = gh$weights
  )
  expect_gt(max(abs(swapped$theta_mean - closed_form$mean)), tol)
  expect_gt(max(abs(swapped$theta_sd - closed_form$sd)), tol)
})

# Section 4.3: functionals of g, here its mean and variance.
test_that("the mean and variance of g agree with Gauss-Hermite quadrature", {
  tol <- 1e-6

  mu <- -0.2
  tau <- 0.9
  gh <- gauss_hermite_normal(n = 151L, mu = mu, tau = tau)
  posterior <- discrete_posterior(
    theta_hat = c(-0.15, 0.3),
    sigma = c(0.4, 0.65),
    grid = gh$nodes,
    g = gh$weights
  )

  mean_reference <- sum(gh$weights * gh$nodes)
  var_reference <- sum(gh$weights * (gh$nodes - mean_reference)^2)

  expect_equal(
    posterior$mean_g,
    mean_reference,
    tolerance = tol
  )
  expect_equal(mean_reference, mu, tolerance = tol)

  expect_equal(
    posterior$var_g,
    var_reference,
    tolerance = tol
  )
  expect_equal(var_reference, tau^2, tolerance = tol)

  shifted_mean <- sum(gh$weights * (gh$nodes + 0.01))
  expect_gt(abs(shifted_mean - mean_reference), tol)

  perturbed_weights <- gh$weights
  center <- which.min(abs(gh$nodes - mu))
  perturbed_weights[[center]] <- perturbed_weights[[center]] + 0.01
  perturbed_weights <- perturbed_weights / sum(perturbed_weights)
  perturbed_mean <- sum(perturbed_weights * gh$nodes)
  perturbed_var <- sum(perturbed_weights * (gh$nodes - perturbed_mean)^2)
  expect_gt(abs(perturbed_var - var_reference), tol)
})

# Section 4.3: the posterior mode of a site effect on the grid.
test_that("the posterior mode agrees with optim on the grid", {
  tol <- 1e-6

  mu <- 0.1
  tau <- 0.6
  sigma <- c(0.3, 0.45)
  desired_modes <- c(-0.25, 0.5)
  theta_hat <- desired_modes * (1 + sigma^2 / tau^2) - mu * sigma^2 / tau^2
  grid <- seq(-1, 1, by = 0.025)
  g <- stats::dnorm(grid, mean = mu, sd = tau)
  g <- g / sum(g)

  objective <- function(site) {
    function(theta) {
      -(
        stats::dnorm(theta, mean = mu, sd = tau, log = TRUE) +
          stats::dnorm(theta_hat[[site]], mean = theta, sd = sigma[[site]], log = TRUE)
      )
    }
  }
  optim_modes <- vapply(
    seq_along(theta_hat),
    function(site) {
      stats::optim(
        par = theta_hat[[site]],
        fn = objective(site),
        method = "BFGS",
        control = list(reltol = .Machine$double.eps)
      )$par
    },
    numeric(1)
  )

  posterior <- discrete_posterior(
    theta_hat = theta_hat,
    sigma = sigma,
    grid = grid,
    g = g
  )

  expect_equal(optim_modes, desired_modes, tolerance = tol)
  expect_equal(posterior$theta_map, optim_modes, tolerance = tol)

  off_grid <- grid[!grid %in% desired_modes]
  off_grid_g <- stats::dnorm(off_grid, mean = mu, sd = tau)
  off_grid_g <- off_grid_g / sum(off_grid_g)
  off_grid_posterior <- discrete_posterior(
    theta_hat = theta_hat,
    sigma = sigma,
    grid = off_grid,
    g = off_grid_g
  )
  expect_gt(
    max(abs(off_grid_posterior$theta_map - optim_modes)),
    tol
  )

  tie_grid <- c(-0.1, 0.1)
  tie_g <- c(0.5, 0.5)
  tie_posterior <- discrete_posterior(
    theta_hat = 0,
    sigma = 0.3,
    grid = tie_grid,
    g = tie_g
  )
  expect_identical(tie_posterior$theta_map, tie_grid[[1L]])
})

# Section 4.3: equal-tailed credible intervals from posterior draws.
test_that("posterior quantiles agree with stats::quantile", {
  tol <- 1e-12
  probs <- bef_internal(".bef_interval_probs")(0.8)

  draws <- c(-1.2, -0.4, 0, 0.35, 0.9, 1.4, 2.1)
  expect_equal(
    posterior::quantile2(draws, probs = probs, names = FALSE),
    stats::quantile(draws, probs = probs, names = FALSE, type = 7),
    tolerance = tol
  )

  theta_rep <- matrix(
    c(
      -1.1, -0.8, -0.2, 0.1, 0.4, 0.9,
      -0.7, -0.1, 0.2, 0.6, 1.0, 1.6,
      -1.4, -0.9, -0.3, 0.05, 0.8, 1.7
    ),
    nrow = 6L,
    ncol = 3L
  )
  base_quantiles <- vapply(
    seq_len(ncol(theta_rep)),
    function(col) {
      stats::quantile(theta_rep[, col], probs = probs, names = FALSE, type = 7)
    },
    numeric(length(probs))
  )

  expect_equal(
    bef_internal(".bef_matrix_quantiles")(theta_rep, probs = probs),
    base_quantiles,
    tolerance = tol
  )
})

# The generated quantities of the Stan program, recomputed in R from the
# stored draws of alpha and of log g with the functions checked above.
test_that("the quantities written by the Stan program agree with the R formulas", {
  fit <- bayesEfron::raudenbush_fit
  data <- fit$metadata$data_list
  draws <- posterior::as_draws_matrix(fit$draws)
  indexed <- function(name, n) {
    x <- draws[, sprintf("%s[%d]", name, seq_len(n)), drop = FALSE]
    matrix(as.numeric(x), nrow = nrow(x))
  }
  scalar <- function(name) as.numeric(draws[, name])

  alpha <- indexed("alpha", data$M)
  log_g <- indexed("log_g", data$L)
  theta_mean <- indexed("theta_mean", data$K)
  theta_sd <- indexed("theta_sd", data$K)
  theta_map <- indexed("theta_map", data$K)
  log_lik <- indexed("log_lik", data$K)
  mean_g <- scalar("mean_g")
  var_g <- scalar("var_g")
  sd_g <- scalar("sd_g")
  effective_params <- scalar("effective_params")
  log_marginal_likelihood <- scalar("log_marginal_likelihood")
  B <- unname(as.matrix(data$B))

  for (s in round(seq(1, nrow(alpha), length.out = 40))) {
    # log g = log softmax(B alpha); alpha is stored to fewer digits than the
    # sampler used, hence the wider tolerance.
    expect_equal(log(softmax(drop(B %*% alpha[s, ]))), log_g[s, ], tolerance = 1e-4)

    g <- exp(log_g[s, ])
    reference <- discrete_posterior(data$theta_hat, data$sigma, data$grid, g)
    expect_equal(reference$mean_g, mean_g[s], tolerance = 1e-6)
    expect_equal(reference$var_g, var_g[s], tolerance = 1e-6)
    expect_equal(sqrt(reference$var_g), sd_g[s], tolerance = 1e-6)
    expect_equal(reference$theta_mean, theta_mean[s, ], tolerance = 1e-6)
    expect_equal(reference$theta_sd, theta_sd[s, ], tolerance = 1e-6)
    expect_equal(reference$theta_map, theta_map[s, ], tolerance = 1e-6)

    site_log_lik <- vapply(seq_len(data$K), function(i) {
      log_post <- log_g[s, ] +
        stats::dnorm(data$theta_hat[i], mean = data$grid, sd = data$sigma[i], log = TRUE)
      max(log_post) + log(sum(exp(log_post - max(log_post))))
    }, numeric(1))
    expect_equal(site_log_lik, log_lik[s, ], tolerance = 1e-6)
    expect_equal(sum(site_log_lik), log_marginal_likelihood[s], tolerance = 1e-6)
    expect_equal(sum(reference$theta_sd^2 / data$sigma^2), effective_params[s],
                 tolerance = 1e-6)
  }

  # The draws of a site effect lie on the grid.
  theta_rep <- unique(as.numeric(indexed("theta_rep", data$K)))
  distance_to_grid <- vapply(theta_rep, function(x) min(abs(x - data$grid)), numeric(1))
  expect_lt(max(distance_to_grid), 1e-8)
})
