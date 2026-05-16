tier1_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

tier1_active_target <- function(target_id) {
  .bef_target(target_id, statuses = "active")
}

tier1_target_any_status <- function(target_id) {
  .bef_target(target_id, statuses = c("active", "deferred"))
}

tier1_tolerance <- function(target) {
  as.numeric(target$tolerance_value)
}

tier1_stan_source <- local({
  source <- NULL

  function() {
    if (is.null(source)) {
      source <<- paste(
        readLines(.bef_source_file("inst", "stan", "efron_re.stan"), warn = FALSE),
        collapse = "\n"
      )
    }
    source
  }
})

tier1_expect_stan_source <- function(pattern) {
  expect_match(tier1_stan_source(), pattern, perl = TRUE)
}

tier1_softmax <- function(log_w) {
  shifted <- log_w - max(log_w)
  exp_shifted <- exp(shifted)
  exp_shifted / sum(exp_shifted)
}

tier1_discrete_posterior <- function(theta_hat, sigma, grid, g) {
  K <- length(theta_hat)
  theta_map <- theta_mean <- theta_sd <- numeric(K)
  for (site in seq_len(K)) {
    log_post <- log(g) + stats::dnorm(
      theta_hat[site],
      mean = grid,
      sd = sigma[site],
      log = TRUE
    )
    weights <- tier1_softmax(log_post)
    theta_map[site] <- grid[which.max(log_post)]
    theta_mean[site] <- sum(weights * grid)
    second_moment <- sum(weights * grid^2)
    theta_sd[site] <- sqrt(max(0, second_moment - theta_mean[site]^2))
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

tier1_gh_normal <- function(n, mu, tau) {
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

tier1_normal_normal_closed_form <- function(theta_hat, sigma, mu, tau) {
  variance <- 1 / (1 / tau^2 + 1 / sigma^2)
  list(
    mean = variance * (mu / tau^2 + theta_hat / sigma^2),
    sd = sqrt(variance)
  )
}

test_that("Tier 1 grid basis agrees with the splines::ns reference", {
  target <- tier1_active_target("tier1_make_efron_grid_ns_agreement")

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
    tolerance = tier1_tolerance(target)
  )
  expect_equal(attr(grid$B, "knots"), attr(reference, "knots"))
  expect_equal(attr(grid$B, "Boundary.knots"), attr(reference, "Boundary.knots"))
  expect_identical(attr(grid$B, "intercept"), attr(reference, "intercept"))
})

test_that("Tier 1 discrete prior g agrees with prop-table reference", {
  target <- tier1_active_target("tier1_g_prop_table_discrete_prior")
  tier1_expect_stan_source("vector\\[L\\] log_w = B \\* alpha")
  tier1_expect_stan_source("simplex\\[L\\] g = softmax\\(log_w\\)")

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
  g_stable <- tier1_softmax(log_w)

  expect_equal(g_stable, g_reference, tolerance = tier1_tolerance(target))
  expect_equal(sum(g_stable), 1, tolerance = tier1_tolerance(target))
  expect_true(all(g_stable > 0))

  perturbed_alpha <- alpha
  perturbed_alpha[[2L]] <- perturbed_alpha[[2L]] + 0.25
  g_perturbed <- tier1_softmax(drop(as.matrix(grid$B) %*% perturbed_alpha))
  expect_gt(max(abs(g_perturbed - g_reference)), tier1_tolerance(target))

  g_permuted <- tier1_softmax(drop(as.matrix(grid$B[nrow(grid$B):1L, ]) %*% alpha))
  expect_gt(max(abs(g_permuted - g_reference)), tier1_tolerance(target))
})

test_that("Tier 1 K2 normal-normal posterior agrees with closed form", {
  target <- tier1_active_target("tier1_normal_normal_K2_posterior")
  tier1_expect_stan_source("log_post\\[j\\] = log_g\\[j\\]")
  tier1_expect_stan_source("normal_lpdf\\(theta_hat\\[i\\] \\| grid\\[j\\], sigma\\[i\\]\\)")
  tier1_expect_stan_source("theta_mean\\[i\\] = dot_product\\(w, grid\\)")
  tier1_expect_stan_source("theta_sd\\[i\\] = sqrt\\(second_moment - square\\(theta_mean\\[i\\]\\)\\)")

  mu <- 0.15
  tau <- 0.7
  theta_hat <- c(-0.35, 0.8)
  sigma <- c(0.2, 0.45)
  gh <- tier1_gh_normal(n = 151L, mu = mu, tau = tau)

  posterior <- tier1_discrete_posterior(
    theta_hat = theta_hat,
    sigma = sigma,
    grid = gh$nodes,
    g = gh$weights
  )
  closed_form <- tier1_normal_normal_closed_form(
    theta_hat = theta_hat,
    sigma = sigma,
    mu = mu,
    tau = tau
  )

  expect_equal(
    posterior$theta_mean,
    closed_form$mean,
    tolerance = tier1_tolerance(target)
  )
  expect_equal(
    posterior$theta_sd,
    closed_form$sd,
    tolerance = tier1_tolerance(target)
  )

  swapped <- tier1_discrete_posterior(
    theta_hat = theta_hat,
    sigma = rev(sigma),
    grid = gh$nodes,
    g = gh$weights
  )
  expect_gt(max(abs(swapped$theta_mean - closed_form$mean)), tier1_tolerance(target))
  expect_gt(max(abs(swapped$theta_sd - closed_form$sd)), tier1_tolerance(target))
})

test_that("Tier 1 prior moments agree with Gauss-Hermite references", {
  mean_target <- tier1_active_target("tier1_mean_g_gauss_hermite")
  var_target <- tier1_active_target("tier1_var_g_gauss_hermite")
  tier1_expect_stan_source("real mean_g = dot_product\\(g, grid\\)")
  tier1_expect_stan_source("real var_g = dot_product\\(g, square\\(grid - mean_g\\)\\)")

  mu <- -0.2
  tau <- 0.9
  gh <- tier1_gh_normal(n = 151L, mu = mu, tau = tau)
  posterior <- tier1_discrete_posterior(
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
    tolerance = tier1_tolerance(mean_target)
  )
  expect_equal(mean_reference, mu, tolerance = tier1_tolerance(mean_target))

  expect_equal(
    posterior$var_g,
    var_reference,
    tolerance = tier1_tolerance(var_target)
  )
  expect_equal(var_reference, tau^2, tolerance = tier1_tolerance(var_target))

  shifted_mean <- sum(gh$weights * (gh$nodes + 0.01))
  expect_gt(abs(shifted_mean - mean_reference), tier1_tolerance(mean_target))

  perturbed_weights <- gh$weights
  center <- which.min(abs(gh$nodes - mu))
  perturbed_weights[[center]] <- perturbed_weights[[center]] + 0.01
  perturbed_weights <- perturbed_weights / sum(perturbed_weights)
  perturbed_mean <- sum(perturbed_weights * gh$nodes)
  perturbed_var <- sum(perturbed_weights * (gh$nodes - perturbed_mean)^2)
  expect_gt(abs(perturbed_var - var_reference), tier1_tolerance(var_target))
})

test_that("Tier 1 theta_map agrees with on-grid optim reference", {
  target <- tier1_active_target("tier1_theta_map_optim_agreement")
  tier1_expect_stan_source("int max_idx = 1")
  tier1_expect_stan_source("if \\(log_post\\[j\\] > log_post\\[max_idx\\]\\)")
  tier1_expect_stan_source("theta_map\\[i\\] = grid\\[max_idx\\]")

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

  posterior <- tier1_discrete_posterior(
    theta_hat = theta_hat,
    sigma = sigma,
    grid = grid,
    g = g
  )

  expect_equal(optim_modes, desired_modes, tolerance = tier1_tolerance(target))
  expect_equal(posterior$theta_map, optim_modes, tolerance = tier1_tolerance(target))

  off_grid <- grid[!grid %in% desired_modes]
  off_grid_g <- stats::dnorm(off_grid, mean = mu, sd = tau)
  off_grid_g <- off_grid_g / sum(off_grid_g)
  off_grid_posterior <- tier1_discrete_posterior(
    theta_hat = theta_hat,
    sigma = sigma,
    grid = off_grid,
    g = off_grid_g
  )
  expect_gt(
    max(abs(off_grid_posterior$theta_map - optim_modes)),
    tier1_tolerance(target)
  )

  tie_grid <- c(-0.1, 0.1)
  tie_g <- c(0.5, 0.5)
  tie_posterior <- tier1_discrete_posterior(
    theta_hat = 0,
    sigma = 0.3,
    grid = tie_grid,
    g = tie_g
  )
  expect_identical(tie_posterior$theta_map, tie_grid[[1L]])
})

test_that("Tier 1 posterior quantiles agree with base R references", {
  target <- tier1_active_target("tier1_posterior_quantile_base_r")
  probs <- tier1_ns(".bef_interval_probs")(0.8)

  draws <- c(-1.2, -0.4, 0, 0.35, 0.9, 1.4, 2.1)
  expect_equal(
    posterior::quantile2(draws, probs = probs, names = FALSE),
    stats::quantile(draws, probs = probs, names = FALSE, type = 7),
    tolerance = tier1_tolerance(target)
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
    tier1_ns(".bef_matrix_quantiles")(theta_rep, probs = probs),
    base_quantiles,
    tolerance = tier1_tolerance(target)
  )
})
