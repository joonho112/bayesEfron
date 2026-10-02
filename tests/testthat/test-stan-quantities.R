# These checks run the Stan program as it is now, and not a stored fit.
# Fixed parameters separate the arithmetic from the behavior of short chains.

stan_quantities_data <- function(theta_hat = c(-0.45, -0.1, 0, 0.25, 0.55),
                                 sigma = c(0.12, 0.18, 0.15, 0.22, 0.2),
                                 L = 51L, M = 3L) {
  grid <- make_efron_grid(theta_hat, sigma, L = L, M = M)
  list(K = length(theta_hat), theta_hat = theta_hat, sigma = sigma,
       L = grid$L, grid = grid$grid, M = grid$M, B = unname(grid$B),
       store_grid_quantities = 1L)
}

stan_quantities_log_sum_exp <- function(x) {
  largest <- max(x)
  largest + log(sum(exp(x - largest)))
}

# `tolerance` is for the summaries and `exact` for the quantities that follow
# from the coefficients without any averaging. Both have to be looser for
# draws read with CmdStan's default number of digits.
expect_stan_quantities <- function(draws, data, tolerance = 1e-7,
                                   exact = 1e-9) {
  draws <- posterior::as_draws_matrix(draws)
  draws <- matrix(as.numeric(draws), nrow = nrow(draws),
                  dimnames = list(NULL, colnames(draws)))
  indexed <- function(row, name, size) {
    as.numeric(draws[row, sprintf("%s[%d]", name, seq_len(size))])
  }
  scalar <- function(row, name) as.numeric(draws[row, name])
  expect_true(all(is.finite(draws)))

  for (row in seq_len(nrow(draws))) {
    alpha <- indexed(row, "alpha", data$M)
    log_w <- drop(data$B %*% alpha)
    log_g <- log_w - stan_quantities_log_sum_exp(log_w)
    expect_lt(max(abs(indexed(row, "log_g", data$L) - log_g)), exact)
    g <- exp(log_g)
    mean_g <- sum(g * data$grid)
    var_g <- sum(g * (data$grid - mean_g)^2)
    # Use absolute errors in effect units, rather than relative tolerances
    # that become permissive when estimates have a large common offset.
    expect_lt(abs(draws[row, "mean_g"] - mean_g), tolerance)
    expect_equal(scalar(row, "var_g"), var_g, tolerance = tolerance)
    expect_equal(scalar(row, "sd_g"), sqrt(var_g), tolerance = tolerance)

    if (data$store_grid_quantities == 1L) {
      expect_equal(indexed(row, "log_w", data$L), log_w, tolerance = exact)
      expect_equal(indexed(row, "g", data$L), g, tolerance = exact)
    }

    theta_mean <- theta_sd <- theta_map <- log_lik <- numeric(data$K)
    for (site in seq_len(data$K)) {
      log_post <- log_g + stats::dnorm(
        data$theta_hat[[site]], data$grid, data$sigma[[site]], log = TRUE
      )
      log_lik[[site]] <- stan_quantities_log_sum_exp(log_post)
      weights <- exp(log_post - max(log_post))
      weights <- weights / sum(weights)
      theta_mean[[site]] <- sum(weights * data$grid)
      theta_sd[[site]] <- sqrt(sum(weights * (data$grid - theta_mean[[site]])^2))
      theta_map[[site]] <- data$grid[[which.max(log_post)]]
    }
    expect_lt(max(abs(indexed(row, "theta_mean", data$K) - theta_mean)),
              tolerance)
    expect_equal(indexed(row, "theta_sd", data$K), theta_sd,
                 tolerance = tolerance)
    expect_lt(max(abs(indexed(row, "theta_map", data$K) - theta_map)), tolerance)
    expect_equal(indexed(row, "log_lik", data$K), log_lik, tolerance = exact)
    expect_equal(scalar(row, "log_marginal_likelihood"), sum(log_lik),
                 tolerance = exact)
    expect_equal(scalar(row, "effective_params"), sum(theta_sd^2 / data$sigma^2),
                 tolerance = tolerance)

    theta_rep <- indexed(row, "theta_rep", data$K)
    distances <- vapply(theta_rep, function(x) min(abs(x - data$grid)), numeric(1))
    expect_lt(max(distances), tolerance)
  }
  invisible(TRUE)
}

stan_quantities_fixed_fit <- function(model, data, alpha) {
  model$sample(data = data, chains = 1L, parallel_chains = 1L,
               iter_warmup = 0L, iter_sampling = 1L, fixed_param = TRUE,
               init = list(list(alpha = alpha, lambda = 0.5)),
               seed = 1985L, refresh = 0L, sig_figs = 18L)
}

test_that("the generated quantities agree with independent arithmetic, also far from zero", {
  skip_if_no_cmdstan()
  model <- bayes_efron_compile()
  data <- stan_quantities_data()
  alpha <- c(0.2, -0.3, 0.5)
  ordinary <- stan_quantities_fixed_fit(model, data, alpha)
  expect_stan_quantities(ordinary$draws(), data)

  offset <- 1e8
  translated_data <- data
  translated_data$theta_hat <- data$theta_hat + offset
  translated_data$grid <- data$grid + offset
  translated <- stan_quantities_fixed_fit(model, translated_data, alpha)
  expect_stan_quantities(translated$draws(), translated_data)

  first <- as.matrix(ordinary$draws(format = "draws_matrix"))
  second <- as.matrix(translated$draws(format = "draws_matrix"))
  means <- c("mean_g", sprintf("theta_mean[%d]", seq_len(data$K)))
  scales <- c("var_g", "sd_g", sprintf("theta_sd[%d]", seq_len(data$K)),
              "effective_params")
  expect_lt(max(abs(second[, means] - offset - first[, means])), 1e-7)
  expect_equal(as.numeric(second[, scales]), as.numeric(first[, scales]),
               tolerance = 1e-7)
})

test_that("a sampled fit agrees with the same arithmetic and records its settings", {
  skip_if_no_cmdstan()
  data <- stan_quantities_data()
  fit <- suppressWarnings(bayes_efron_fit(
    theta_hat = data$theta_hat, sigma = data$sigma,
    L = data$L, M = data$M, chains = 1L, parallel_chains = 1L,
    iter_warmup = 150L, iter_sampling = 50L, seed = 4701L
  ))
  # The draws are read with eight significant digits.
  expect_stan_quantities(fit$draws, fit$metadata$data_list,
                         tolerance = 1e-5, exact = 1e-5)
  settings <- attr(fit$metadata, "sampler_settings", exact = TRUE)
  expect_identical(settings$iter_sampling, 50L)
  expect_identical(settings$seed, 4701L)
  expect_identical(settings$init, 0.5)
  expect_identical(attr(fit$metadata, "package_version", exact = TRUE),
                   as.character(getNamespaceVersion("bayesEfron")))
  expect_identical(attr(fit$metadata, "cmdstanr_version", exact = TRUE),
                   as.character(getNamespaceVersion("cmdstanr")))
})
