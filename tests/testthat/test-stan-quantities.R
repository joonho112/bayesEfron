# These checks run the Stan program as it is now, and not a stored fit.
# Fixed parameters separate the arithmetic from the behavior of short chains.

stan_quantities_data <- function(theta_hat = c(-0.45, -0.1, 0, 0.25, 0.55),
                                 sigma = c(0.12, 0.18, 0.15, 0.22, 0.2),
                                 L = 51L, M = 3L, tau_c = NULL) {
  grid <- make_efron_grid(theta_hat, sigma, L = L, M = M, grid_method = "paper_realdata")
  c(
    list(K = length(theta_hat), theta_hat = theta_hat, sigma = sigma,
         L = grid$L, grid = grid$grid, M = grid$M, B = unname(grid$B)),
    # The program of version 0.3.0 has no tau_c.
    if (!is.null(tau_c)) list(tau_c = tau_c),
    list(store_grid_quantities = 1L)
  )
}

stan_quantities_log_sum_exp <- function(x) {
  largest <- max(x)
  largest + log(sum(exp(x - largest)))
}

# `tolerance` is for the summaries and `exact` for the quantities that follow
# from the coefficients without any averaging. Both have to be looser for
# draws read with CmdStan's default number of digits. With `sites = FALSE`
# the draws hold nothing for the sites, as the program efron_mix.stan writes
# them.
expect_stan_quantities <- function(draws, data, tolerance = 1e-7,
                                   exact = 1e-9, sites = TRUE) {
  draws <- posterior::as_draws_matrix(draws)
  draws <- matrix(as.numeric(draws), nrow = nrow(draws),
                  dimnames = list(NULL, colnames(draws)))
  indexed <- function(row, name, size) {
    as.numeric(draws[row, sprintf("%s[%d]", name, seq_len(size))])
  }
  scalar <- function(row, name) as.numeric(draws[row, name])
  expect_true(all(is.finite(draws)))
  tau_c <- if (is.null(data$tau_c)) 0 else data$tau_c
  # Within a component the effect is normal about a point between the grid
  # point and the estimate, with variance `shrink * sigma^2`.
  shrink <- tau_c^2 / (tau_c^2 + data$sigma^2)

  for (row in seq_len(nrow(draws))) {
    alpha <- indexed(row, "alpha", data$M)
    log_w <- drop(data$B %*% alpha)
    log_g <- log_w - stan_quantities_log_sum_exp(log_w)
    expect_lt(max(abs(indexed(row, "log_g", data$L) - log_g)), exact)
    g <- exp(log_g)
    mean_g <- sum(g * data$grid)
    var_g <- sum(g * (data$grid - mean_g)^2) + tau_c^2
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
        data$theta_hat[[site]], data$grid,
        sqrt(data$sigma[[site]]^2 + tau_c^2), log = TRUE
      )
      log_lik[[site]] <- stan_quantities_log_sum_exp(log_post)
      weights <- exp(log_post - max(log_post))
      weights <- weights / sum(weights)
      centers <- data$grid + shrink[[site]] * (data$theta_hat[[site]] - data$grid)
      theta_mean[[site]] <- sum(weights * centers)
      theta_sd[[site]] <- sqrt(
        shrink[[site]] * data$sigma[[site]]^2 +
          sum(weights * (centers - theta_mean[[site]])^2)
      )
      theta_map[[site]] <- centers[[which.max(log_post)]]
    }
    expect_equal(scalar(row, "log_marginal_likelihood"), sum(log_lik),
                 tolerance = exact)
    if (!isTRUE(sites)) {
      next
    }
    expect_lt(max(abs(indexed(row, "theta_mean", data$K) - theta_mean)),
              tolerance)
    expect_equal(indexed(row, "theta_sd", data$K), theta_sd,
                 tolerance = tolerance)
    expect_lt(max(abs(indexed(row, "theta_map", data$K) - theta_map)), tolerance)
    expect_equal(indexed(row, "log_lik", data$K), log_lik, tolerance = exact)
    expect_equal(scalar(row, "effective_params"), sum(theta_sd^2 / data$sigma^2),
                 tolerance = tolerance)

    theta_rep <- indexed(row, "theta_rep", data$K)
    if (tau_c == 0) {
      distances <- vapply(theta_rep, function(x) min(abs(x - data$grid)), numeric(1))
      expect_lt(max(distances), tolerance)
    }
  }
  invisible(TRUE)
}

stan_quantities_fixed_fit <- function(model, data, alpha) {
  # Without warmup there is nothing to adapt, and CmdStan 2.34 refuses a
  # run without warmup unless adaptation is turned off.
  model$sample(data = data, chains = 1L, parallel_chains = 1L,
               iter_warmup = 0L, iter_sampling = 1L, fixed_param = TRUE,
               adapt_engaged = FALSE,
               init = list(list(alpha = alpha, lambda = 0.5)),
               seed = 1985L, refresh = 0L, sig_figs = 18L)
}

stan_quantities_model <- function(program) {
  getFromNamespace(".bef_model", "bayesEfron")(program = program)
}

# The draws of efron_mix.stan with the site-level quantities that the package
# computes in R added to them.
stan_quantities_completed <- function(fit, data) {
  getFromNamespace(".bef_complete_mix_draws", "bayesEfron")(
    fit$draws(format = "draws_array"), stan_data = data, site_draws = TRUE,
    seed = 77L
  )
}

test_that("the generated quantities agree with independent arithmetic, also far from zero", {
  skip_if_no_cmdstan()
  model <- stan_quantities_model("efron_re")
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

test_that("the program of this version writes g only, and R adds the same site quantities", {
  skip_if_no_cmdstan()
  model <- stan_quantities_model("efron_mix")
  alpha <- c(0.2, -0.3, 0.5)

  for (tau_c in c(0, 0.04)) {
    data <- stan_quantities_data(tau_c = tau_c)
    fit <- stan_quantities_fixed_fit(model, data, alpha)
    written <- dimnames(fit$draws())[[3L]]
    # A run with fixed parameters has no lp__.
    expect_setequal(
      setdiff(unique(sub("\\[.*$", "", written)), "lp__"),
      c("alpha", "lambda", "log_g", "mean_g", "var_g", "sd_g",
        "log_marginal_likelihood", "log_w", "g")
    )
    expect_stan_quantities(fit$draws(), data, sites = FALSE)

    completed <- stan_quantities_completed(fit, data)
    expect_s3_class(completed, "draws_array")
    expect_stan_quantities(completed, data)

    offset <- 1e8
    translated_data <- data
    translated_data$theta_hat <- data$theta_hat + offset
    translated_data$grid <- data$grid + offset
    translated <- stan_quantities_fixed_fit(model, translated_data, alpha)
    expect_stan_quantities(stan_quantities_completed(translated, translated_data),
                           translated_data)
  }
})

test_that("with point masses the two programs give the same quantities for the same coefficients", {
  skip_if_no_cmdstan()
  alpha <- c(0.2, -0.3, 0.5)
  data <- stan_quantities_data()
  old <- stan_quantities_fixed_fit(stan_quantities_model("efron_re"), data, alpha)
  mix_data <- stan_quantities_data(tau_c = 0)
  new <- stan_quantities_completed(
    stan_quantities_fixed_fit(stan_quantities_model("efron_mix"), mix_data, alpha),
    mix_data
  )
  old <- posterior::as_draws_matrix(old$draws())
  new <- posterior::as_draws_matrix(new)

  # The variables and their order are those of the program of version 0.3.0.
  expect_identical(colnames(new), colnames(old))
  same <- setdiff(colnames(old), c("lp__", sprintf("theta_rep[%d]", seq_len(data$K))))
  expect_equal(as.numeric(new[, same]), as.numeric(old[, same]), tolerance = 1e-12)
})

test_that("a sampled fit agrees with the same arithmetic and records its settings", {
  skip_if_no_cmdstan()
  data <- stan_quantities_data()
  fit <- suppressWarnings(bayes_efron_fit(
    theta_hat = data$theta_hat, sigma = data$sigma,
    L = data$L, M = data$M, chains = 1L, parallel_chains = 1L,
    iter_warmup = 150L, iter_sampling = 50L, seed = 4701L, site_draws = TRUE
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

test_that("with one normal component the program of this version is the normal-normal model", {
  skip_if_no_cmdstan()
  # One grid point, at 0.3, with a normal component of standard deviation
  # 0.25: every effect is N(0.3, 0.25^2), whatever the coefficient is. The
  # grid functions of the package do not make such a grid; the program takes
  # it. The chains are sampled, and not held fixed.
  theta_hat <- c(-0.45, -0.1, 0, 0.25, 0.55, 0.9)
  sigma <- c(0.12, 0.18, 0.15, 0.22, 0.2, 0.3)
  mu <- 0.3
  tau <- 0.25
  data <- list(K = length(theta_hat), theta_hat = theta_hat, sigma = sigma, L = 1L,
               grid = array(mu, dim = 1L), M = 1L, B = matrix(0, 1L, 1L), tau_c = tau,
               store_grid_quantities = 0L)
  fit <- suppressWarnings(stan_quantities_model("efron_mix")$sample(
    data = data, chains = 2L, parallel_chains = 2L, iter_warmup = 200L, iter_sampling = 100L,
    seed = 4046L, refresh = 0L, init = 0.5, sig_figs = 18L, show_messages = FALSE,
    show_exceptions = FALSE, diagnostics = NULL
  ))
  completed <- getFromNamespace(".bef_complete_mix_draws", "bayesEfron")(
    fit$draws(format = "draws_array"), stan_data = data, site_draws = TRUE, seed = 9L
  )
  draws <- posterior::as_draws_matrix(completed)
  draws <- matrix(as.numeric(draws), nrow = nrow(draws), dimnames = list(NULL, colnames(draws)))
  n <- nrow(draws)
  expect_identical(n, 200L)
  # The coefficient moves; nothing that describes g or the sites does.
  expect_gt(stats::sd(draws[, "alpha[1]"]), 0.1)
  expect_equal(draws[, "log_g[1]"], rep(0, n), ignore_attr = TRUE)
  expect_equal(draws[, "mean_g"], rep(mu, n), ignore_attr = TRUE, tolerance = 1e-12)
  expect_equal(draws[, "sd_g"], rep(tau, n), ignore_attr = TRUE, tolerance = 1e-12)
  expect_equal(draws[, "var_g"], rep(tau^2, n), ignore_attr = TRUE, tolerance = 1e-12)
  marginal <- sum(stats::dnorm(theta_hat, mu, sqrt(sigma^2 + tau^2), log = TRUE))
  expect_equal(draws[, "log_marginal_likelihood"], rep(marginal, n), ignore_attr = TRUE,
               tolerance = 1e-12)

  # The normal-normal posterior of each site, in closed form.
  shrink <- tau^2 / (tau^2 + sigma^2)
  post_mean <- mu + shrink * (theta_hat - mu)
  post_sd <- sqrt(shrink) * sigma
  for (i in seq_along(theta_hat)) {
    expect_equal(draws[, sprintf("theta_mean[%d]", i)], rep(post_mean[[i]], n),
                 ignore_attr = TRUE, tolerance = 1e-12)
    expect_equal(draws[, sprintf("theta_sd[%d]", i)], rep(post_sd[[i]], n),
                 ignore_attr = TRUE, tolerance = 1e-12)
    expect_equal(draws[, sprintf("theta_map[%d]", i)], rep(post_mean[[i]], n),
                 ignore_attr = TRUE, tolerance = 1e-12)
    expect_equal(draws[, sprintf("log_lik[%d]", i)],
                 rep(stats::dnorm(theta_hat[[i]], mu, sqrt(sigma[[i]]^2 + tau^2), log = TRUE), n),
                 ignore_attr = TRUE, tolerance = 1e-12)
    # The draws of the effect are draws from that normal distribution.
    rep_i <- draws[, sprintf("theta_rep[%d]", i)]
    expect_lt(abs(mean(rep_i) - post_mean[[i]]), 4 * post_sd[[i]] / sqrt(n))
    expect_gt(stats::ks.test(rep_i, "pnorm", post_mean[[i]], post_sd[[i]])$p.value, 1e-4)
  }
  # The weight of the data in each posterior mean, summed over the sites.
  expect_equal(draws[, "effective_params"], rep(sum(shrink), n), ignore_attr = TRUE,
               tolerance = 1e-12)
})
