# Small fits for the tests of the plot functions, made without CmdStan: the
# draws of the spline coefficients are drawn here around coefficients that
# give a plausible g, and the package turns them into a fit as it does with
# the output of CmdStan.

plot_fit_ns <- function(name) getFromNamespace(name, "bayesEfron")

plot_fit_cache <- new.env(parent = emptyenv())

# A fit of `K` sites with `n_draws` draws in two chains. `labels` gives the
# sites names; `component = "normal"` puts a normal distribution at each
# grid point.
plot_fit <- function(K = 12L, labels = FALSE, component = "point", n_draws = 200L,
                     seed = 2026L) {
  key <- paste(K, labels, component, n_draws, seed)
  if (!is.null(plot_fit_cache[[key]])) {
    return(plot_fit_cache[[key]])
  }
  fit <- withr::with_seed(seed, {
    theta <- c(stats::rnorm(K %/% 3L, -1, 0.3), stats::rnorm(K - K %/% 3L, 1, 0.5))
    sigma <- stats::runif(K, 0.2, 0.9)
    theta_hat <- theta + stats::rnorm(K, 0, sigma)
    if (isTRUE(labels)) {
      names(theta_hat) <- sprintf("Site %s", formatC(seq_len(K), width = 3L, flag = "0"))
    }
    context <- plot_fit_ns(".bef_fit_prepare_context")(
      call = quote(bayes_efron_fit(theta_hat, sigma)), theta_hat = theta_hat, sigma = sigma,
      grid_method = "paper_realdata", component = component, L = 101L, seed = seed,
      chains = 2L, iter_sampling = n_draws %/% 2L, check_installed = FALSE
    )
    data <- context$stan_data
    # Coefficients whose g is close to a two-group shape, and draws around
    # them.
    target <- log(0.35 * stats::dnorm(data$grid, -1, 0.45) + 0.65 * stats::dnorm(data$grid, 1, 0.6) +
                    1e-4)
    center <- stats::coef(stats::lm.fit(cbind(1, data$B), target))[-1L]
    alpha <- matrix(center, nrow = n_draws, ncol = data$M, byrow = TRUE) +
      matrix(stats::rnorm(n_draws * data$M, sd = 0.25), nrow = n_draws)
    context$cmdstan_fit <- plot_fit_cmdstan(alpha, data)
    context$runtime_seconds <- 1.5
    context$sampler_settings <- plot_fit_ns(".bef_fit_sampler_settings")(
      context, interactive_fun = function() FALSE
    )
    suppressWarnings(plot_fit_ns(".bef_fit_assemble")(
      context, cmdstan_version_fun = function() "2.38.0"
    ))
  })
  plot_fit_cache[[key]] <- fit
  fit
}

# A stand-in for the object that CmdStan returns for the program of this
# version, for given draws of the coefficients.
plot_fit_cmdstan <- function(alpha, data, n_chains = 2L) {
  n_draws <- nrow(alpha)
  n_iter <- n_draws %/% n_chains
  log_w <- alpha %*% t(data$B)
  top <- apply(log_w, 1L, max)
  log_g <- log_w - (top + log(rowSums(exp(log_w - top))))
  g <- exp(log_g)
  mean_g <- drop(g %*% data$grid)
  var_g <- rowSums(g * outer(mean_g, data$grid, function(m, x) (x - m)^2)) + data$tau_c^2
  s_total <- sqrt(data$sigma^2 + data$tau_c^2)
  marginal <- vapply(seq_len(n_draws), function(s) {
    sum(vapply(seq_len(data$K), function(i) {
      log(sum(g[s, ] * stats::dnorm(data$theta_hat[[i]], data$grid, s_total[[i]])))
    }, numeric(1L)))
  }, numeric(1L))
  values <- cbind(-marginal, alpha, stats::rgamma(n_draws, 2, 2), log_g, mean_g, var_g,
                  sqrt(var_g), marginal)
  variables <- c("lp__", sprintf("alpha[%d]", seq_len(data$M)), "lambda",
                 sprintf("log_g[%d]", seq_len(data$L)), "mean_g", "var_g", "sd_g",
                 "log_marginal_likelihood")
  draws <- array(NA_real_, dim = c(n_iter, n_chains, length(variables)),
                 dimnames = list(NULL, NULL, variables))
  for (j in seq_along(variables)) {
    draws[, , j] <- matrix(values[, j], nrow = n_iter, ncol = n_chains)
  }
  nuts <- array(
    c(rep(0, n_draws), rep(3, n_draws), stats::rnorm(n_draws, 10), stats::runif(n_draws, 0.7, 1),
      rep(0.2, n_draws), rep(7, n_draws)),
    dim = c(n_iter, n_chains, 6L),
    dimnames = list(NULL, NULL, c("divergent__", "treedepth__", "energy__", "accept_stat__",
                                  "stepsize__", "n_leapfrog__"))
  )
  list(
    draws = function(variables = NULL, format = "draws_array") posterior::as_draws_array(draws),
    diagnostic_summary = function(...) {
      list(num_divergent = rep(0, n_chains), num_max_treedepth = rep(0, n_chains),
           ebfmi = rep(0.9, n_chains))
    },
    sampler_diagnostics = function(format = "draws_array") posterior::as_draws_array(nuts)
  )
}

# The data of a layer of a plot, and the classes of the geoms of its layers.
plot_layer_data <- function(plot, i) ggplot2::layer_data(plot, i)
plot_geoms <- function(plot) {
  vapply(plot$layers, function(layer) class(layer$geom)[[1L]], character(1L))
}

# The caption of a plot on one line.
plot_caption <- function(plot) gsub("\\s+", " ", plot$labels$caption)

# The labels on an axis of a panel of a plot: those of the breaks that fall
# in the range of the panel, from the bottom up (or from the left).
plot_axis_labels <- function(plot, axis = "y", panel = 1L) {
  view <- ggplot2::ggplot_build(plot)$layout$panel_params[[panel]][[axis]]
  labels <- view$get_labels()
  if (length(labels) == 0L) {
    return(character())
  }
  breaks <- view$get_breaks()
  shown <- !is.na(breaks)
  labels <- as.character(labels[shown])
  # The breaks of a discrete axis are already in its order.
  if (is.numeric(breaks)) labels[order(breaks[shown])] else labels
}

# A comparison of three ranges made without CmdStan: bayes_efron_fit() is
# replaced by a function that returns a small fit that depends on
# `expansion`, or fails for the value in `fail`.
plot_comparison <- function(fail = NULL, env = parent.frame()) {
  fits <- list(
    "0.25" = plot_fit(12L, labels = TRUE),
    "0.5" = plot_fit(12L, labels = TRUE, n_draws = 100L),
    "1" = plot_fit(12L, labels = TRUE, component = "normal")
  )
  testthat::local_mocked_bindings(
    bayes_efron_fit = function(...) {
      args <- list(...)
      if (!is.null(fail) && isTRUE(all.equal(args$expansion, fail))) {
        stop("the sampler stopped")
      }
      fits[[format(args$expansion)]]
    },
    .package = "bayesEfron",
    .env = env
  )
  data <- fits[[1L]]$metadata$data_list
  suppressWarnings(suppressMessages(compare_efron_grids(
    stats::setNames(data$theta_hat, site_summary(fits[[1L]])$label), data$sigma,
    L = 101L, M = 6L, expansion = c(0.25, 0.5, 1), sites = c(2, 9)
  )))
}
