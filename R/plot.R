#' Plot a fitted model
#'
#' @description
#' The `plot()` method for a fitted model shows the plot selected by
#' `type`. Each type is produced by one of the package's plot functions,
#' e.g., `"sites"` by [plot_sites()] and `"g"` by [plot_g()], and further
#' arguments are passed on to that function. See Details for all types.
#'
#' @details
#' `type` can be one of the following. The help page of each function
#' describes the plot and its arguments.
#'
#' * `"sites"` (the default): [plot_sites()], a caterpillar plot of the
#'   site effects with their posterior intervals
#' * `"g"`: [plot_g()], the estimated distribution of effects
#' * `"shrinkage"`: [plot_shrinkage()], a shrinkage plot of posterior means
#'   against estimates
#' * `"exceedance"`: [plot_exceedance()], the probability that each effect
#'   exceeds a threshold
#' * `"site_posterior"`: [plot_site_posterior()], the posterior
#'   distribution of the effect of each site specified in `sites`
#' * `"diagnostics"`: [plot_diagnostics()], sampler and grid diagnostics
#' * `"loo"`: [plot_loo()], the Pareto k values of the sites from
#'   leave-one-out cross-validation
#' * `"ppc"`: [plot_ppc()], a posterior predictive check of the model
#'   against the estimates
#' * `"mcmc"`: [plot_mcmc()], trace plots and other plots of the posterior
#'   draws
#' * `"grid"`: [plot.bef_grid()], the grid in relation to the estimates
#' * `"recovery"`: [plot_recovery()], a comparison of the fit with the
#'   known true effects specified in `truth`
#'
#' `"caterpillar"` and `"diagnostic"` are accepted as aliases for `"sites"`
#' and `"diagnostics"`.
#'
#' `type = "g"` shows the density of \eqn{g}. To plot the probabilities of
#' the grid points instead, call `plot_g(x, type = "mass")`. Note that
#' their heights depend on the number of grid points.
#'
#' @param x A `bef_fit` object created by [bayes_efron_fit()].
#' @param type The plot to show, e.g., `"sites"` (the default), `"g"` or
#'   `"shrinkage"`. See Details for all types.
#' @param ... Further arguments passed to the plot function for `type`,
#'   such as `level`, `sort_by` or `scale`.
#' @param backend Deprecated. `backend = "base"` plots the types
#'   `"caterpillar"`, `"g"`, `"diagnostic"` and `"loo"` with base graphics
#'   instead of ggplot2 and issues a warning. This option will be removed.
#'
#' @return A `ggplot` object, which can be modified with ggplot2 functions,
#'   e.g., by adding layers. For `type = "mcmc"`, the value returned by
#'   [plot_mcmc()]. With `backend = "base"`, the plot is shown on the
#'   current graphics device and `NULL` is returned invisibly.
#'
#' @seealso [plot_sites()], [plot_g()], [bef_fit], [diagnose()],
#'   `vignette("plots")`
#'
#' @examples
#' plot(krw_fit, n_extreme = 10, scale = 100)
#' plot(krw_fit, type = "g", scale = 100)
#' plot(krw_fit, type = "diagnostics")
#'
#' # The result is a ggplot object, so labels and layers can be added.
#' plot(krw_fit, type = "shrinkage", scale = 100) +
#'   ggplot2::labs(x = "Estimated gap (percentage points)")
#'
#' @export
plot.bef_fit_re <- function(x, type = "sites", ..., backend = NULL) {
  x <- .bef_prepare_fit(x)
  types <- .bef_plot_types()
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c(names(types), names(.bef_plot_type_aliases()))) {
    .bef_abort_invalid_args(
      sprintf("`type` must be one of: %s.",
              paste(sprintf("\"%s\"", names(types)), collapse = ", ")),
      arg = "type"
    )
  }
  args <- list(...)
  if (!is.null(backend)) {
    backend <- .bef_validate_method_choice(backend, c("ggplot2", "base"), arg = "backend")
    if (identical(backend, "base")) {
      return(.bef_plot_base_backend(x, type, args))
    }
  }
  old_name <- type %in% c(names(.bef_plot_type_aliases()), "g", "loo")
  if (type %in% names(.bef_plot_type_aliases())) {
    type <- .bef_plot_type_aliases()[[type]]
  }
  if (identical(type, "grid")) {
    return(.bef_plot_fit_grid(x, args))
  }
  target <- types[[type]]
  if (old_name) {
    # Version 0.3 took `level` and `sort_by` for every type and ignored them
    # where they had no meaning.
    ignored <- setdiff(c("level", "sort_by"), names(formals(target)))
    args <- args[!names(args) %in% ignored | !nzchar(names(args))]
  }
  do.call(target, c(list(x), args))
}

# The types of plot() and the functions that draw them.
.bef_plot_types <- function() {
  list(
    sites = plot_sites, g = plot_g, shrinkage = plot_shrinkage,
    exceedance = plot_exceedance, diagnostics = plot_diagnostics, loo = plot_loo,
    grid = NULL, site_posterior = plot_site_posterior, ppc = plot_ppc,
    recovery = plot_recovery, mcmc = plot_mcmc
  )
}

# The names that version 0.3 gave to two of the types.
.bef_plot_type_aliases <- function() {
  c(caterpillar = "sites", diagnostic = "diagnostics")
}

# The grid of a fit, drawn with the data of the fit.
.bef_plot_fit_grid <- function(x, args) {
  grid <- .bef_metadata_attr(x$metadata, "grid")
  if (!inherits(grid, "bef_grid")) {
    .bef_abort_invalid_fit(
      paste("This fit does not hold its grid as an object, which fits made by bayesEfron 0.4",
            "or later do. Build the grid with `make_efron_grid()` and plot that.")
    )
  }
  data <- x$metadata$data_list
  do.call(plot, c(list(grid, theta_hat = unname(data$theta_hat), sigma = data$sigma), args))
}

# The four plots of version 0.3 in base graphics. The code below this
# function is that of version 0.3 and is not developed further.
.bef_plot_base_backend <- function(x, type, args) {
  old <- c("caterpillar", "g", "diagnostic", "loo")
  if (identical(type, "sites")) {
    type <- "caterpillar"
  }
  if (identical(type, "diagnostics")) {
    type <- "diagnostic"
  }
  if (!type %in% old) {
    .bef_abort_invalid_args(
      sprintf("`backend = \"base\"` draws only the types %s; `type = \"%s\"` is drawn with ggplot2.",
              paste(sprintf("\"%s\"", old), collapse = ", "), type),
      arg = "backend"
    )
  }
  extra <- setdiff(names(args), c("level", "sort_by"))
  if (length(args) > 0L && (is.null(names(args)) || length(extra) > 0L)) {
    .bef_abort_invalid_args(
      "With `backend = \"base\"` the only further arguments are `level` and `sort_by`.",
      arg = "..."
    )
  }
  level <- .bef_validate_summary_level(if (is.null(args$level)) 0.9 else args$level)
  sort_by <- .bef_validate_method_choice(
    if (is.null(args$sort_by)) "mean" else args$sort_by, c("mean", "sigma", "none"),
    arg = "sort_by"
  )
  .bef_warn_base_backend()
  payload <- .bef_plot_payload_bef_fit_re(x = x, type = type, level = level, sort_by = sort_by)
  .bef_plot_base(payload)
  invisible(NULL)
}

# Said once in a session.
.bef_warn_base_backend <- function() {
  key <- "plot_backend_base"
  if (isTRUE(.bayesEfron_msgs_emitted[[key]])) {
    return(invisible(FALSE))
  }
  .bayesEfron_msgs_emitted[[key]] <- TRUE
  .bef_warn(
    paste("`backend = \"base\"` is deprecated and will be removed in bayesEfron 0.5. The plots",
          "are drawn with ggplot2; `g_summary()` and `site_summary()` give their numbers for",
          "plots of your own. This warning is shown once in a session."),
    "bef_deprecated"
  )
  invisible(TRUE)
}

.bef_plot_payload_bef_fit_re <- function(x, type, level, sort_by) {
  switch(
    type,
    caterpillar = list(
      type = type,
      level = level,
      sort_by = sort_by,
      reference = x$metadata$mean_g_summary$mean,
      data = .bef_caterpillar_data(x, level = level, sort_by = sort_by)
    ),
    g = list(
      type = type,
      level = level,
      data = .bef_g_plot_data(x, level = level)
    ),
    diagnostic = list(
      type = type,
      level = level,
      data = .bef_diagnostic_plot_data(summary(x)$diagnostics)
    ),
    loo = list(
      type = type,
      level = level,
      data = .bef_loo_plot_data(x)
    )
  )
}

.bef_caterpillar_data <- function(x, level = 0.9, sort_by = "mean") {
  # The plot keeps the column names of version 0.3.0.
  theta <- as.data.frame(confint(x, level = level, type = "theta"))
  theta <- theta[, c("site", "lower", "upper", "point")]
  engine <- .bef_engine(x)
  inner <- .bef_site_quantile(engine, c(0.25, 0.75))
  theta$inner_lower <- inner[, 1L]
  theta$inner_upper <- inner[, 2L]
  theta$sd <- x$metadata$theta_summary$sd
  theta$sigma <- x$metadata$data_list$sigma
  ord <- switch(
    sort_by,
    mean = order(theta$point),
    sigma = order(theta$sigma),
    none = seq_len(nrow(theta))
  )
  theta <- theta[ord, , drop = FALSE]
  theta$position <- seq_len(nrow(theta))
  theta
}

.bef_g_plot_data <- function(x, level) {
  density <- .bef_g_density_plot_data(x, level = level)
  if (!is.null(density)) {
    return(density)
  }
  data <- as.data.frame(confint(x, level = level, type = "g"))
  data <- data[, c("site", "lower", "upper", "point")]
  data$kind <- "moment"
  data$position <- seq_len(nrow(data))
  data
}

.bef_g_density_plot_data <- function(x, level) {
  # g is stored only when the fit was run with store_grid_quantities = TRUE;
  # otherwise it is recovered from log_g. A fit without either has no density
  # to draw.
  g_draws <- tryCatch(
    .bef_plain_draw_matrix(.bef_grid_draws(x, "g")),
    error = function(err) NULL
  )
  if (is.null(g_draws)) {
    return(NULL)
  }
  probs <- .bef_interval_probs(level)
  intervals <- apply(
    g_draws,
    2L,
    function(value) {
      posterior::quantile2(value, probs = probs, names = FALSE)
    }
  )
  data.frame(
    kind = "density",
    grid = x$metadata$data_list$grid,
    lower = intervals[1L, ],
    upper = intervals[2L, ],
    point = colMeans(g_draws),
    row.names = NULL,
    check.names = FALSE
  )
}

# One value per diagnostic for the diagnostic plot: the largest Rhat, the
# smallest effective sample sizes and the total counts, as in
# summary.bef_diagnostic().
.bef_diagnostic_plot_data <- function(diagnostics) {
  worst <- .bef_diagnostic_worst
  data.frame(
    metric = c("rhat", "ess_bulk", "ess_tail", "divergences", "max_treedepth"),
    value = c(
      worst(diagnostics$rhat, max),
      worst(diagnostics$ess_bulk, min),
      worst(diagnostics$ess_tail, min),
      worst(diagnostics$divergences, sum),
      worst(diagnostics$max_treedepth, sum)
    ),
    row.names = NULL,
    check.names = FALSE
  )
}

.bef_plot_base <- function(payload) {
  switch(
    payload$type,
    caterpillar = .bef_plot_caterpillar_base(payload),
    g = .bef_plot_g_base(payload),
    diagnostic = .bef_plot_diagnostic_base(payload),
    loo = .bef_plot_loo_base(payload)
  )
  invisible(NULL)
}

.bef_plot_caterpillar_base <- function(payload, main = "Site effects") {
  data <- payload$data
  xlim <- range(c(data$lower, data$upper, data$inner_lower, data$inner_upper, payload$reference))
  graphics::plot(
    NA,
    xlim = xlim,
    ylim = c(0.5, nrow(data) + 0.5),
    xlab = "Effect",
    ylab = "Site",
    yaxt = "n",
    main = main
  )
  graphics::axis(2, at = data$position, labels = data$site, las = 1)
  graphics::segments(data$lower, data$position, data$upper, data$position)
  graphics::segments(
    data$inner_lower, data$position, data$inner_upper, data$position,
    lwd = 3
  )
  graphics::points(data$point, data$position, pch = 19)
  graphics::abline(v = payload$reference, lty = 2, col = "grey50")
}

.bef_plot_g_base <- function(payload) {
  data <- payload$data
  if (identical(data$kind[[1L]], "density")) {
    ylim <- range(c(data$lower, data$upper))
    graphics::plot(
      data$grid,
      data$point,
      type = "l",
      ylim = ylim,
      xlab = "Effect",
      ylab = "Probability",
      main = "Distribution of effects"
    )
    graphics::lines(data$grid, data$lower, lty = 2, col = "grey50")
    graphics::lines(data$grid, data$upper, lty = 2, col = "grey50")
    return(invisible(NULL))
  }

  ylim <- range(c(data$lower, data$upper))
  graphics::plot(
    data$position,
    data$point,
    ylim = ylim,
    xaxt = "n",
    xlab = "",
    ylab = "Posterior interval",
    pch = 19,
    main = "Summaries of the distribution of effects"
  )
  graphics::axis(1, at = data$position, labels = data$site)
  graphics::segments(data$position, data$lower, data$position, data$upper)
}

.bef_plot_diagnostic_base <- function(payload) {
  data <- payload$data
  old_par <- graphics::par(mfrow = c(2, 2))
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::barplot(
    data$value[data$metric == "rhat"],
    names.arg = "Rhat",
    ylab = "max",
    main = "Rhat"
  )
  graphics::barplot(
    data$value[data$metric %in% c("ess_bulk", "ess_tail")],
    names.arg = c("bulk", "tail"),
    ylab = "min",
    main = "ESS"
  )
  graphics::barplot(
    data$value[data$metric == "divergences"],
    names.arg = "divergences",
    ylab = "count",
    main = "Divergences"
  )
  graphics::barplot(
    data$value[data$metric == "max_treedepth"],
    names.arg = "treedepth",
    ylab = "count",
    main = "Max treedepth hits"
  )
}

# Pareto k by site for the loo plot. The bands at 0.5 and 0.7 follow the loo
# package.
.bef_loo_plot_data <- function(x) {
  .bef_require_loo()
  fit_loo <- loo::loo(x)
  k <- fit_loo$diagnostics$pareto_k
  data.frame(
    site = seq_along(k),
    pareto_k = as.numeric(k),
    band = cut(
      as.numeric(k),
      breaks = c(-Inf, 0.5, 0.7, Inf),
      labels = c("good", "ok", "bad"),
      right = TRUE
    ),
    row.names = NULL
  )
}

.bef_plot_loo_base <- function(payload, ...) {
  data <- payload$data
  # The range takes in the two reference lines, as the ggplot2 version does.
  graphics::plot(
    data$site, data$pareto_k,
    ylim = range(c(data$pareto_k, 0.5, 0.7), finite = TRUE),
    xlab = "Site", ylab = "Pareto k",
    main = "Pareto k by site",
    pch = 16, ...
  )
  graphics::abline(h = c(0.5, 0.7), lty = c(3L, 2L))
  invisible(data)
}
