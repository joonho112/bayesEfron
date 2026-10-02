#' Plot a fitted model
#'
#' @description
#' `plot()` draws one of four plots of a fitted model: the site effects
#' with their intervals, the estimated distribution of effects, the
#' convergence diagnostics, or the Pareto k values of leave-one-out
#' cross-validation.
#'
#' @param x A `bef_fit_re` object from [bayes_efron_fit()].
#' @param type Character string, the plot to draw. `"caterpillar"`
#'   (the default) shows for each site the posterior mean of
#'   \eqn{\theta_i}, an interval of probability `level` as a thin line and
#'   the central 50% interval as a thick line; a dashed vertical line marks
#'   the posterior mean of the mean of \eqn{g}. `"g"` shows the posterior
#'   mean of the probability that \eqn{g} gives to each grid point, with a
#'   pointwise band of probability `level`. `"diagnostic"` shows the
#'   largest R-hat, the smallest bulk and tail effective sample sizes and
#'   the numbers of divergent transitions and of iterations at the maximum
#'   tree depth, as [diagnose()] reports them; E-BFMI is not drawn.
#'   `"loo"` shows the Pareto k value of each site from [loo.bef_fit()],
#'   with reference lines at 0.5 and 0.7; it needs the loo package.
#' @param level Number between 0 and 1, the probability content of the
#'   intervals and of the band. Defaults to `0.9`.
#' @param sort_by Character string, the order of the sites in the
#'   caterpillar plot: `"mean"` (the default) for the order of the
#'   posterior means, `"sigma"` for the order of the standard errors, or
#'   `"none"` for the order of the data.
#' @param ... Not used.
#' @param backend Character string, `"ggplot2"` or `"base"`. If it is not
#'   given, ggplot2 is used when it is installed and base graphics
#'   otherwise.
#'
#' @return With ggplot2, a `ggplot` object, to which layers can be added.
#'   With base graphics the plot is drawn and `NULL` is returned invisibly.
#'
#' @seealso [bef_fit], [diagnose()], [loo.bef_fit()]
#'
#' @examples
#' plot(raudenbush_fit)
#' plot(raudenbush_fit, sort_by = "sigma", level = 0.95)
#' plot(raudenbush_fit, type = "g")
#' plot(raudenbush_fit, type = "diagnostic", backend = "base")
#'
#' # With ggplot2 the result can be extended.
#' p <- plot(raudenbush_fit, type = "g")
#' if (inherits(p, "ggplot")) {
#'   p + ggplot2::labs(title = "Teacher expectancy effects")
#' }
#'
#' @examplesIf requireNamespace("loo", quietly = TRUE)
#' plot(raudenbush_fit, type = "loo")
#'
#' @export
plot.bef_fit_re <- function(x,
                            type = c("caterpillar", "g", "diagnostic", "loo"),
                            level = 0.9,
                            sort_by = c("mean", "sigma", "none"),
                            ...,
                            backend = c("ggplot2", "base")) {
  x <- .bef_prepare_fit(x)
  type <- .bef_validate_method_choice(
    type,
    c("caterpillar", "g", "diagnostic", "loo"),
    arg = "type"
  )
  level <- .bef_validate_summary_level(level)
  sort_by <- .bef_validate_method_choice(
    sort_by,
    c("mean", "sigma", "none"),
    arg = "sort_by"
  )

  payload <- .bef_plot_payload_bef_fit_re(
    x = x,
    type = type,
    level = level,
    sort_by = sort_by
  )

  if (.bef_plot_use_ggplot2(backend, explicit = !missing(backend))) {
    return(.bef_plot_ggplot2(payload))
  }
  .bef_plot_base(payload)
  invisible(NULL)
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
  theta <- confint(x, level = level, type = "theta")
  theta_rep <- .bef_site_draws(x, "theta_rep")
  inner <- vapply(
    seq_len(ncol(theta_rep)),
    function(site) {
      posterior::quantile2(
        theta_rep[, site],
        probs = c(0.25, 0.75),
        names = FALSE
      )
    },
    numeric(2L)
  )
  theta$inner_lower <- inner[1L, ]
  theta$inner_upper <- inner[2L, ]
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
  data <- confint(x, level = level, type = "g")
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

.bef_plot_use_ggplot2 <- function(backend = c("ggplot2", "base"),
                                  explicit = FALSE) {
  backend <- .bef_validate_method_choice(
    backend, c("ggplot2", "base"), arg = "backend"
  )
  if (identical(backend, "base")) {
    return(FALSE)
  }
  if (isTRUE(explicit)) {
    rlang::check_installed("ggplot2", reason = "for `backend = \"ggplot2\"`")
  }
  requireNamespace("ggplot2", quietly = TRUE)
}

.bef_plot_ggplot2 <- function(payload) {
  plot <- switch(
    payload$type,
    caterpillar = .bef_plot_caterpillar_ggplot2(payload),
    g = .bef_plot_g_ggplot2(payload),
    diagnostic = .bef_plot_diagnostic_ggplot2(payload),
    loo = .bef_plot_loo_ggplot2(payload)
  )
  attr(plot, "bef_plot_payload") <- payload
  plot
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

.bef_plot_caterpillar_ggplot2 <- function(payload,
                                          title = "Site effects") {
  data <- payload$data
  ggplot2::ggplot(data, ggplot2::aes(y = .data$position)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$lower, xend = .data$upper, yend = .data$position)
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$inner_lower,
        xend = .data$inner_upper,
        yend = .data$position
      ),
      linewidth = 1.1
    ) +
    ggplot2::geom_point(ggplot2::aes(x = .data$point), size = 1.8) +
    ggplot2::geom_vline(xintercept = payload$reference, linetype = 2, color = "grey50") +
    ggplot2::scale_y_continuous(breaks = data$position, labels = data$site) +
    ggplot2::labs(
      title = title,
      x = "Effect",
      y = "Site"
    ) +
    ggplot2::theme_minimal()
}

.bef_plot_g_ggplot2 <- function(payload) {
  if (identical(payload$data$kind[[1L]], "density")) {
    return(
      ggplot2::ggplot(
        payload$data,
        ggplot2::aes(x = .data$grid, y = .data$point)
      ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
          alpha = 0.2
        ) +
        ggplot2::geom_line() +
        ggplot2::labs(
          title = "Distribution of effects",
          x = "Effect",
          y = "Probability"
        ) +
        ggplot2::theme_minimal()
    )
  }

  ggplot2::ggplot(
    payload$data,
    ggplot2::aes(x = .data$site, y = .data$point)
  ) +
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper)
    ) +
    ggplot2::labs(
      title = "Summaries of the distribution of effects",
      x = NULL,
      y = "Posterior interval"
    ) +
    ggplot2::theme_minimal()
}

.bef_plot_diagnostic_ggplot2 <- function(payload) {
  ggplot2::ggplot(
    payload$data,
    ggplot2::aes(x = .data$metric, y = .data$value)
  ) +
    ggplot2::geom_col() +
    # The diagnostics are on different scales, so each has its own panel.
    ggplot2::facet_wrap(ggplot2::vars(.data$metric), scales = "free", nrow = 1L) +
    ggplot2::labs(
      title = "Convergence diagnostics",
      x = NULL,
      y = "Value"
    ) +
    ggplot2::theme_minimal()
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

.bef_plot_loo_ggplot2 <- function(payload) {
  data <- payload$data
  ggplot2::ggplot(data, ggplot2::aes(x = .data$site, y = .data$pareto_k)) +
    ggplot2::geom_hline(yintercept = c(0.5, 0.7), linetype = c(3L, 2L)) +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Site", y = "Pareto k", title = "Pareto k by site"
    )
}
