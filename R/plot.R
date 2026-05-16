#' Plot a bayesEfron random-effects fit
#'
#' @description
#' Draw one of four diagnostic-or-summary views of a fitted
#' `bef_fit_re` object. The four views support visual answers to the
#' four most common questions a meta-analyst asks of a deconvolution
#' fit: where do the per-site effects sit (caterpillar), what does
#' the underlying mixing distribution look like (prior summary), how
#' sensitive are the per-site intervals to alternative credible
#' levels (sensitivity), and is the sampler healthy (diagnostic).
#'
#' Both `ggplot2` and base-graphics backends are available. The
#' `ggplot2` backend is used when the package is installed and the
#' environment variable `BAYESEFRON_NO_GGPLOT2` is not exactly `"1"`;
#' otherwise the base backend is used.
#'
#' @details
#' # Plot types
#'
#' | `type` | Shows | Notes |
#' |:-------|:------|:------|
#' | `"caterpillar"` | Per-site posterior means with credible intervals. The intervals are at the requested `level`. | Use `sort_by = "mean"` (default) to order by point estimate, `"sigma"` to order by within-study uncertainty, or `"none"` to keep input order. |
#' | `"g"` | Posterior summary of the mixing distribution \eqn{g}: posterior mean ± credible band over the discrete grid. | Reads `mean_g_summary` and `var_g_summary` from `fit$metadata`; reads `sd_g_summary` from the `sd_g_summary` attribute when available. |
#' | `"sensitivity"` | Caterpillar overlay at the requested `level` against the default 90% intervals. | Useful for checking interval-width sensitivity to `level`. |
#' | `"diagnostic"` | Sampler-diagnostic summary view: `rhat`, `ess_bulk`, `ess_tail`, divergences, max-treedepth flags. | Reads from `fit$metadata`'s diagnostic attribute. Pairs with [diagnose()]. |
#'
#' # Backend selection
#'
#' By default, the function returns a `ggplot2::ggplot` object when
#' `ggplot2` is available, allowing downstream `+ theme(...)` and
#' `+ labs(...)` chaining. Setting `BAYESEFRON_NO_GGPLOT2=1` in the
#' environment forces the base-graphics path, which writes to the
#' active graphics device and returns `invisible(NULL)`.
#'
#' Code that needs to chain `ggplot2` layers should guard against the
#' base path with `inherits(p, "ggplot")` before applying `+ theme()`
#' or `+ labs()` (see the `@examples` block).
#'
#' @param x A `bef_fit_re` object returned by [bayes_efron_fit()].
#' @param type Plot type. One of `"caterpillar"` (default), `"g"`,
#'   `"sensitivity"`, or `"diagnostic"`.
#' @param level Numeric credible level in \eqn{(0, 1)}. Defaults to
#'   `0.9`.
#' @param sort_by Caterpillar ordering. One of `"mean"` (default),
#'   `"sigma"`, or `"none"`. Ignored by `"g"` and `"diagnostic"`
#'   views.
#' @param ... Reserved for future expansion; must be empty in v0.1.
#'
#' @return Backend-dependent:
#'
#'   * **ggplot2 backend** (`ggplot2` installed and
#'     `BAYESEFRON_NO_GGPLOT2 != "1"`): returns a
#'     `ggplot2::ggplot` object visibly so it can be assigned and
#'     extended.
#'   * **Base backend**: draws on the active graphics device and
#'     returns `invisible(NULL)`.
#'
#' @seealso
#'   * [bayes_efron_fit()] for producing the fit.
#'   * [diagnose()] for a structured (non-graphical) diagnostic
#'     producer that pairs with `type = "diagnostic"`.
#'   * [confint.bef_fit_re()][bayesEfron-methods] for the numeric
#'     intervals shown in the caterpillar view.
#'
#' @examples
#' # Load the cached five-site smoke fit shipped with the package.
#' fit <- readRDS(system.file(
#'   "examples", "cached_fit_re_smoke.rds",
#'   package = "bayesEfron"
#' ))
#'
#' \donttest{
#' plot(fit, type = "caterpillar")
#' plot(fit, type = "caterpillar", sort_by = "sigma", level = 0.95)
#' plot(fit, type = "g")
#' plot(fit, type = "sensitivity", level = 0.95)
#' plot(fit, type = "diagnostic")
#'
#' # Chain ggplot2 layers when the ggplot2 backend is active.
#' p <- plot(fit, type = "caterpillar")
#' if (inherits(p, "ggplot") && requireNamespace("ggplot2", quietly = TRUE)) {
#'   p + ggplot2::theme_minimal() +
#'       ggplot2::labs(title = "Per-site posterior intervals")
#' }
#' }
#'
#' @export
plot.bef_fit_re <- function(x,
                            type = c("caterpillar", "g", "sensitivity", "diagnostic"),
                            level = 0.9,
                            sort_by = c("mean", "sigma", "none"),
                            ...) {
  type <- .bef_validate_method_choice(
    type,
    c("caterpillar", "g", "sensitivity", "diagnostic"),
    arg = "type",
    module = "plot.bef_fit_re"
  )
  level <- .bef_validate_summary_level(level)
  sort_by <- .bef_validate_method_choice(
    sort_by,
    c("mean", "sigma", "none"),
    arg = "sort_by",
    module = "plot.bef_fit_re"
  )

  payload <- .bef_plot_payload_bef_fit_re(
    x = x,
    type = type,
    level = level,
    sort_by = sort_by
  )

  if (.bef_plot_use_ggplot2()) {
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
    sensitivity = list(
      type = type,
      level = level,
      sort_by = sort_by,
      reference = x$metadata$mean_g_summary$mean,
      data = .bef_caterpillar_data(x, level = level, sort_by = sort_by)
    ),
    diagnostic = list(
      type = type,
      level = level,
      data = .bef_diagnostic_plot_data(summary(x)$diagnostics)
    )
  )
}

.bef_caterpillar_data <- function(x, level = 0.9, sort_by = "mean") {
  theta <- confint(x, level = level, type = "theta")
  inner <- vapply(
    seq_len(ncol(x$metadata$theta_rep_draws)),
    function(site) {
      posterior::quantile2(
        x$metadata$theta_rep_draws[, site],
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
  variables <- dimnames(x$draws)[[3L]]
  L <- as.integer(x$metadata$data_list$L)
  fields <- paste0("g[", seq_len(L), "]")
  if (!is.character(variables) || !all(fields %in% variables)) {
    return(NULL)
  }

  draw_matrix <- posterior::as_draws_matrix(x$draws)
  g_draws <- as.matrix(draw_matrix[, fields, drop = FALSE])
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

.bef_diagnostic_plot_data <- function(diagnostics) {
  data.frame(
    metric = c("rhat", "ess_bulk", "ess_tail", "divergences", "max_treedepth"),
    value = c(
      diagnostics$rhat,
      diagnostics$ess_bulk,
      diagnostics$ess_tail,
      diagnostics$divergences,
      diagnostics$max_treedepth
    ),
    row.names = NULL,
    check.names = FALSE
  )
}

.bef_plot_use_ggplot2 <- function() {
  requireNamespace("ggplot2", quietly = TRUE) &&
    !identical(Sys.getenv("BAYESEFRON_NO_GGPLOT2"), "1")
}

.bef_plot_ggplot2 <- function(payload) {
  plot <- switch(
    payload$type,
    caterpillar = .bef_plot_caterpillar_ggplot2(payload),
    g = .bef_plot_g_ggplot2(payload),
    sensitivity = .bef_plot_caterpillar_ggplot2(
      payload,
      title = "bayesEfron sensitivity"
    ),
    diagnostic = .bef_plot_diagnostic_ggplot2(payload)
  )
  attr(plot, "bef_plot_payload") <- payload
  plot
}

.bef_plot_base <- function(payload) {
  switch(
    payload$type,
    caterpillar = .bef_plot_caterpillar_base(payload),
    g = .bef_plot_g_base(payload),
    sensitivity = .bef_plot_caterpillar_base(
      payload,
      main = "bayesEfron sensitivity"
    ),
    diagnostic = .bef_plot_diagnostic_base(payload)
  )
  invisible(NULL)
}

.bef_plot_caterpillar_base <- function(payload, main = "bayesEfron caterpillar") {
  data <- payload$data
  xlim <- range(c(data$lower, data$upper, data$inner_lower, data$inner_upper, payload$reference))
  graphics::plot(
    NA,
    xlim = xlim,
    ylim = c(0.5, nrow(data) + 0.5),
    xlab = "theta",
    ylab = "site",
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
      xlab = "theta",
      ylab = "g(theta)",
      main = "bayesEfron prior density"
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
    xlab = "prior moment",
    ylab = "posterior interval",
    pch = 19,
    main = "bayesEfron prior summaries"
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
    main = "Max treedepth"
  )
}

.bef_plot_caterpillar_ggplot2 <- function(payload,
                                          title = "bayesEfron caterpillar") {
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
    ggplot2::geom_vline(xintercept = payload$reference, linetype = 2, colour = "grey50") +
    ggplot2::scale_y_continuous(breaks = data$position, labels = data$site) +
    ggplot2::labs(
      title = title,
      x = "theta",
      y = "site"
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
          title = "bayesEfron prior density",
          x = "theta",
          y = "g(theta)"
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
      title = "bayesEfron prior summaries",
      x = "prior moment",
      y = "posterior interval"
    ) +
    ggplot2::theme_minimal()
}

.bef_plot_diagnostic_ggplot2 <- function(payload) {
  ggplot2::ggplot(
    payload$data,
    ggplot2::aes(x = .data$metric, y = .data$value)
  ) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title = "bayesEfron diagnostics",
      x = "diagnostic",
      y = "value"
    ) +
    ggplot2::theme_minimal()
}
