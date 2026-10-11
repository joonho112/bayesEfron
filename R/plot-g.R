#' Density and CDF plots of the estimated distribution of effects
#'
#' @description
#' `plot_g()` plots the estimate of \eqn{g}, the distribution of the true
#' effects: its density, its cumulative distribution function (CDF) or its
#' probabilities on the grid. The posterior mean is shown as a line and a
#' pointwise posterior interval as a band. The estimates, a normal
#' distribution, known true effects, individual posterior draws of \eqn{g}
#' and other estimates of \eqn{g} can be added to the plot.
#'
#' @details
#' The line and the band are the columns `mean`, `lower` and `upper`
#' returned by [g_summary()]. For point masses, the density is the
#' probability of a grid point divided by the grid spacing, and the CDF
#' spreads that probability evenly over the cell around the grid point.
#'
#' The band is pointwise: the stated posterior probability applies to the
#' curve at each effect value separately, not to the whole curve at once.
#'
#' `normal = TRUE` adds the distribution of effects under the normal
#' random-effects model, with its mean and standard deviation estimated by
#' restricted maximum likelihood (REML) from the same estimates. This is a
#' plug-in estimate, shown as a single curve without a band, because the
#' uncertainty in the estimated mean and standard deviation is ignored.
#'
#' Note that a histogram of the estimates is not an estimate of \eqn{g}.
#' The estimates include sampling error, so they are more spread out than
#' the true effects. The histogram shows the data to which the model was
#' fitted.
#'
#' The horizontal axis covers the grid unless the grid is much wider than
#' the region where \eqn{g} has its mass, which would leave long empty
#' stretches on both sides. If the central 99.9% of the posterior mean of
#' \eqn{g} covers less than 60% of the grid, the axis is limited to that
#' region plus a small margin. Use `xlim` to set other limits.
#'
#' @param x A `bef_fit` object created by [bayes_efron_fit()].
#' @param type What to plot: `"density"` (the default), `"cdf"` for the
#'   cumulative distribution function, or `"mass"` for the probabilities of
#'   the grid points. With `"mass"`, only a rug of the estimates and
#'   posterior draws of \eqn{g} can be added to the plot.
#' @param level Posterior probability of the pointwise band, a number
#'   between 0 and 1.
#' @param inner_level Posterior probability of a second, darker band inside
#'   the first, a number between 0 and `level`. If `NULL` (the default), no
#'   inner band is shown.
#' @param observed How to show the estimates: `"none"` (the default),
#'   `"histogram"` or `"rug"`. When `type = "cdf"`, `"histogram"` shows
#'   their empirical distribution function.
#' @param normal If `TRUE`, the normal distribution described in Details is
#'   added. Can also be a numeric vector `c(mean = , sd = )`, which adds
#'   the normal distribution with that mean and standard deviation, in the
#'   units of the estimates.
#' @param truth For simulated data, the true effects. A numeric vector is
#'   shown as a histogram (with `type = "cdf"`, as an empirical
#'   distribution function). Can also be a function of the effect that
#'   returns the true density (with `type = "cdf"`, the true distribution
#'   function). If `NULL` (the default), no true effects are shown.
#' @param draws Number of individual posterior draws of \eqn{g} to add as
#'   thin lines, in addition to the posterior mean and the band. At most
#'   the number of draws in the fit can be shown.
#' @param compare An optional named list of other estimates of \eqn{g},
#'   each plotted as one line. An element can be a `bef_fit` object or a
#'   data frame with a column `x` of effect values, in the units of the
#'   estimates, and a column `density` or `cdf`.
#' @param xlim Limits of the horizontal axis, a numeric vector of length 2
#'   in the units of the plot, i.e., after multiplication by `scale`. If
#'   `NULL` (the default), the limits are chosen as described in Details.
#' @param scale Positive number by which effect values are multiplied and
#'   densities divided, e.g., `scale = 100` to show differences in
#'   proportions as percentage points.
#' @param ... These dots are for future extensions and must be empty.
#'
#' @return A `ggplot` object.
#'
#' @seealso [g_summary()], [plot_sites()], `vignette("plots")`
#'
#' @examples
#' # The callback gaps in percentage points, with a histogram of the
#' # estimates and the normal distribution of the random-effects model.
#' plot_g(krw_fit, observed = "histogram", normal = TRUE, scale = 100)
#'
#' # The distribution function, with a 50% band inside the 90% band.
#' plot_g(krw_fit, type = "cdf", inner_level = 0.5, scale = 100)
#'
#' # Simulated sites with known true effects.
#' plot_g(twin_towers_fit, truth = twin_towers$theta_true, normal = TRUE)
#'
#' # The posterior mean and band with 40 posterior draws of g.
#' plot_g(raudenbush_fit, draws = 40)
#'
#' @export
plot_g <- function(x,
                   type = c("density", "cdf", "mass"),
                   level = 0.9,
                   inner_level = NULL,
                   observed = c("none", "histogram", "rug"),
                   normal = FALSE,
                   truth = NULL,
                   draws = 0,
                   compare = NULL,
                   xlim = NULL,
                   scale = 1,
                   ...) {
  .bef_check_dots_empty(list(...), "plot_g")
  x <- .bef_plot_fit(x)
  type <- .bef_validate_method_choice(type, c("density", "cdf", "mass"), arg = "type")
  observed <- .bef_validate_method_choice(observed, c("none", "histogram", "rug"),
                                          arg = "observed")
  level <- .bef_plot_level(level)
  inner_level <- .bef_plot_level(inner_level, "inner_level", null_ok = TRUE)
  if (!is.null(inner_level) && inner_level >= level) {
    .bef_abort_invalid_args("`inner_level` must be smaller than `level`.", arg = "inner_level")
  }
  draws <- .bef_plot_count(draws, "draws")
  scale <- .bef_validate_scale(scale)
  engine <- .bef_engine(x)
  normal <- .bef_plot_normal(normal, engine$theta_hat, engine$sigma)
  if (identical(type, "mass") &&
      (identical(observed, "histogram") || !is.null(normal) || !is.null(truth) ||
         !is.null(compare))) {
    .bef_abort_invalid_args(
      "With `type = \"mass\"` nothing but a rug of the estimates and draws of g can be laid over the plot: the height of the probabilities depends on the number of grid points. Use `type = \"density\"`.",
      arg = "type"
    )
  }
  if (!is.null(xlim) && !(is.numeric(xlim) && length(xlim) == 2L && all(is.finite(xlim)) &&
                            xlim[1L] < xlim[2L])) {
    .bef_abort_invalid_args("`xlim` must be two finite numbers in increasing order, or `NULL`.",
                            arg = "xlim")
  }

  colors <- .bef_plot_colors()
  at <- .bef_plot_g_values(engine, type)
  band <- g_summary(x, type = type, level = level, at = at$values, scale = scale)
  plot <- ggplot2::ggplot()

  # What lies under the estimate of g: histograms of the true effects and
  # of the estimates. With the distribution function they are step lines,
  # which join the other lines further down.
  fills <- character()
  steps <- list()
  if (is.numeric(truth)) {
    if (length(truth) == 0L || any(!is.finite(truth))) {
      .bef_abort_invalid_args("`truth` must hold finite numbers.", arg = "truth")
    }
    if (identical(type, "cdf")) {
      steps[[length(steps) + 1L]] <- .bef_plot_g_ecdf(truth * scale, "True effects")
    } else {
      plot <- plot + .bef_plot_g_histogram(truth * scale, "True effects")
      fills <- c(fills, "True effects")
    }
  } else if (!is.null(truth) && !is.function(truth)) {
    .bef_abort_invalid_args("`truth` must be a numeric vector, a function, or `NULL`.",
                            arg = "truth")
  }
  if (identical(observed, "histogram")) {
    if (identical(type, "cdf")) {
      steps[[length(steps) + 1L]] <- .bef_plot_g_ecdf(engine$theta_hat * scale, "Estimates")
    } else {
      plot <- plot + .bef_plot_g_histogram(engine$theta_hat * scale, "Estimates",
                                           outline = "True effects" %in% fills)
      fills <- c(fills, "Estimates")
    }
  }
  if (length(fills) > 0L) {
    plot <- plot +
      ggplot2::scale_fill_manual(
        values = c("True effects" = "#B3B3B3", "Estimates" = "#B3B3B3")[fills],
        breaks = fills
      )
  }

  # The bands and single draws of g.
  plot <- plot +
    ggplot2::geom_ribbon(
      data = band,
      ggplot2::aes(x = .data$value, ymin = .data$lower, ymax = .data$upper),
      fill = colors[["fit"]], alpha = 0.2
    )
  if (!is.null(inner_level)) {
    inner <- g_summary(x, type = type, level = inner_level, at = at$values, scale = scale)
    plot <- plot +
      ggplot2::geom_ribbon(
        data = inner,
        ggplot2::aes(x = .data$value, ymin = .data$lower, ymax = .data$upper),
        fill = colors[["fit"]], alpha = 0.25
      )
  }
  drawn_count <- 0L
  if (draws > 0L) {
    single_draws <- .bef_plot_g_draws(engine, type, at$values, draws, scale)
    drawn_count <- length(unique(single_draws$draw))
    plot <- plot +
      ggplot2::geom_line(
        data = single_draws,
        ggplot2::aes(x = .data$value, y = .data$y, group = .data$draw),
        color = colors[["fit"]], alpha = 0.3, linewidth = 0.25
      )
  }

  # The lines: this fit, the normal distribution, the truth, other fits.
  lines <- list(data.frame(series = "Posterior mean of g", value = band$value, y = band$mean))
  values <- at$curve
  if (!is.null(normal)) {
    lines[[length(lines) + 1L]] <- data.frame(
      series = if (identical(normal$source, "REML")) "Normal (REML)" else "Normal",
      value = values * scale,
      y = if (identical(type, "cdf")) {
        stats::pnorm(values, normal$mean, normal$sd)
      } else {
        stats::dnorm(values, normal$mean, normal$sd) / scale
      }
    )
  }
  if (is.function(truth)) {
    y <- truth(values)
    if (!is.numeric(y) || length(y) != length(values)) {
      .bef_abort_invalid_args(
        "The function `truth` must return one number for each effect value that it is given.",
        arg = "truth"
      )
    }
    lines[[length(lines) + 1L]] <- data.frame(
      series = "True distribution", value = values * scale,
      y = if (identical(type, "cdf")) y else y / scale
    )
  }
  lines <- c(lines, steps, .bef_plot_g_compare(compare, type, level, scale))
  series <- do.call(rbind, lines)
  names_in_order <- unique(series$series)
  series$series <- factor(series$series, levels = names_in_order)
  line_colors <- .bef_plot_g_line_colors(names_in_order, colors)
  line_types <- .bef_plot_g_line_types(names_in_order)
  plot <- plot +
    ggplot2::geom_line(
      data = series,
      ggplot2::aes(x = .data$value, y = .data$y, color = .data$series,
                   linetype = .data$series),
      linewidth = 0.8
    ) +
    ggplot2::scale_color_manual(values = line_colors) +
    ggplot2::scale_linetype_manual(values = line_types)
  if (identical(observed, "rug")) {
    plot <- plot +
      ggplot2::geom_rug(
        data = data.frame(value = engine$theta_hat * scale),
        ggplot2::aes(x = .data$value), color = colors[["data"]], alpha = 0.7
      )
  }

  limits <- if (is.null(xlim)) .bef_plot_g_limits(engine, scale) else xlim
  plot <- plot +
    ggplot2::coord_cartesian(xlim = limits) +
    ggplot2::labs(
      x = "Effect",
      y = switch(type, density = "Density", cdf = "Cumulative probability",
                 mass = "Probability of a grid point"),
      caption = .bef_plot_g_caption(type, level, inner_level, observed, normal, truth, drawn_count,
                                    compare)
    ) +
    theme_bef()
  if (length(names_in_order) == 1L && length(fills) == 0L) {
    plot <- plot + ggplot2::theme(legend.position = "none")
  }
  plot
}

# The effect values at which the band of g is evaluated, and those at which
# smooth curves are drawn. With point masses the distribution function is
# linear between the edges of the cells, so the edges are enough.
.bef_plot_g_values <- function(engine, type) {
  grid <- engine$grid
  d <- engine$spacing
  smooth <- seq(grid[1L] - d / 2, grid[engine$L] + d / 2, length.out = 401L)
  if (engine$tau_c > 0 && !identical(type, "mass")) {
    reach <- 3 * engine$tau_c
    values <- seq(grid[1L] - reach, grid[engine$L] + reach, length.out = 401L)
    return(list(values = values, curve = values))
  }
  values <- if (identical(type, "cdf")) c(grid[1L] - d / 2, grid + d / 2) else grid
  list(values = values, curve = smooth)
}

# A histogram of a sample of effects or estimates. The bars are computed
# here over the range of the sample, so that they do not depend on what else
# the plot holds.
.bef_plot_g_histogram <- function(sample, label, outline = FALSE) {
  breaks <- seq(min(sample), max(sample), length.out = .bef_plot_bins(sample) + 1L)
  bars <- graphics::hist(sample, breaks = breaks, plot = FALSE)
  ggplot2::geom_rect(
    data = data.frame(
      xmin = utils::head(bars$breaks, -1L), xmax = bars$breaks[-1L], ymin = 0,
      ymax = bars$density, sample = label
    ),
    ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax,
                 fill = .data$sample),
    color = if (isTRUE(outline)) "grey20" else "white",
    alpha = if (isTRUE(outline)) 0 else 0.6, linewidth = 0.3
  )
}

# The empirical distribution function of a sample as the corners of its
# steps, so that it can be drawn as a line.
.bef_plot_g_ecdf <- function(sample, label) {
  sorted <- sort(sample)
  n <- length(sorted)
  data.frame(
    series = label,
    value = rep(sorted, each = 2L),
    y = as.numeric(rbind((seq_len(n) - 1) / n, seq_len(n) / n))
  )
}

# The number of bins of a histogram: the rule of Freedman and Diaconis,
# between 10 and 60.
.bef_plot_bins <- function(sample) {
  width <- 2 * stats::IQR(sample) / length(sample)^(1 / 3)
  if (!is.finite(width) || width <= 0) {
    return(30L)
  }
  as.integer(min(max(ceiling(diff(range(sample)) / width), 10), 60))
}

# Single draws of g, evenly spaced over the draws of the fit.
.bef_plot_g_draws <- function(engine, type, values, n, scale) {
  g <- exp(engine$log_g)
  g <- g / rowSums(g)
  chosen <- unique(round(seq(1, nrow(g), length.out = min(n, nrow(g)))))
  g <- g[chosen, , drop = FALSE]
  y <- switch(
    type,
    density = .bef_g_density_draws(engine, g, values) / scale,
    cdf = .bef_g_cdf_draws(engine, g, values, .bef_point_mass_reading()),
    mass = g
  )
  data.frame(
    draw = rep(seq_along(chosen), times = length(values)),
    value = rep(values * scale, each = length(chosen)),
    y = as.numeric(y)
  )
}

# Other estimates of the distribution, each as one line.
.bef_plot_g_compare <- function(compare, type, level, scale) {
  if (is.null(compare)) {
    return(list())
  }
  if (!is.list(compare) || is.data.frame(compare) || is.null(names(compare)) ||
      any(!nzchar(names(compare))) || anyDuplicated(names(compare)) > 0L) {
    .bef_abort_invalid_args(
      "`compare` must be a list whose elements have different names, each a fitted model or a data frame.",
      arg = "compare"
    )
  }
  lapply(names(compare), function(name) {
    other <- compare[[name]]
    if (inherits(other, "bef_fit")) {
      other <- .bef_prepare_fit(other)
      at <- .bef_plot_g_values(.bef_engine(other), type)
      summary <- g_summary(other, type = type, level = level, at = at$values, scale = scale)
      return(data.frame(series = name, value = summary$value, y = summary$mean))
    }
    column <- if (identical(type, "cdf")) "cdf" else "density"
    if (!is.data.frame(other) || !"x" %in% names(other) || !is.numeric(other$x) ||
        nrow(other) < 2L) {
      .bef_abort_invalid_args(
        sprintf("The element `%s` of `compare` must be a fitted model or a data frame with a column `x` and a column `density` or `cdf`.", name),
        arg = "compare"
      )
    }
    other <- other[order(other$x), , drop = FALSE]
    if (column %in% names(other)) {
      y <- other[[column]]
    } else if (identical(type, "cdf") && "density" %in% names(other)) {
      # The distribution function by the trapezoid rule, made to end at one.
      steps <- diff(other$x) * (utils::head(other$density, -1L) + other$density[-1L]) / 2
      y <- c(0, cumsum(steps))
      y <- y / y[length(y)]
    } else {
      .bef_abort_invalid_args(
        sprintf("The element `%s` of `compare` has no column `%s`. A distribution function with jumps has no density to draw; use `type = \"cdf\"`.",
                name, column),
        arg = "compare"
      )
    }
    data.frame(series = name, value = other$x * scale,
               y = if (identical(type, "cdf")) y else y / scale)
  })
}

# The color of each line: fixed for the fit and the reference curves, and
# the colors for other fits in turn.
.bef_plot_g_line_colors <- function(series, colors) {
  out <- stats::setNames(rep(NA_character_, length(series)), series)
  fixed <- c(
    "Posterior mean of g" = colors[["fit"]], "Normal (REML)" = colors[["reference"]],
    "Normal" = colors[["reference"]], "True distribution" = "#404040",
    "True effects" = "#404040", "Estimates" = "#9E9E9E"
  )
  known <- intersect(series, names(fixed))
  out[known] <- fixed[known]
  others <- unname(colors[c("other1", "other2", "other3", "other4")])
  rest <- which(is.na(out))
  out[rest] <- others[(seq_along(rest) - 1L) %% length(others) + 1L]
  out
}

# The type of each line: solid for this fit and for samples, dashed for a
# normal distribution, and other types in turn for the rest.
.bef_plot_g_line_types <- function(series) {
  out <- stats::setNames(rep(NA_character_, length(series)), series)
  fixed <- c(
    "Posterior mean of g" = "solid", "Normal (REML)" = "dashed", "Normal" = "dashed",
    "True distribution" = "solid", "True effects" = "solid", "Estimates" = "solid"
  )
  known <- intersect(series, names(fixed))
  out[known] <- fixed[known]
  others <- c("dotdash", "longdash", "twodash", "dotted")
  rest <- which(is.na(out))
  out[rest] <- others[(seq_along(rest) - 1L) %% length(others) + 1L]
  out
}

# The limits of the horizontal axis: the cells of the grid, or the part that
# holds the central 99.9% of the posterior mean of g, with a margin, where
# that part is less than 60% of the grid.
.bef_plot_g_limits <- function(engine, scale) {
  grid <- engine$grid
  d <- engine$spacing
  reach <- if (engine$tau_c > 0) 3 * engine$tau_c else d / 2
  whole <- c(grid[1L] - reach, grid[engine$L] + reach)
  g <- colMeans(exp(engine$log_g))
  cumulative <- cumsum(g / sum(g))
  first <- which(cumulative >= 0.0005)[1L]
  last <- which(cumulative >= 0.9995)[1L]
  if (is.na(last)) {
    last <- engine$L
  }
  core <- c(grid[first], grid[last])
  if (diff(core) >= 0.6 * diff(whole)) {
    return(whole * scale)
  }
  pad <- 0.15 * diff(core)
  c(max(core[1L] - pad, whole[1L]), min(core[2L] + pad, whole[2L])) * scale
}

.bef_plot_g_caption <- function(type, level, inner_level, observed, normal, truth, draws,
                                compare) {
  what <- switch(type, density = "density", cdf = "distribution function",
                 mass = "probabilities")
  .bef_caption(
    sprintf("Solid line: posterior mean of the %s of g. Band: pointwise %s posterior interval%s.",
            what, .bef_percent(level),
            if (is.null(inner_level)) "" else sprintf(", with the %s interval inside",
                                                      .bef_percent(inner_level))),
    if (draws > 0L) sprintf("Thin lines: %d posterior draws of g.", draws),
    if (!is.null(normal)) {
      if (identical(normal$source, "REML")) {
        "Dashed line: normal distribution with its mean and standard deviation estimated by REML; the uncertainty of those estimates is not shown."
      } else {
        "Dashed line: the normal distribution that was given."
      }
    },
    if (identical(observed, "histogram") && is.numeric(truth)) {
      if (identical(type, "cdf")) {
        "Dark steps: the true effects. Gray steps: the estimates, which spread more widely than the effects."
      } else {
        "Gray bars: the true effects. Outlined bars: the estimates, which spread more widely than the effects."
      }
    } else if (identical(observed, "histogram")) {
      "Gray: the estimates, which spread more widely than the effects, so that their distribution is not an estimate of g."
    } else if (is.numeric(truth)) {
      if (identical(type, "cdf")) "Dark steps: the true effects." else "Gray bars: the true effects."
    },
    if (identical(observed, "rug")) "Marks on the axis: the estimates.",
    if (is.function(truth)) "Dark line: the true distribution.",
    if (!is.null(compare)) "Other lines: the distributions named in the legend, without bands."
  )
}
