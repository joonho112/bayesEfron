#' Plot a grid and its spline basis with the estimates
#'
#' @description
#' The `plot()` method for a grid created by [make_efron_grid()] shows
#' where the grid lies in relation to the data. The spline functions are
#' plotted over the range of the grid, with the knots and the grid points,
#' and the estimates are shown below them on the same axis, with 2
#' standard errors on each side. The plot shows whether estimates lie
#' outside the range of the grid, whether the range is much wider than the
#' estimates, and whether there are knots in regions without estimates.
#'
#' @details
#' Up to 100 estimates are shown one by one, in increasing order, each as
#' a point with a line from 2 standard errors below it to 2 standard
#' errors above it. More than 100 estimates are shown as a histogram.
#'
#' A grid created by [make_efron_grid()] contains the estimates from which
#' it was built, but the grid stored in a fit does not. To plot the grid
#' of a fit, use `plot(fit, type = "grid")`, which adds the data of the
#' fit, or supply the data as `theta_hat` and `sigma`. If there are no
#' estimates, only the spline functions are shown.
#'
#' The spline functions are the natural cubic spline basis computed by
#' [splines::ns()]. The boundary knots are at the limits of the grid,
#' which the caption of the plot calls the ends of the grid. The interior
#' knots are at quantiles of the grid points, so they are equally spaced.
#'
#' @param x A `bef_grid` object created by [make_efron_grid()].
#' @param theta_hat,sigma Estimates and their standard errors, as numeric
#'   vectors. As in [make_efron_grid()], `theta_hat` can also be a data
#'   frame or a `bef_data` object. If `NULL` (the default), the estimates
#'   stored in the grid are used.
#' @param scale Multiplier for the values on the effect axis, a positive
#'   number. For example, `scale = 100` shows differences in proportions
#'   as percentage points.
#' @param ... These dots are for future extensions and must be empty.
#'
#' @return A `ggplot` object.
#'
#' @seealso [make_efron_grid()], [plot.bef_grid_comparison()],
#'   `vignette("choosing-a-grid")`
#'
#' @examples
#' grid <- make_efron_grid(krw_discrimination, grid_method = "paper_realdata")
#' plot(grid, scale = 100)
#'
#' # The grid of a fit, with the data of the fit.
#' plot(twin_towers_fit, type = "grid")
#'
#' @export
plot.bef_grid <- function(x, theta_hat = NULL, sigma = NULL, scale = 1, ...) {
  .bef_check_dots_empty(list(...), "plot")
  .bef_require_inherits(x, "bef_grid", "x", "bef_invalid_args")
  scale <- .bef_validate_scale(scale)
  data <- x$data
  if (!is.null(theta_hat)) {
    given <- .bef_resolve_data_args(theta_hat, sigma)
    data <- list(
      theta_hat = unname(.bef_assert_numeric_vector(given$theta_hat, "theta_hat", min_len = 1L)),
      sigma = unname(.bef_assert_numeric_vector(given$sigma, "sigma", min_len = 1L))
    )
    if (length(data$sigma) != length(data$theta_hat) || any(data$sigma <= 0)) {
      .bef_abort_invalid_args(
        "`sigma` must be strictly positive and the same length as `theta_hat`.", arg = "sigma"
      )
    }
  } else if (!is.null(sigma)) {
    .bef_abort_invalid_args("`sigma` needs the estimates in `theta_hat`.", arg = "sigma")
  }
  .bef_plot_grid(x, data, scale)
}

# The number of estimates up to which each is drawn with its own line.
.bef_plot_grid_rule <- function() {
  list(marks = 100L)
}

.bef_plot_grid <- function(grid, data, scale) {
  colors <- .bef_plot_colors()
  points <- grid$grid * scale
  top <- "Spline functions"
  K <- length(data$theta_hat)
  marks <- K > 0L && K <= .bef_plot_grid_rule()$marks
  bottom <- if (K == 0L) {
    NULL
  } else if (marks) {
    "Estimates, with two standard errors on each side"
  } else {
    "Estimates"
  }
  panels <- c(top, bottom)
  as_panel <- function(name) factor(name, levels = panels)

  basis <- data.frame(
    panel = as_panel(top),
    value = rep(points, times = ncol(grid$B)),
    y = as.numeric(grid$B),
    fn = rep(seq_len(ncol(grid$B)), each = length(points))
  )
  knots <- rbind(
    data.frame(value = grid$knots$boundary * scale, kind = "boundary"),
    if (length(grid$knots$interior) > 0L) {
      data.frame(value = grid$knots$interior * scale, kind = "interior")
    }
  )
  ticks <- data.frame(panel = as_panel(top), value = points)

  plot <- ggplot2::ggplot() +
    ggplot2::geom_vline(data = knots,
                        ggplot2::aes(xintercept = .data$value, linetype = .data$kind,
                                     color = .data$kind)) +
    ggplot2::scale_linetype_manual(values = c(boundary = "solid", interior = "dashed"),
                                   guide = "none") +
    ggplot2::scale_color_manual(values = c(boundary = "grey25", interior = "grey60"),
                                guide = "none") +
    ggplot2::geom_line(data = basis,
                       ggplot2::aes(x = .data$value, y = .data$y, group = .data$fn),
                       color = colors[["fit"]], linewidth = 0.5) +
    ggplot2::geom_rug(data = ticks, ggplot2::aes(x = .data$value), sides = "b",
                      color = "grey35", linewidth = 0.2,
                      length = ggplot2::unit(0.025, "npc"))

  if (marks) {
    estimates <- data.frame(
      panel = as_panel(bottom),
      value = sort(data$theta_hat) * scale,
      lower = (data$theta_hat - 2 * data$sigma)[order(data$theta_hat)] * scale,
      upper = (data$theta_hat + 2 * data$sigma)[order(data$theta_hat)] * scale,
      y = seq_len(K)
    )
    plot <- plot +
      ggplot2::geom_segment(data = estimates,
                            ggplot2::aes(x = .data$lower, xend = .data$upper, y = .data$y,
                                         yend = .data$y),
                            color = colors[["data"]], linewidth = 0.3) +
      ggplot2::geom_point(data = estimates, ggplot2::aes(x = .data$value, y = .data$y),
                          color = "grey25", size = if (K > 40L) 0.7 else 1.4)
  } else if (K > 0L) {
    values <- data$theta_hat * scale
    bins <- graphics::hist(values, breaks = .bef_plot_bins(values), plot = FALSE)
    bars <- data.frame(
      panel = as_panel(bottom),
      xmin = bins$breaks[-length(bins$breaks)], xmax = bins$breaks[-1L],
      ymin = 0, ymax = bins$counts
    )
    plot <- plot +
      ggplot2::geom_rect(data = bars,
                         ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin,
                                      ymax = .data$ymax),
                         fill = "grey75", color = "white", linewidth = 0.2)
  }
  if (length(panels) > 1L) {
    plot <- plot + ggplot2::facet_wrap(ggplot2::vars(.data$panel), ncol = 1L, scales = "free_y")
  }
  plot +
    ggplot2::labs(
      x = "Effect", y = NULL,
      caption = .bef_plot_grid_caption(grid, data, K, marks)
    ) +
    theme_bef() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   axis.text.y = if (marks) ggplot2::element_blank())
}

.bef_plot_grid_caption <- function(grid, data, K, marks) {
  outside <- if (K > 0L) {
    sum(data$theta_hat < grid$support[1L] | data$theta_hat > grid$support[2L])
  }
  .bef_caption(
    sprintf("Curves: the %d basis functions of the natural cubic spline.", grid$M),
    sprintf("Ticks: the %d grid points.", grid$L),
    "Solid vertical lines: the ends of the grid, which are the boundary knots.",
    if (length(grid$knots$interior) > 0L) "Dashed vertical lines: the interior knots.",
    if (K == 0L) {
      "The grid holds no estimates; give them in `theta_hat` and `sigma` to see them."
    } else if (marks) {
      sprintf("Points and lines: the %d estimates, in increasing order, and two standard errors on each side.",
              K)
    } else {
      sprintf("Bars: the %s estimates.", format(K, big.mark = ","))
    },
    if (!is.null(outside) && outside > 0L) {
      sprintf("%d estimate%s outside the grid.", outside,
              if (outside == 1L) " lies" else "s lie")
    }
  )
}

#' Plot the sensitivity of the results to the grid
#'
#' @description
#' The `plot()` method for a comparison created by [compare_efron_grids()]
#' shows one row per specification, so that you can see how the results
#' change with the grid. Each row
#' shows the estimated density of \eqn{g} with its pointwise posterior
#' band, the posterior means of the mean and of the standard deviation of
#' \eqn{g}, and the posterior means and intervals of the effects of the
#' selected sites. Read down a column to see what changes when the range,
#' the number of grid points or the number of spline functions is varied.
#'
#' @details
#' The specifications are shown in the order in which they were fitted,
#' with the first at the top. They are labeled by the settings that differ
#' between them; in these labels, "rule" means the grid method. A
#' specification whose fit failed still has a row, which is labeled as
#' failed and is otherwise empty.
#'
#' All densities are plotted on the same vertical scale, so their heights
#' can be compared between rows, but the scale itself is not shown. To
#' avoid exaggerating small differences between specifications, the axis
#' for the mean of \eqn{g} is at least as wide as the standard deviation
#' of \eqn{g}, and the axis for the standard deviation starts at zero.
#'
#' The sites are those selected with the `sites` argument of
#' [compare_efron_grids()]. By default these are the sites with the
#' smallest and the largest estimates, whose effects depend most on the
#' limits of the grid.
#'
#' @param x A `bef_grid_comparison` object created by
#'   [compare_efron_grids()].
#' @inheritParams plot.bef_grid
#'
#' @return A `ggplot` object.
#'
#' @seealso [compare_efron_grids()], [plot.bef_grid()],
#'   `vignette("choosing-a-grid")`
#'
#' @examples
#' # Three ranges for the 19 teacher expectancy studies, stored with the
#' # package. The examples of compare_efron_grids() show the call that made
#' # the comparison.
#' comparison <- readRDS(system.file("examples", "raudenbush_ranges.rds",
#'                                   package = "bayesEfron"))
#' comparison[, c("expansion", "lower", "upper", "sd_g", "site_4")]
#'
#' # As the range widens, the mean of study 4 moves toward the others.
#' plot(comparison)
#'
#' @export
plot.bef_grid_comparison <- function(x, scale = 1, ...) {
  .bef_check_dots_empty(list(...), "plot")
  .bef_require_inherits(x, "bef_grid_comparison", "x", "bef_invalid_args")
  scale <- .bef_validate_scale(scale)
  g <- attr(x, "g", exact = TRUE)
  kept <- attr(x, "sites", exact = TRUE)
  if (is.null(g) || !"spec" %in% names(x)) {
    .bef_abort_invalid_args(
      paste("This comparison does not hold the estimates of g of its fits, which version 0.4",
            "keeps. Run `compare_efron_grids()` again."),
      arg = "x"
    )
  }
  specs <- as.data.frame(x)[order(x$spec), , drop = FALSE]
  n <- nrow(specs)
  specs$failed <- vapply(as.character(specs$spec), function(id) is.null(g[[id]]), logical(1L))
  if (all(specs$failed)) {
    .bef_abort_invalid_args("Every fit of this comparison failed, so there is nothing to draw.",
                            arg = "x")
  }
  # The first specification is the top row.
  specs$position <- rev(seq_len(n))
  specs$row_label <- paste0(.bef_comparison_labels(specs),
                            ifelse(specs$failed, " (fit failed)", ""))
  position_of <- stats::setNames(specs$position, as.character(specs$spec))

  site_labels <- unique(kept$label)
  panels <- c("Density of g", "Mean of g", "SD of g", sprintf("Site %s", site_labels))
  as_panel <- function(name) factor(name, levels = panels)

  # The densities on one vertical scale: the highest band fills 85% of a row.
  drawn <- specs[!specs$failed, , drop = FALSE]
  top <- max(vapply(as.character(drawn$spec), function(id) max(g[[id]]$upper), numeric(1L)))
  ridges <- do.call(rbind, lapply(as.character(drawn$spec), function(id) {
    curve <- g[[id]]
    base <- position_of[[id]]
    data.frame(panel = as_panel("Density of g"), spec = id, value = curve$value * scale,
               y = base + 0.85 * curve$mean / top, ymin = base + 0.85 * curve$lower / top,
               ymax = base + 0.85 * curve$upper / top, base = base)
  }))
  summaries <- rbind(
    data.frame(panel = as_panel("Mean of g"), position = drawn$position,
               value = drawn$mean_g * scale, lower = NA_real_, upper = NA_real_),
    data.frame(panel = as_panel("SD of g"), position = drawn$position,
               value = drawn$sd_g * scale, lower = NA_real_, upper = NA_real_),
    data.frame(panel = as_panel(sprintf("Site %s", kept$label)),
               position = unname(position_of[as.character(kept$spec)]),
               value = kept$mean * scale, lower = kept$lower * scale,
               upper = kept$upper * scale)
  )
  intervals <- summaries[!is.na(summaries$lower), , drop = FALSE]
  # The axis of the mean of g is at least as wide as the standard deviation
  # of g, and that of the standard deviation starts at zero, so that small
  # differences between the specifications are not magnified.
  spread <- mean(drawn$sd_g) * scale
  room <- rbind(
    data.frame(panel = as_panel("Mean of g"), position = 1,
               value = mean(range(drawn$mean_g)) * scale + c(-0.5, 0.5) * spread),
    data.frame(panel = as_panel("SD of g"), position = 1,
               value = c(0, 1.05 * max(drawn$sd_g) * scale))
  )
  colors <- .bef_plot_colors()
  level <- attr(g[[as.character(drawn$spec[[1L]])]], "level", exact = TRUE)
  if (is.null(level)) {
    level <- 0.9
  }

  plot <- ggplot2::ggplot() +
    ggplot2::geom_blank(data = room, ggplot2::aes(x = .data$value, y = .data$position)) +
    ggplot2::geom_ribbon(data = ridges,
                         ggplot2::aes(x = .data$value, ymin = .data$ymin, ymax = .data$ymax,
                                      group = .data$spec),
                         fill = colors[["fit"]], alpha = 0.25) +
    ggplot2::geom_line(data = ridges,
                       ggplot2::aes(x = .data$value, y = .data$y, group = .data$spec),
                       color = colors[["fit"]], linewidth = 0.5) +
    ggplot2::geom_segment(data = intervals,
                          ggplot2::aes(x = .data$lower, xend = .data$upper, y = .data$position,
                                       yend = .data$position),
                          color = colors[["fit"]], linewidth = 0.4) +
    ggplot2::geom_point(data = summaries, ggplot2::aes(x = .data$value, y = .data$position),
                        color = colors[["fit"]], size = 1.8) +
    ggplot2::facet_grid(cols = ggplot2::vars(.data$panel), scales = "free_x") +
    ggplot2::scale_x_continuous(n.breaks = 4, expand = ggplot2::expansion(mult = 0.08),
                                guide = ggplot2::guide_axis(check.overlap = TRUE)) +
    ggplot2::scale_y_continuous(breaks = specs$position, labels = specs$row_label,
                                limits = c(0.75, n + 1)) +
    ggplot2::labs(
      x = "Effect", y = NULL,
      caption = .bef_caption(
        "Each row is a specification.",
        sprintf("Curves: posterior mean of the density of g with a pointwise %s posterior interval, all on one vertical scale.",
                .bef_percent(level)),
        "Points: posterior means of the mean and the standard deviation of g and of the effects of the sites.",
        sprintf("Lines: %s posterior intervals of those effects.", .bef_percent(level)),
        if (any(specs$failed)) "A row without marks is a specification whose fit failed."
      )
    ) +
    theme_bef() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = "grey92"),
                   panel.spacing.x = ggplot2::unit(1, "lines"))
  if (utils::packageVersion("ggplot2") >= "4.0.0") {
    # The panel of the densities is given twice the width of the others.
    plot <- plot +
      ggplot2::theme(panel.widths = ggplot2::unit(c(2, rep(1, length(panels) - 1L)), "null"))
  }
  plot
}

# The specifications of a comparison in words: the settings that differ
# between them, or the numbers of grid points and spline functions when
# none does.
.bef_comparison_labels <- function(specs) {
  varies <- function(column) length(unique(specs[[column]])) > 1L
  given_range <- specs$grid_method == "support"
  parts <- list(
    if (varies("L")) sprintf("L = %d", specs$L),
    if (varies("M")) sprintf("M = %d", specs$M),
    if (varies("grid_method")) sprintf("rule \"%s\"", specs$grid_method),
    if (varies("expansion")) {
      ifelse(is.na(specs$expansion), "", sprintf("expansion %s", format(specs$expansion)))
    },
    if (any(given_range) && (varies("lower") || varies("upper"))) {
      ifelse(given_range & !is.na(specs$lower),
             sprintf("%s to %s", format(signif(specs$lower, 3L)), format(signif(specs$upper, 3L))),
             "")
    },
    if (varies("component")) {
      ifelse(specs$component == "normal", "normal components", "point masses")
    }
  )
  parts <- parts[!vapply(parts, is.null, logical(1L))]
  if (length(parts) == 0L) {
    return(sprintf("L = %d, M = %d", specs$L, specs$M))
  }
  labels <- do.call(paste, c(parts, list(sep = ", ")))
  labels <- gsub("(, )+", ", ", labels)
  labels <- gsub("^, |, $", "", labels)
  ifelse(nzchar(labels), labels, sprintf("specification %d", specs$spec))
}
