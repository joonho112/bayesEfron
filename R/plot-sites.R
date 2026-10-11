#' Caterpillar (forest) plot of the site effects
#'
#' @description
#' `plot_sites()` makes a caterpillar plot, or forest plot, of the site
#' effects. Each row is a labeled site with the posterior mean of its
#' effect, an outer and an inner posterior interval, and its estimate.
#' When there are many sites, the plot can be limited to selected sites or
#' to the sites with the lowest and the highest posterior means.
#'
#' @details
#' The point in each row is the posterior mean, the thick line is the
#' inner interval and the thin line is the outer interval. The open circle
#' is the estimate, so its distance from the point is the shrinkage of the
#' site, i.e., how far the model moves the estimate ([plot_shrinkage()]
#' shows this for all sites together). The dashed line is the reference
#' line selected with `ref`.
#'
#' The plot shows the values returned by [site_summary()], so its intervals
#' are the intervals of [confint.bef_fit_re()] at the same level.
#'
#' If more than `max_sites` sites are to be plotted, they are shown as
#' thin lines without labels, inner intervals or estimates, and a message
#' explains how to select fewer sites.
#'
#' @inheritParams plot_g
#' @param sites Sites to plot, specified as a numeric vector of positions
#'   in the data or a character vector of site labels. If `NULL` (the
#'   default), all sites are plotted.
#' @param n_extreme An optional integer. If specified, only the `n_extreme`
#'   sites with the lowest and the `n_extreme` sites with the highest
#'   posterior means are plotted, in two panels. Cannot be used together
#'   with `sites`.
#' @param sort_by Variable by which the sites are ordered: `"mean"` (the
#'   default) for the posterior mean, `"sigma"` for the standard error,
#'   `"observed"` for the estimate, or `"prob"` for the probability that
#'   the effect is positive. `"none"` keeps the order of the data.
#' @param level Posterior probability of the outer interval, a number
#'   between 0 and 1.
#' @param inner_level Posterior probability of the inner interval, a number
#'   between 0 and `level`. The inner interval is shown as a thicker line.
#'   Use `NULL` to omit it.
#' @param highlight Sites to show in a different color, specified in the
#'   same way as `sites`. If `NULL` (the default), no site is highlighted.
#' @param show_observed If `TRUE` (the default), the estimates are shown as
#'   open circles.
#' @param ref Where to place the dashed reference line: at the posterior
#'   mean of the mean of \eqn{g} (`"mean_g"`, the default) or at zero
#'   (`"zero"`). Use `"none"` for no line.
#' @param truth An optional numeric vector with the true effect of every
#'   site, for simulated data. The true effects are shown as crosses.
#' @param max_sites Maximum number of sites that are plotted with labels.
#'   See Details for how a larger number of sites is shown.
#' @param scale Positive number by which the effects are multiplied in the
#'   plot, e.g., `scale = 100` to show differences in proportions as
#'   percentage points.
#'
#' @return A `ggplot` object.
#'
#' @seealso [site_summary()], [plot_shrinkage()], [plot_exceedance()],
#'   `vignette("plots")`
#'
#' @examples
#' # The 10 firms with the smallest and the 10 with the largest gaps, in
#' # percentage points.
#' plot_sites(krw_fit, n_extreme = 10, scale = 100)
#'
#' # All 19 studies, ordered by standard error.
#' plot_sites(raudenbush_fit, sort_by = "sigma", ref = "zero")
#'
#' # 1,500 sites: lines without labels.
#' plot_sites(twin_towers_fit)
#'
#' @export
plot_sites <- function(x,
                       sites = NULL,
                       n_extreme = NULL,
                       sort_by = c("mean", "sigma", "observed", "prob", "none"),
                       level = 0.9,
                       inner_level = 0.5,
                       highlight = NULL,
                       show_observed = TRUE,
                       ref = c("mean_g", "zero", "none"),
                       truth = NULL,
                       max_sites = 100,
                       scale = 1,
                       ...) {
  .bef_check_dots_empty(list(...), "plot_sites")
  x <- .bef_plot_fit(x)
  sort_by <- .bef_validate_method_choice(
    sort_by, c("mean", "sigma", "observed", "prob", "none"), arg = "sort_by"
  )
  ref <- .bef_validate_method_choice(ref, c("mean_g", "zero", "none"), arg = "ref")
  level <- .bef_plot_level(level)
  inner_level <- .bef_plot_level(inner_level, "inner_level", null_ok = TRUE)
  if (!is.null(inner_level) && inner_level >= level) {
    .bef_abort_invalid_args("`inner_level` must be smaller than `level`.", arg = "inner_level")
  }
  show_observed <- .bef_plot_flag(show_observed, "show_observed")
  max_sites <- .bef_plot_count(max_sites, "max_sites", lower = 1L)
  scale <- .bef_validate_scale(scale)
  if (!is.null(sites) && !is.null(n_extreme)) {
    .bef_abort_invalid_args("Give `sites` or `n_extreme`, not both.", arg = "n_extreme")
  }

  table <- site_summary(x, level = level, inner_level = inner_level, threshold = 0,
                        scale = scale)
  K <- nrow(table)
  truth <- .bef_plot_truth(truth, K)
  if (!is.null(truth)) {
    table$truth <- truth * scale
  }
  table$highlight <- FALSE
  marked <- .bef_plot_sites(highlight, table$label, arg = "highlight")
  table$highlight[marked] <- TRUE
  table$panel <- ""

  chosen <- .bef_plot_sites(sites, table$label)
  if (!is.null(n_extreme)) {
    n_extreme <- .bef_plot_count(n_extreme, "n_extreme", lower = 1L)
    if (2L * n_extreme > K) {
      .bef_abort_invalid_args(
        sprintf("`n_extreme` must be at most %d, half the number of sites.", K %/% 2L),
        arg = "n_extreme"
      )
    }
    by_mean <- order(table$mean)
    lowest <- by_mean[seq_len(n_extreme)]
    highest <- by_mean[seq.int(K - n_extreme + 1L, K)]
    table$panel[lowest] <- sprintf("Lowest %d", n_extreme)
    table$panel[highest] <- sprintf("Highest %d", n_extreme)
    chosen <- c(lowest, highest)
  }
  if (!is.null(chosen)) {
    table <- table[chosen, , drop = FALSE]
  }
  labeled <- nrow(table) <= max_sites
  if (!labeled) {
    .bef_inform_many_sites(nrow(table), max_sites)
  }

  # The sites are placed from the top down: the largest posterior mean,
  # estimate or probability first, the smallest standard error first, or
  # the order of the data (the order of `sites`, when sites were chosen).
  key <- switch(
    sort_by,
    mean = table$mean, sigma = -table$se, observed = table$estimate,
    prob = table$prob_above,
    none = if (is.null(sites)) -table$site else -seq_len(nrow(table))
  )
  in_upper_panel <- startsWith(table$panel, "Highest")
  table <- table[order(in_upper_panel, key, -table$site), , drop = FALSE]
  table$position <- seq_len(nrow(table))
  panels <- unique(table$panel)
  if (length(panels) > 1L) {
    table$panel <- factor(table$panel, levels = rev(panels))
  }

  colors <- .bef_plot_colors()
  table$color <- ifelse(table$highlight, "Highlighted", "Site")
  plot <- ggplot2::ggplot(table, ggplot2::aes(y = .data$position))
  reference <- switch(
    ref,
    mean_g = mean(x$posterior$mean_g) * scale, zero = 0, none = NULL
  )
  reference_line <- if (!is.null(reference)) {
    ggplot2::geom_vline(xintercept = reference, linetype = "dashed", color = "grey35")
  }
  intervals <- ggplot2::geom_segment(
    ggplot2::aes(x = .data$lower, xend = .data$upper, yend = .data$position,
                 color = .data$color),
    linewidth = if (labeled) 0.4 else 0.2, alpha = if (labeled) 1 else 0.6
  )
  # With labels the reference line lies behind the intervals; among the
  # many lines of a plot without labels it would not be seen there.
  plot <- if (labeled) plot + reference_line + intervals else plot + intervals + reference_line
  if (!is.null(inner_level) && labeled) {
    plot <- plot +
      ggplot2::geom_segment(
        ggplot2::aes(x = .data$inner_lower, xend = .data$inner_upper, yend = .data$position,
                     color = .data$color),
        linewidth = 1.2
      )
  }
  if (show_observed && labeled) {
    plot <- plot +
      ggplot2::geom_point(ggplot2::aes(x = .data$estimate), shape = 1, size = 1.8,
                          color = colors[["data"]])
  }
  if (!is.null(truth)) {
    plot <- plot +
      ggplot2::geom_point(ggplot2::aes(x = .data$truth), shape = 4,
                          size = if (labeled) 1.8 else 0.6, color = "#404040")
  }
  plot <- plot +
    ggplot2::geom_point(ggplot2::aes(x = .data$mean, color = .data$color),
                        size = if (labeled) 1.8 else 0.4) +
    ggplot2::scale_color_manual(
      values = c(Site = colors[["fit"]], Highlighted = colors[["highlight"]]), guide = "none"
    )

  if (labeled) {
    plot <- plot +
      ggplot2::scale_y_continuous(breaks = table$position,
                                  labels = .bef_short_labels(table$label),
                                  expand = ggplot2::expansion(add = 0.6))
  } else {
    plot <- plot +
      ggplot2::scale_y_continuous(breaks = NULL, expand = ggplot2::expansion(mult = 0.01))
  }
  if (length(panels) > 1L) {
    plot <- plot + ggplot2::facet_wrap(ggplot2::vars(.data$panel), ncol = 1L, scales = "free_y")
  }
  rows <- nrow(table) %/% length(panels)
  plot +
    ggplot2::labs(
      x = "Effect", y = NULL,
      caption = .bef_plot_sites_caption(level, inner_level, show_observed, ref, truth, sort_by,
                                        labeled, any(table$highlight))
    ) +
    theme_bef() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = .bef_plot_label_size(rows))
    )
}

# The size of the labels of the sites, smaller when a panel holds many.
.bef_plot_label_size <- function(rows) {
  if (rows <= 30L) 9 else if (rows <= 60L) 7 else 5.5
}

.bef_plot_sites_caption <- function(level, inner_level, show_observed, ref, truth, sort_by,
                                    labeled, highlighted) {
  order <- switch(
    sort_by,
    mean = "by posterior mean, the largest at the top",
    sigma = "by standard error, the smallest at the top",
    observed = "by estimate, the largest at the top",
    prob = "by the probability of a positive effect, the largest at the top",
    none = NULL
  )
  .bef_caption(
    if (labeled) {
      sprintf("Points: posterior means. %s",
              if (is.null(inner_level)) {
                sprintf("Lines: %s posterior intervals.", .bef_percent(level))
              } else {
                sprintf("Thick and thin lines: %s and %s posterior intervals.",
                        .bef_percent(inner_level), .bef_percent(level))
              })
    } else {
      sprintf("Each line is the %s posterior interval of one site, with its posterior mean.",
              .bef_percent(level))
    },
    if (show_observed && labeled) "Open circles: the estimates.",
    if (!is.null(truth)) "Crosses: the true effects.",
    if (identical(ref, "mean_g")) "Dashed line: posterior mean of the mean of g.",
    if (identical(ref, "zero")) "Dashed line: zero.",
    if (highlighted) "Orange: the highlighted sites.",
    if (!is.null(order)) sprintf("Sites are ordered %s.", order)
  )
}
