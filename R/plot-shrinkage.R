#' Shrinkage plot of posterior means against estimates
#'
#' @description
#' The shrinkage of a site is the difference between its estimate and the
#' posterior mean of its effect. `plot_shrinkage()` shows it by plotting
#' the posterior means against the estimates. An imprecise estimate far
#' from the values supported by the estimated distribution of effects is
#' shrunk strongly toward them, while a precise estimate is hardly changed.
#'
#' @details
#' In the scatterplot (`style = "scatter"`), each site is a point with its
#' estimate on the horizontal axis and its posterior mean on the vertical
#' axis. The two axes have the same scale, and the color of a point shows
#' the standard error of its estimate. A point on the dashed diagonal has
#' not been moved, and the vertical distance of a point from the diagonal
#' is the shrinkage of the site.
#'
#' In the arrow plot (`style = "arrows"`), each site is an arrow from its
#' estimate to its posterior mean, plotted at the height of its standard
#' error, so color is not needed to show the standard error. Sites are
#' labeled if there are 30 or fewer. If there are more than 100 sites, the
#' arrows have no heads and a point marks the posterior mean.
#'
#' If `truth` is specified, the color shows the true effects unless
#' `color_by` is set, and `style = "arrows"` also marks them with crosses
#' if there are at most 100 sites. In both styles, the points and arrows
#' are semi-transparent if there are more than 100 sites, so that regions
#' where they overlap appear darker.
#'
#' The panel added by `normal` shows the same plot for a normal
#' random-effects model: the mean of each effect given its estimate. These
#' are plug-in means, i.e., the mean and the standard deviation of the
#' model are fixed at their estimates. The normal model shrinks every
#' estimate toward a single mean, by a fraction that depends only on the
#' standard error. This package's model instead shrinks an estimate toward
#' the values supported by \eqn{g}.
#'
#' @inheritParams plot_sites
#' @inheritParams plot_g
#' @param style One of `"scatter"` (the default) for a scatterplot of
#'   posterior means against estimates, or `"arrows"` for an arrow from the
#'   estimate to the posterior mean of each site.
#' @param color_by Variable shown by the color of a site: `"sigma"` (the
#'   default) for the standard error, `"truth"` for the true effect
#'   specified in `truth`, or `"none"`. If `truth` is specified and
#'   `color_by` is not, the color shows the true effects.
#' @param normal If `TRUE`, a panel is added for the normal random-effects
#'   model fitted to the same estimates by restricted maximum likelihood.
#'   Can also be a numeric vector `c(mean = , sd = )`, which adds a panel
#'   for the normal model with that mean and standard deviation.
#' @param highlight Sites to show in a different color and to label,
#'   specified as a numeric vector of positions in the data or a character
#'   vector of site labels. If `NULL` (the default), no site is
#'   highlighted.
#' @param truth An optional numeric vector with the true effect of every
#'   site, for simulated data.
#'
#' @return A `ggplot` object.
#'
#' @seealso [site_summary()], [plot_sites()], `vignette("plots")`
#'
#' @examples
#' # The callback gaps in percentage points: imprecise estimates far from
#' # the bulk are shrunk the most.
#' plot_shrinkage(krw_fit, scale = 100)
#'
#' # The arrows style, for 19 studies.
#' plot_shrinkage(raudenbush_fit, style = "arrows")
#'
#' # Simulated sites with two groups of true effects, with a panel for the
#' # normal model.
#' plot_shrinkage(twin_towers_fit, normal = TRUE, truth = twin_towers$theta_true)
#'
#' @export
plot_shrinkage <- function(x,
                           style = c("scatter", "arrows"),
                           color_by = c("sigma", "truth", "none"),
                           normal = FALSE,
                           highlight = NULL,
                           truth = NULL,
                           scale = 1,
                           ...) {
  .bef_check_dots_empty(list(...), "plot_shrinkage")
  x <- .bef_plot_fit(x)
  style <- .bef_validate_method_choice(style, c("scatter", "arrows"), arg = "style")
  color_given <- !missing(color_by)
  color_by <- .bef_validate_method_choice(color_by, c("sigma", "truth", "none"),
                                          arg = "color_by")
  scale <- .bef_validate_scale(scale)

  table <- site_summary(x, scale = scale)
  K <- nrow(table)
  truth <- .bef_plot_truth(truth, K)
  if (!color_given && !is.null(truth)) {
    color_by <- "truth"
  }
  if (identical(color_by, "truth") && is.null(truth)) {
    .bef_abort_invalid_args("`color_by = \"truth\"` needs the true effects in `truth`.",
                            arg = "color_by")
  }
  # In the arrows style the standard error is the vertical axis.
  if (identical(style, "arrows") && identical(color_by, "sigma")) {
    color_by <- "none"
  }
  normal <- .bef_plot_normal(normal, table$estimate / scale, table$se / scale)

  sites <- table[, c("site", "label", "estimate", "se", "mean")]
  sites$truth <- if (is.null(truth)) NA_real_ else truth * scale
  sites$highlight <- FALSE
  sites$highlight[.bef_plot_sites(highlight, table$label, arg = "highlight")] <- TRUE
  sites$model <- "This fit"
  centers <- data.frame(model = "This fit", center = mean(x$posterior$mean_g) * scale)
  if (!is.null(normal)) {
    name <- if (identical(normal$source, "REML")) "Normal (REML)" else "Normal"
    other <- sites
    other$mean <- .bef_normal_site_means(table$estimate / scale, table$se / scale, normal) * scale
    other$model <- name
    sites <- rbind(sites, other)
    centers <- rbind(centers, data.frame(model = name, center = normal$mean * scale))
    sites$model <- factor(sites$model, levels = c("This fit", name))
    centers$model <- factor(centers$model, levels = c("This fit", name))
  }
  sites$color <- switch(color_by, sigma = sites$se, truth = sites$truth, none = NA_real_)

  plot <- if (identical(style, "scatter")) {
    .bef_plot_shrinkage_scatter(sites, K, color_by)
  } else {
    .bef_plot_shrinkage_arrows(sites, centers, K, color_by, !is.null(truth))
  }
  if (!identical(color_by, "none")) {
    plot <- plot +
      ggplot2::scale_color_viridis_c(
        name = if (identical(color_by, "sigma")) "Standard error" else "True effect",
        end = 0.9
      )
  }
  if (!is.null(normal)) {
    plot <- plot + ggplot2::facet_wrap(ggplot2::vars(.data$model), nrow = 1L)
  }
  plot +
    ggplot2::labs(
      caption = .bef_plot_shrinkage_caption(style, K, color_by, normal, !is.null(truth),
                                            any(sites$highlight))
    ) +
    theme_bef() +
    ggplot2::theme(legend.title = ggplot2::element_text(size = ggplot2::rel(0.9)),
                   legend.key.width = ggplot2::unit(1.6, "lines"))
}

# The number of sites up to which arrows have heads and estimates and true
# effects their own marks, and the number up to which sites are labeled.
.bef_plot_shrinkage_rule <- function() {
  list(marks = 100L, labels = 30L)
}

# The transparency of the points of `K` sites.
.bef_plot_alpha <- function(K) {
  max(0.15, min(1, 100 / K))
}

# The mean of each effect given its estimate under a normal model of the
# effects with a known mean and standard deviation.
.bef_normal_site_means <- function(theta_hat, sigma, normal) {
  share <- normal$sd^2 / (normal$sd^2 + sigma^2)
  normal$mean + share * (theta_hat - normal$mean)
}

.bef_plot_shrinkage_scatter <- function(sites, K, color_by) {
  colors <- .bef_plot_colors()
  limits <- range(sites$estimate, sites$mean)
  marked <- sites[sites$highlight, , drop = FALSE]
  plot <- ggplot2::ggplot(sites, ggplot2::aes(x = .data$estimate, y = .data$mean)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50")
  plot <- if (identical(color_by, "none")) {
    plot + ggplot2::geom_point(color = colors[["fit"]], alpha = .bef_plot_alpha(K),
                               size = if (K > 100L) 0.9 else 1.8)
  } else {
    plot + ggplot2::geom_point(ggplot2::aes(color = .data$color), alpha = .bef_plot_alpha(K),
                               size = if (K > 100L) 0.9 else 1.8)
  }
  if (nrow(marked) > 0L) {
    plot <- plot +
      ggplot2::geom_point(data = marked, color = colors[["highlight"]], shape = 21,
                          fill = colors[["highlight"]], size = 2.2) +
      ggplot2::geom_text(data = marked, ggplot2::aes(label = .bef_short_labels(.data$label)),
                         hjust = 1.12, vjust = -0.5, size = 3, color = "grey20")
  }
  plot +
    ggplot2::coord_fixed(xlim = limits, ylim = limits) +
    ggplot2::labs(x = "Estimate", y = "Posterior mean")
}

.bef_plot_shrinkage_arrows <- function(sites, centers, K, color_by, has_truth) {
  colors <- .bef_plot_colors()
  rule <- .bef_plot_shrinkage_rule()
  marks <- K <= rule$marks
  head <- if (marks) {
    ggplot2::arrow(length = ggplot2::unit(0.012, "npc"), type = "closed")
  }
  width <- if (marks) 0.4 else 0.25
  # Thin lines need to be less transparent than points to be seen.
  alpha <- if (marks) 1 else 0.4
  plot <- ggplot2::ggplot(sites, ggplot2::aes(y = .data$se)) +
    ggplot2::geom_vline(data = centers, ggplot2::aes(xintercept = .data$center),
                        linetype = "dashed", color = "grey50")
  plain <- identical(color_by, "none")
  plot <- if (plain) {
    plot + ggplot2::geom_segment(
      ggplot2::aes(x = .data$estimate, xend = .data$mean, yend = .data$se),
      color = colors[["fit"]], alpha = alpha, linewidth = width, arrow = head
    )
  } else {
    plot + ggplot2::geom_segment(
      ggplot2::aes(x = .data$estimate, xend = .data$mean, yend = .data$se,
                   color = .data$color),
      alpha = alpha, linewidth = width, arrow = head
    )
  }
  if (marks) {
    plot <- plot +
      ggplot2::geom_point(ggplot2::aes(x = .data$estimate), shape = 1, size = 1.8,
                          color = colors[["data"]])
    if (has_truth) {
      plot <- plot +
        ggplot2::geom_point(ggplot2::aes(x = .data$truth), shape = 4, size = 1.8,
                            color = "#404040")
    }
  } else if (plain) {
    plot <- plot + ggplot2::geom_point(ggplot2::aes(x = .data$mean), size = 0.5,
                                       color = colors[["fit"]])
  } else {
    plot <- plot + ggplot2::geom_point(ggplot2::aes(x = .data$mean, color = .data$color),
                                       size = 0.5)
  }
  marked <- sites[sites$highlight, , drop = FALSE]
  if (nrow(marked) > 0L) {
    plot <- plot +
      ggplot2::geom_segment(
        data = marked,
        ggplot2::aes(x = .data$estimate, xend = .data$mean, yend = .data$se),
        color = colors[["highlight"]], linewidth = 0.8, arrow = head
      )
  }
  labeled <- if (K <= rule$labels) sites else marked
  if (nrow(labeled) > 0L) {
    plot <- plot +
      ggplot2::geom_text(data = labeled,
                         ggplot2::aes(x = .data$estimate, label = .bef_short_labels(.data$label)),
                         vjust = -0.7, size = 2.8, color = "grey20")
  }
  plot + ggplot2::labs(x = "Effect", y = "Standard error")
}

.bef_plot_shrinkage_caption <- function(style, K, color_by, normal, has_truth, highlighted) {
  marks <- K <= .bef_plot_shrinkage_rule()$marks
  .bef_caption(
    if (identical(style, "scatter")) {
      c("Each point is a site: its estimate and the posterior mean of its effect.",
        "Dashed line: the two are equal.")
    } else if (marks) {
      c("Each arrow runs from the estimate of a site (open circle) to the posterior mean of",
        "its effect, at the height of its standard error.",
        "Dashed line: mean of the distribution of effects.")
    } else {
      c("Each line runs from the estimate of a site to the posterior mean of its effect,",
        "which a point marks, at the height of its standard error.",
        "Dashed line: mean of the distribution of effects.")
    },
    if (identical(color_by, "sigma")) "Color: standard error of the estimate.",
    if (identical(color_by, "truth")) "Color: true effect.",
    if (identical(style, "arrows") && has_truth && marks) "Crosses: the true effects.",
    if (!is.null(normal)) {
      if (identical(normal$source, "REML")) {
        paste("Second panel: a normal model of the effects with its mean and standard deviation",
              "estimated by REML and held there.")
      } else {
        "Second panel: a normal model of the effects with the mean and standard deviation given."
      }
    },
    if (highlighted) "Orange: the highlighted sites."
  )
}
