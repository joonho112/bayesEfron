#' Plot the recovery of true effects in simulated data
#'
#' @description
#' For simulated data, where the true effects are known, `plot_recovery()`
#' compares a fit with the truth. It shows the posterior intervals of the
#' sites together with the true effects, the proportion of intervals that
#' contain the true effect at each level, or the posterior means against
#' the true effects.
#'
#' @details
#' `type` selects one of three plots.
#'
#' * `"intervals"` shows the posterior interval and the true effect of each
#'   site, with the sites ordered by posterior mean, and marks the sites
#'   where the true effect lies outside the interval. If there are more
#'   than `n` sites, `n` of them are chosen at random; specify `seed` to
#'   make the choice reproducible.
#' * `"calibration"` shows, for levels from 5% to 95%, the proportion of
#'   all sites with the true effect inside the posterior interval of that
#'   level. The intervals are the ones that [confint.bef_fit_re()] returns.
#'   They are computed from the posterior distribution of each effect and
#'   not from a normal approximation, which would give wrong intervals for
#'   sites with a bimodal posterior distribution.
#' * `"scatter"` shows the posterior means against the true effects. The
#'   caption gives the root mean squared error of the posterior means and
#'   of the estimates.
#'
#' The proportions in these plots describe one data set. They are not the
#' coverage of the intervals, which is the average of such a proportion
#' over repeated data sets. The intervals of all sites depend on the same
#' estimate of \eqn{g}, so they tend to miss together, and the proportion
#' in one data set varies around the nominal level more than it would for
#' independent sites.
#'
#' @inheritParams plot_g
#' @param truth Numeric vector of true effects, one per site, in the same
#'   order as the sites in the data.
#' @param type Type of plot: `"intervals"` (the default), `"calibration"`
#'   or `"scatter"`. See Details.
#' @param level Posterior probability of the intervals, a number between 0
#'   and 1. Only used with `type = "intervals"`.
#' @param n Maximum number of sites shown with `type = "intervals"`, an
#'   integer of at least 1.
#' @param seed Integer seed for the random choice of sites. If `NULL` (the
#'   default), the current state of R's random number generator is used.
#' @param scale Positive number by which the effects are multiplied in the
#'   plot, e.g., `100` to show differences of proportions in percentage
#'   points. Not used with `type = "calibration"`.
#'
#' @return A `ggplot` object.
#'
#' @seealso [site_summary()], [confint.bef_fit_re()],
#'   `vignette("twin-towers")`
#'
#' @examples
#' truth <- twin_towers$theta_true
#' plot_recovery(twin_towers_fit, truth, seed = 1)
#' plot_recovery(twin_towers_fit, truth, type = "calibration")
#' plot_recovery(twin_towers_fit, truth, type = "scatter")
#'
#' @export
plot_recovery <- function(x,
                          truth,
                          type = c("intervals", "calibration", "scatter"),
                          level = 0.9,
                          n = 150,
                          seed = NULL,
                          scale = 1,
                          ...) {
  .bef_check_dots_empty(list(...), "plot_recovery")
  x <- .bef_plot_fit(x)
  type <- .bef_validate_method_choice(type, c("intervals", "calibration", "scatter"),
                                      arg = "type")
  level <- .bef_plot_level(level)
  n <- .bef_plot_count(n, "n", lower = 1L)
  if (!is.null(seed)) {
    seed <- .bef_plot_count(seed, "seed")
  }
  scale <- .bef_validate_scale(scale)
  table <- site_summary(x, level = level)
  K <- nrow(table)
  if (missing(truth) || is.null(truth)) {
    .bef_abort_invalid_args("`truth` must hold the true effect of every site.", arg = "truth")
  }
  truth <- .bef_plot_truth(truth, K)
  table$truth <- truth
  switch(
    type,
    intervals = .bef_plot_recovery_intervals(table, level, n, seed, scale),
    calibration = .bef_plot_recovery_calibration(x, truth),
    scatter = .bef_plot_recovery_scatter(table, scale)
  )
}

# What is said under every plot of recovery.
.bef_recovery_note <- function() {
  "These numbers describe this one data set; they are not the coverage over repeated samples."
}

# `n` of `K` sites chosen at random, in increasing order. With a seed the
# state of R's generator is put back.
.bef_sample_sites <- function(K, n, seed = NULL) {
  if (K <= n) {
    return(seq_len(K))
  }
  if (is.null(seed)) {
    return(sort(sample.int(K, n)))
  }
  old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
  on.exit({
    if (is.null(old_seed)) {
      suppressWarnings(rm(".Random.seed", envir = globalenv(), inherits = FALSE))
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  sort(sample.int(K, n))
}

.bef_plot_recovery_intervals <- function(table, level, n, seed, scale) {
  colors <- .bef_plot_colors()
  K <- nrow(table)
  table$inside <- table$truth >= table$lower & table$truth <= table$upper
  share_all <- mean(table$inside)
  drawn <- table[.bef_sample_sites(K, n, seed), , drop = FALSE]
  drawn <- drawn[order(drawn$mean, -drawn$site), , drop = FALSE]
  drawn$position <- seq_len(nrow(drawn))
  drawn$mark <- ifelse(drawn$inside, "inside", "outside")
  for (column in c("mean", "lower", "upper", "truth")) {
    drawn[[column]] <- drawn[[column]] * scale
  }
  few <- nrow(drawn) <= 60L
  ggplot2::ggplot(drawn, ggplot2::aes(y = .data$position)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$lower, xend = .data$upper, yend = .data$position,
                   color = .data$mark),
      linewidth = if (few) 0.5 else 0.3
    ) +
    ggplot2::geom_point(ggplot2::aes(x = .data$mean, color = .data$mark),
                        size = if (few) 1.4 else 0.6) +
    ggplot2::geom_point(ggplot2::aes(x = .data$truth, shape = .data$mark, fill = .data$mark),
                        color = "#404040", size = if (few) 1.8 else 1.2) +
    ggplot2::scale_color_manual(
      values = c(inside = colors[["fit"]], outside = colors[["highlight"]]), guide = "none"
    ) +
    ggplot2::scale_shape_manual(values = c(inside = 4, outside = 24), guide = "none") +
    ggplot2::scale_fill_manual(
      values = c(inside = "#404040", outside = colors[["highlight"]]), guide = "none"
    ) +
    ggplot2::scale_y_continuous(breaks = NULL, expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::labs(
      x = "Effect", y = NULL,
      caption = .bef_caption(
        if (nrow(drawn) < K) {
          sprintf("Lines: %s posterior intervals of %d of the %s sites, chosen at random, in the order of their posterior means (points).",
                  .bef_percent(level), nrow(drawn), format(K, big.mark = ","))
        } else {
          sprintf("Lines: %s posterior intervals of the %d sites, in the order of their posterior means (points).",
                  .bef_percent(level), K)
        },
        sprintf("Crosses: true effects inside the interval. Orange triangles: true effects outside it (%d of %d).",
                sum(!drawn$inside), nrow(drawn)),
        if (nrow(drawn) < K) {
          sprintf("Of all %s sites, %s have their true effect inside the interval.",
                  format(K, big.mark = ","), .bef_percent(round(share_all, 3L)))
        },
        .bef_recovery_note()
      )
    ) +
    theme_bef()
}

# The levels at which the calibration plot counts the true effects inside
# the intervals.
.bef_recovery_levels <- function() {
  seq(0.05, 0.95, by = 0.05)
}

# The share of the sites whose true effect is inside the posterior interval
# of each level, the intervals being those of confint().
.bef_recovery_inclusion <- function(fit, truth, levels = .bef_recovery_levels()) {
  engine <- .bef_engine(fit)
  weights <- .bef_site_weights(engine)
  lower <- .bef_site_quantile(engine, (1 - levels) / 2, weights = weights)
  upper <- .bef_site_quantile(engine, (1 + levels) / 2, weights = weights)
  data.frame(
    level = levels,
    share = colMeans(truth >= lower & truth <= upper)
  )
}

.bef_plot_recovery_calibration <- function(fit, truth) {
  colors <- .bef_plot_colors()
  inclusion <- .bef_recovery_inclusion(fit, truth)
  K <- length(truth)
  ggplot2::ggplot(inclusion, ggplot2::aes(x = .data$level, y = .data$share)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_line(color = colors[["fit"]], linewidth = 0.6) +
    ggplot2::geom_point(color = colors[["fit"]], size = 1.6) +
    ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
    ggplot2::labs(
      x = "Level of the posterior interval",
      y = "Share of the true effects inside it",
      caption = .bef_caption(
        sprintf("Points: for each level, the share of the %s sites whose true effect is inside the posterior interval of that level.",
                format(K, big.mark = ",")),
        "Dashed line: the share equals the level.",
        .bef_recovery_note()
      )
    ) +
    theme_bef()
}

.bef_plot_recovery_scatter <- function(table, scale) {
  colors <- .bef_plot_colors()
  K <- nrow(table)
  rmse <- function(values) sqrt(mean((values - table$truth)^2)) * scale
  data <- data.frame(truth = table$truth * scale, mean = table$mean * scale)
  limits <- range(data$truth, data$mean)
  ggplot2::ggplot(data, ggplot2::aes(x = .data$truth, y = .data$mean)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(color = colors[["fit"]], alpha = .bef_plot_alpha(K),
                        size = if (K > 100L) 0.9 else 1.8) +
    ggplot2::coord_fixed(xlim = limits, ylim = limits) +
    ggplot2::labs(
      x = "True effect", y = "Posterior mean",
      caption = .bef_caption(
        "Each point is a site. Dashed line: the posterior mean equals the true effect.",
        sprintf("Root mean squared error: %s for the posterior means and %s for the estimates.",
                format(signif(rmse(table$mean), 3L)), format(signif(rmse(table$estimate), 3L))),
        .bef_recovery_note()
      )
    ) +
    theme_bef()
}
