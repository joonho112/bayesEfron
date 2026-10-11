#' Plot the posterior distribution of individual site effects
#'
#' @description
#' `plot_site_posterior()` shows the posterior distribution of the effect
#' of each selected site, together with the two components that determine
#' it: the likelihood of the site's estimate and the estimated distribution
#' of effects \eqn{g}. A posterior distribution with two modes, or with
#' substantial mass at a limit of the grid, can be seen in this plot but
#' not in a posterior mean and an interval.
#'
#' @details
#' The posterior distribution of an effect is proportional to the
#' likelihood times \eqn{g}, averaged over the posterior draws of
#' \eqn{g}. If the likelihood is wide and \eqn{g} has two groups, the
#' posterior distribution has two modes and its mean lies between them, at
#' a value that the effect is unlikely to have.
#'
#' Each site is shown in its own panel. All panels have the same
#' horizontal axis, but their vertical axes differ.
#'
#' For point masses, the curves for the posterior distribution and for
#' \eqn{g} are the probabilities of the grid points divided by the grid
#' spacing, as in [plot_g()].
#'
#' @inheritParams plot_g
#' @param sites Sites to plot, at most 12, specified as a numeric vector of
#'   positions in the data or a character vector of site labels.
#' @param likelihood If `TRUE` (the default), the likelihood of the
#'   estimate is shown. It is a normal density with the estimate as its
#'   mean and the standard error as its standard deviation.
#' @param prior If `TRUE` (the default), the posterior mean of \eqn{g} is
#'   shown.
#'
#' @return A `ggplot` object.
#'
#' @seealso [site_summary()], [plot_g()], `vignette("plots")`
#'
#' @examples
#' # Two imprecise estimates, whose posterior distributions have a mode in
#' # each group, and the most precise one.
#' plot_site_posterior(twin_towers_fit, sites = c(268, 959, 304))
#'
#' # The two sites with the most posterior mass at the lower and at the
#' # upper end of the grid.
#' ends <- diagnose(twin_towers_fit)$grid$ends
#' plot_site_posterior(twin_towers_fit, sites = ends$site)
#'
#' @export
plot_site_posterior <- function(x, sites, likelihood = TRUE, prior = TRUE, scale = 1, ...) {
  .bef_check_dots_empty(list(...), "plot_site_posterior")
  x <- .bef_plot_fit(x)
  likelihood <- .bef_plot_flag(likelihood, "likelihood")
  prior <- .bef_plot_flag(prior, "prior")
  scale <- .bef_validate_scale(scale)
  table <- site_summary(x)
  if (missing(sites) || is.null(sites)) {
    .bef_abort_invalid_args("`sites` must name the sites to draw, by position or by label.",
                            arg = "sites")
  }
  chosen <- .bef_plot_sites(sites, table$label)
  rule <- .bef_plot_site_posterior_rule()
  if (length(chosen) > rule$panels) {
    .bef_abort_invalid_args(
      sprintf("`sites` must name at most %d sites; it names %d.", rule$panels, length(chosen)),
      arg = "sites"
    )
  }

  engine <- .bef_engine(x)
  estimate <- engine$theta_hat[chosen]
  se <- engine$sigma[chosen]
  # The axis holds the part of the grid where g lies and three standard
  # errors on each side of every estimate.
  core <- .bef_plot_g_limits(engine, 1)
  limits <- c(min(core[1L], estimate - 3 * se), max(core[2L], estimate + 3 * se))
  fine <- seq(limits[1L], limits[2L], length.out = rule$points)
  weights <- .bef_site_weights(engine, sites = chosen)
  g <- colMeans(exp(engine$log_g))
  smooth <- engine$tau_c > 0
  inside <- engine$grid >= limits[1L] & engine$grid <= limits[2L]
  components <- if (smooth) .bef_engine_components(engine, chosen)

  names <- ifelse(table$label[chosen] == as.character(table$site[chosen]),
                  sprintf("Site %s", table$label[chosen]), table$label[chosen])
  number <- function(values) vapply(values, function(v) format(signif(v, 3L)), character(1L))
  panels <- sprintf("%s: estimate %s, standard error %s", .bef_short_labels(names),
                    number(estimate * scale), number(se * scale))
  curves <- c("Posterior of the effect", "Likelihood of the estimate", "Estimate of g")
  rows <- lapply(seq_along(chosen), function(i) {
    posterior <- if (smooth) {
      data.frame(value = fine, density = vapply(fine, function(at) {
        sum(weights[i, ] * stats::dnorm(at, components$mean[i, ], components$sd[i]))
      }, numeric(1L)))
    } else {
      data.frame(value = engine$grid[inside], density = weights[i, inside] / engine$spacing)
    }
    parts <- list(cbind(curve = curves[[1L]], posterior))
    if (likelihood) {
      parts[[2L]] <- data.frame(curve = curves[[2L]], value = fine,
                                density = stats::dnorm(fine, estimate[[i]], se[[i]]))
    }
    if (prior) {
      parts[[3L]] <- if (smooth) {
        data.frame(curve = curves[[3L]], value = fine, density = vapply(fine, function(at) {
          sum(g * stats::dnorm(at, engine$grid, engine$tau_c))
        }, numeric(1L)))
      } else {
        data.frame(curve = curves[[3L]], value = engine$grid[inside],
                   density = g[inside] / engine$spacing)
      }
    }
    cbind(panel = panels[[i]], do.call(rbind, parts))
  })
  data <- do.call(rbind, rows)
  data$panel <- factor(data$panel, levels = unique(panels))
  data$curve <- factor(data$curve, levels = curves)
  data$value <- data$value * scale
  data$density <- data$density / scale
  colors <- .bef_plot_colors()
  shown <- curves[c(TRUE, likelihood, prior)]

  ggplot2::ggplot(data, ggplot2::aes(x = .data$value, y = .data$density)) +
    ggplot2::geom_area(data = data[data$curve == curves[[1L]], , drop = FALSE],
                       fill = colors[["fit"]], alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(color = .data$curve, linetype = .data$curve),
                       linewidth = 0.6) +
    ggplot2::scale_color_manual(
      values = stats::setNames(c(colors[["fit"]], colors[["data"]], colors[["other1"]]), curves),
      breaks = shown
    ) +
    ggplot2::scale_linetype_manual(
      values = stats::setNames(c("solid", "dashed", "longdash"), curves), breaks = shown
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data$panel), scales = "free_y",
                        ncol = if (length(chosen) <= 4L) 1L else if (length(chosen) <= 8L) 2L else 3L) +
    ggplot2::labs(
      x = "Effect", y = "Density",
      caption = .bef_caption(
        "Solid line and area: posterior distribution of the effect of the site, averaged over the draws of g.",
        if (likelihood) {
          "Dashed line: likelihood of the estimate, a normal density with the standard error as its standard deviation."
        },
        if (prior) "Long dashes: posterior mean of g."
      )
    ) +
    theme_bef()
}

# The largest number of panels, and the number of points at which a curve
# is evaluated.
.bef_plot_site_posterior_rule <- function() {
  list(panels = 12L, points = 401L)
}
