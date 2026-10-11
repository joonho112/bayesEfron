#' Posterior predictive check
#'
#' @description
#' `plot_ppc()` is a graphical posterior predictive check. It compares the
#' estimates with the distribution that the fitted model predicts for
#' them, through the probability integral transform (PIT) of each estimate
#' or through the marginal density of the estimates.
#'
#' @details
#' Under the model, the estimate of a site is an effect from \eqn{g} plus
#' normal sampling error, so its distribution is \eqn{g} convolved with a
#' normal distribution whose standard deviation is the site's standard
#' error. The PIT of an estimate is the probability that the model assigns
#' to a smaller estimate at that site. If \eqn{g} were known and the model
#' were correct, the PIT values would be uniformly distributed between 0
#' and 1, so the check compares the PIT values of the fitted model with
#' the uniform distribution.
#'
#' For site \eqn{i} with estimate \eqn{\hat\theta_i}{theta_hat_i} and
#' standard error \eqn{\sigma_i}{sigma_i}, and for one posterior draw of
#' the probabilities \eqn{g_j} at the grid points \eqn{\tau_j}{tau_j},
#' \deqn{\mathrm{PIT}_i = \sum_j g_j \,
#'   \Phi\left(\frac{\hat\theta_i - \tau_j}
#'   {\sqrt{\sigma_i^2 + \tau_c^2}}\right),}{
#'   PIT_i = sum_j g_j Phi((theta_hat_i - tau_j) / sqrt(sigma_i^2 + tau_c^2)),}
#' where \eqn{\Phi}{Phi} is the standard normal distribution function and
#' \eqn{\tau_c}{tau_c} is the standard deviation of a normal component
#' (zero for point masses). The plots show the mean of
#' \eqn{\mathrm{PIT}_i}{PIT_i} over the posterior draws.
#'
#' There are three types of plot.
#'
#' * `"pit_ecdf"` shows the empirical distribution function of the PIT
#'   values minus the distribution function of the uniform distribution,
#'   with pointwise reference intervals for independent uniform values.
#' * `"pit_hist"` shows a histogram of the PIT values, with the
#'   corresponding reference intervals for the height of a bar. If there
#'   are fewer than 30 sites, the bars would contain too few values, and
#'   the `"pit_ecdf"` plot is shown instead.
#' * `"marginal"` shows a histogram of the estimates together with the
#'   density that the model assigns to an estimate, averaged over the
#'   standard errors of the sites.
#'
#' The reference intervals are computed for independent uniform values,
#' which is what the PIT values would be if \eqn{g} were known. However,
#' the PIT values of a fitted model all depend on the same estimate of
#' \eqn{g}, and the band ignores this dependence. Its level is therefore
#' not the probability that the line of a fitted model lies inside the
#' band at a given point, with or without `loo = TRUE`.
#'
#' Because the model was fitted to the same estimates that it is checked
#' against, the check is lenient. With `loo = TRUE`, the PIT of each site
#' is computed from the posterior distribution without that site, using
#' Pareto smoothed importance sampling. This requires the loo package.
#'
#' The check can show that the estimates are not what a normal likelihood
#' with the given standard errors and some \eqn{g} on this grid would
#' produce, e.g., because the standard errors are too small, the grid is
#' too narrow, or there are too few spline functions. Note that it cannot
#' show that the shape of \eqn{g} is correct. Distributions of effects
#' with different shapes give nearly the same distribution of noisy
#' estimates, so a model can pass the check with a wrong \eqn{g}.
#'
#' @inheritParams plot_g
#' @param type Type of plot: `"pit_ecdf"` (the default), `"pit_hist"` or
#'   `"marginal"`. See Details.
#' @param loo If `TRUE`, the PIT of each site is computed from the posterior
#'   distribution without that site (leave-one-out). Not used with
#'   `type = "marginal"`.
#' @param level Level of the pointwise reference intervals, a number
#'   between 0 and 1. Not used with `type = "marginal"`.
#' @param scale Positive number by which the estimates are multiplied on
#'   the horizontal axis. Only used with `type = "marginal"`.
#'
#' @return A `ggplot` object.
#'
#' @seealso [plot_diagnostics()], [loo.bef_fit()], `vignette("diagnostics")`
#'
#' @examples
#' plot_ppc(krw_fit)
#' plot_ppc(twin_towers_fit, type = "pit_hist")
#' plot_ppc(krw_fit, type = "marginal", scale = 100)
#'
#' @export
plot_ppc <- function(x,
                     type = c("pit_ecdf", "pit_hist", "marginal"),
                     loo = FALSE,
                     level = 0.95,
                     scale = 1,
                     ...) {
  .bef_check_dots_empty(list(...), "plot_ppc")
  x <- .bef_plot_fit(x)
  type <- .bef_validate_method_choice(type, c("pit_ecdf", "pit_hist", "marginal"), arg = "type")
  loo <- .bef_plot_flag(loo, "loo")
  level <- .bef_plot_level(level)
  scale <- .bef_validate_scale(scale)
  engine <- .bef_engine(x)
  if (identical(type, "marginal")) {
    return(.bef_plot_ppc_marginal(engine, scale))
  }
  if (loo) {
    .bef_require_loo()
  }
  pit <- .bef_pit(engine, loo = loo)
  rule <- .bef_plot_ppc_rule()
  if (identical(type, "pit_hist") && engine$K < rule$histogram) {
    rlang::inform(
      sprintf(
        "With %d sites the bars of a histogram hold too few values; the difference of the distribution functions is drawn in its place.",
        engine$K
      ),
      class = "bef_plot_few_sites"
    )
    type <- "pit_ecdf"
  }
  if (identical(type, "pit_ecdf")) {
    .bef_plot_ppc_ecdf(pit, level, loo)
  } else {
    .bef_plot_ppc_hist(pit, level, loo)
  }
}

# The number of sites below which a histogram of the PIT values is not
# drawn, and the number of points at which curves are evaluated.
.bef_plot_ppc_rule <- function() {
  list(histogram = 30L, points = 201L)
}

# The probability integral transform of each estimate under the fitted
# marginal model: the mean over the draws, or, with `loo`, the mean with the
# Pareto smoothed importance weights that leave the site out.
.bef_pit <- function(engine, loo = FALSE, max_cells = 2e6) {
  s_total <- sqrt(engine$sigma^2 + engine$tau_c^2)
  below <- stats::pnorm(outer(engine$theta_hat, engine$grid, "-") / s_total)
  g <- exp(engine$log_g)
  g <- g / rowSums(g)
  if (!loo) {
    return(as.numeric(below %*% colMeans(g)))
  }
  out <- numeric(engine$K)
  for (block in .bef_engine_blocks(seq_len(engine$K), nrow(g), max_cells)) {
    log_lik <- .bef_log_lik(engine, sites = block)
    smoothed <- suppressWarnings(loo::psis(-log_lik, r_eff = NA))
    weights <- stats::weights(smoothed, log = FALSE, normalize = TRUE)
    out[block] <- colSums(weights * (g %*% t(below[block, , drop = FALSE])))
  }
  out
}

.bef_plot_ppc_what <- function(loo) {
  if (loo) {
    "PIT of the estimate, computed without its site"
  } else {
    "PIT of the estimate"
  }
}

.bef_plot_ppc_ecdf <- function(pit, level, loo) {
  colors <- .bef_plot_colors()
  K <- length(pit)
  at <- seq(0, 1, length.out = .bef_plot_ppc_rule()$points)
  tail <- (1 - level) / 2
  band <- data.frame(
    value = at,
    lower = stats::qbinom(tail, K, at) / K - at,
    upper = stats::qbinom(1 - tail, K, at) / K - at
  )
  curve <- data.frame(value = at, difference = stats::ecdf(pit)(at) - at)
  ggplot2::ggplot(curve, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_ribbon(data = band, ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = "grey85") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_step(ggplot2::aes(y = .data$difference), color = colors[["fit"]],
                       linewidth = 0.6) +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +
    ggplot2::labs(
      x = .bef_plot_ppc_what(loo), y = "Distribution function less the uniform one",
      caption = .bef_caption(
        sprintf("Line: the empirical distribution function of the PIT values of the %s estimates, less that of the uniform distribution.",
                format(K, big.mark = ",")),
        sprintf("Band: pointwise %s reference intervals for independent uniform values.",
                .bef_percent(level)),
        "The reference does not account for fitting g or dependence among PIT values.",
        if (loo) "Each value is computed from the posterior distribution without its site."
      )
    ) +
    theme_bef()
}

.bef_plot_ppc_hist <- function(pit, level, loo) {
  colors <- .bef_plot_colors()
  K <- length(pit)
  bars <- as.integer(min(20, max(5, K %/% 15)))
  breaks <- seq(0, 1, length.out = bars + 1L)
  counts <- graphics::hist(pmin(pmax(pit, 0), 1), breaks = breaks, plot = FALSE)$counts
  tail <- (1 - level) / 2
  data <- data.frame(xmin = breaks[-length(breaks)], xmax = breaks[-1L], count = counts)
  band <- data.frame(xmin = 0, xmax = 1, lower = stats::qbinom(tail, K, 1 / bars),
                     upper = stats::qbinom(1 - tail, K, 1 / bars))
  ggplot2::ggplot(data) +
    ggplot2::geom_rect(data = band,
                       ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$lower,
                                    ymax = .data$upper),
                       fill = "grey85") +
    ggplot2::geom_rect(ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = 0,
                                    ymax = .data$count),
                       fill = colors[["fit"]], alpha = 0.55, color = "white", linewidth = 0.2) +
    ggplot2::geom_hline(yintercept = K / bars, linetype = "dashed", color = "grey35") +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +
    ggplot2::labs(
      x = .bef_plot_ppc_what(loo), y = "Number of sites",
      caption = .bef_caption(
        sprintf("Bars: the PIT values of the %s estimates.", format(K, big.mark = ",")),
        "Dashed line: the height of every bar under the uniform distribution.",
        sprintf("Band: pointwise %s reference intervals for bin counts from independent uniform values.",
                .bef_percent(level)),
        "The reference does not account for fitting g or dependence among PIT values.",
        if (loo) "Each value is computed from the posterior distribution without its site."
      )
    ) +
    theme_bef()
}

.bef_plot_ppc_marginal <- function(engine, scale) {
  colors <- .bef_plot_colors()
  K <- engine$K
  s_total <- sqrt(engine$sigma^2 + engine$tau_c^2)
  g <- colMeans(exp(engine$log_g))
  g <- g / sum(g)
  reach <- range(engine$theta_hat) + c(-1, 1) * 2 * stats::median(s_total)
  at <- seq(reach[1L], reach[2L], length.out = .bef_plot_ppc_rule()$points)
  # The density of an estimate at each site, then the mean over the sites.
  density <- numeric(length(at))
  distance <- outer(at, engine$grid, "-")
  for (i in seq_len(K)) {
    density <- density + as.numeric(stats::dnorm(distance, sd = s_total[[i]]) %*% g)
  }
  curve <- data.frame(value = at * scale, density = density / K / scale)
  values <- engine$theta_hat * scale
  bins <- graphics::hist(values, breaks = .bef_plot_bins(values), plot = FALSE)
  bars <- data.frame(xmin = bins$breaks[-length(bins$breaks)], xmax = bins$breaks[-1L],
                     density = bins$density)
  ggplot2::ggplot(curve) +
    ggplot2::geom_rect(data = bars,
                       ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = 0,
                                    ymax = .data$density),
                       fill = "grey75", color = "white", linewidth = 0.2) +
    ggplot2::geom_line(ggplot2::aes(x = .data$value, y = .data$density),
                       color = colors[["fit"]], linewidth = 0.7) +
    ggplot2::labs(
      x = "Estimate", y = "Density",
      caption = .bef_caption(
        sprintf("Bars: the %s estimates.", format(K, big.mark = ",")),
        "Line: the density that the fitted model gives to an estimate, averaged over the standard errors of the sites and the posterior draws of g."
      )
    ) +
    theme_bef()
}
