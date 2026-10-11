#' Summarize the estimated distribution of effects
#'
#' @description
#' Computes the posterior mean and a pointwise posterior interval for a
#' feature of \eqn{g}, the distribution of true effects: its density, its
#' distribution function, its probabilities on the grid, its quantiles,
#' the probability of an interval of effects, or its mean, variance and
#' standard deviation. The feature is computed for each posterior draw of
#' \eqn{g} and then summarized over the draws.
#'
#' @details
#' The model estimates \eqn{g} as probabilities on the grid points. For a
#' fit with point masses (`component = "point"`, the default), there are
#' two ways to compute a quantile or a probability from them, and
#' `reading` selects one.
#'
#' * `"interpolated"` spreads the probability of each grid point evenly
#'   over an interval of one grid spacing centered at the point. The
#'   distribution function is then continuous, and the results change
#'   little when the number of grid points changes. For the first and the
#'   last grid point this interval extends half a grid spacing beyond the
#'   limits of the grid, so a quantile can lie up to half a spacing
#'   outside the grid, and a small amount of probability lies beyond each
#'   limit.
#' * `"step"` leaves the probability on the grid point, as the model
#'   specifies, so all results stay on the grid. The distribution function
#'   then jumps at the grid points, and a probability such as
#'   \eqn{P(\theta > 0)}{P(theta > 0)} depends on where the grid points
#'   lie relative to zero.
#'
#' For `type = "density"`, the probabilities of point masses are always
#' divided by the grid spacing, which corresponds to `"interpolated"`.
#'
#' The posterior interval is pointwise: at each value it has posterior
#' probability `level` for that value, not for the whole curve at once.
#'
#' `at`, `lower` and `upper` are specified in the units of the estimates.
#' `scale` only changes the units of the result: effect values are
#' multiplied by it, densities are divided by it, and probabilities are
#' unchanged. For example, if the estimates are differences in
#' proportions, `scale = 100` gives percentage points.
#'
#' @param x A `bef_fit` object created by [bayes_efron_fit()].
#' @param type Feature of \eqn{g} to summarize: `"density"` (the default),
#'   `"cdf"` for the distribution function, `"mass"` for the probabilities
#'   at the grid points, `"quantile"`, `"prob"` for the probability of an
#'   interval of effects, or `"moments"` for the mean, variance and
#'   standard deviation. The standard deviation of \eqn{g} measures the
#'   heterogeneity of the effects.
#' @param level Posterior probability of the pointwise interval, a number
#'   between 0 and 1.
#' @param at Numeric vector of effect values at which the density or the
#'   distribution function is evaluated. If `NULL` (the default), the grid
#'   points are used.
#' @param probs Numeric vector of probabilities for `type = "quantile"`,
#'   each greater than 0 and less than 1. If `NULL` (the default), 0.05,
#'   0.25, 0.5, 0.75 and 0.95 are used.
#' @param lower,upper Limits of the interval of effects for
#'   `type = "prob"`. The result is the probability that an effect is
#'   greater than `lower` and less than or equal to `upper`.
#' @param scale Positive number by which effect values are multiplied in
#'   the result, e.g., `100` to express differences in proportions as
#'   percentage points. See Details.
#' @param reading How quantiles and probabilities are computed for a fit
#'   with point masses. With `"interpolated"` (the default), the
#'   probability of each grid point is spread evenly over one grid spacing
#'   around it. With `"step"`, it stays on the grid point. Use `"step"` if
#'   the effects must stay within given bounds that were specified as the
#'   limits of the grid (see `support` in [make_efron_grid()]). See
#'   Details.
#'
#' @return A data frame with the following columns. The interval level is
#'   stored in the attribute `level`.
#'   \describe{
#'     \item{`value`}{The effect value, for `"density"`, `"cdf"` and
#'       `"mass"`.}
#'     \item{`prob`}{The probability, for `"quantile"`.}
#'     \item{`from`, `to`}{The limits of the interval of effects, for
#'       `"prob"`.}
#'     \item{`parameter`}{The name of the moment, for `"moments"`.}
#'     \item{`mean`}{The posterior mean.}
#'     \item{`lower`, `upper`}{The limits of the pointwise posterior
#'       interval.}
#'   }
#'
#' @seealso [site_summary()], [exceedance_prob()], [plot_g()],
#'   [confint.bef_fit_re()]
#'
#' @examples
#' fit <- krw_fit
#'
#' # The proportion of firms whose callback gap favors white applicants,
#' # and the proportion with a gap above 2 percentage points.
#' g_summary(fit, type = "prob", lower = 0)
#' g_summary(fit, type = "prob", lower = 0.02)
#'
#' # Quantiles of the gaps, in percentage points.
#' g_summary(fit, type = "quantile", probs = c(0.05, 0.5, 0.95), scale = 100)
#'
#' # The density at the grid points, which is what plot_g() shows.
#' head(g_summary(fit, scale = 100))
#'
#' @export
g_summary <- function(x,
                      type = c("density", "cdf", "mass", "quantile", "prob", "moments"),
                      level = 0.9,
                      at = NULL,
                      probs = NULL,
                      lower = -Inf,
                      upper = Inf,
                      scale = 1,
                      reading = c("interpolated", "step")) {
  .bef_require_inherits(x, "bef_fit", "x", "bef_invalid_fit")
  type <- .bef_validate_method_choice(
    type, c("density", "cdf", "mass", "quantile", "prob", "moments"), arg = "type"
  )
  reading <- .bef_validate_method_choice(
    reading, c("interpolated", "step"), arg = "reading"
  )
  level <- .bef_validate_summary_level(level)
  scale <- .bef_validate_scale(scale)
  x <- .bef_prepare_fit(x)

  if (identical(type, "moments")) {
    out <- .bef_g_moments(x, level = level, scale = scale)
    attr(out, "level") <- level
    return(out)
  }

  engine <- .bef_engine(x)
  # CmdStan writes log_g with a fixed number of digits, so the probabilities
  # of a draw sum to one only to that precision; they are made to sum to one.
  g <- exp(engine$log_g)
  g <- g / rowSums(g)
  out <- switch(
    type,
    mass = {
      .bef_summarize_draws_matrix(g, level, data.frame(value = engine$grid * scale))
    },
    density = {
      values <- .bef_validate_at(at, engine$grid)
      density <- .bef_g_density_draws(engine, g, values) / scale
      .bef_summarize_draws_matrix(density, level, data.frame(value = values * scale))
    },
    cdf = {
      values <- .bef_validate_at(at, engine$grid)
      cdf <- .bef_g_cdf_draws(engine, g, values, reading)
      .bef_summarize_draws_matrix(cdf, level, data.frame(value = values * scale))
    },
    quantile = {
      probs <- .bef_validate_probs(probs)
      quantiles <- .bef_g_quantile_draws(engine, g, probs, reading) * scale
      .bef_summarize_draws_matrix(quantiles, level, data.frame(prob = probs))
    },
    prob = {
      limits <- .bef_validate_prob_limits(lower, upper)
      cdf <- .bef_g_cdf_draws(engine, g, limits, reading)
      .bef_summarize_draws_matrix(
        cdf[, 2L, drop = FALSE] - cdf[, 1L, drop = FALSE], level,
        data.frame(from = limits[1L] * scale, to = limits[2L] * scale)
      )
    }
  )
  attr(out, "level") <- level
  out
}

# The posterior mean and the equal-tailed interval of each column of a
# matrix with one row per draw, next to the columns of `labels`.
.bef_summarize_draws_matrix <- function(draws, level, labels) {
  probs <- .bef_interval_probs(level)
  limits <- apply(draws, 2L, stats::quantile, probs = probs, names = FALSE)
  if (is.null(dim(limits))) limits <- matrix(limits, nrow = 2L)
  data.frame(
    labels,
    mean = colMeans(draws),
    lower = limits[1L, ],
    upper = limits[2L, ],
    row.names = NULL,
    check.names = FALSE
  )
}

# The density of g at `values` for every draw: one row per draw and one
# column per value. Point masses are divided by the spacing of the grid and
# are zero outside the cells; normal components give a mixture density.
.bef_g_density_draws <- function(engine, g, values) {
  if (engine$tau_c > 0) {
    kernel <- vapply(values, function(value) {
      stats::dnorm(value, engine$grid, engine$tau_c)
    }, numeric(engine$L))
    return(g %*% kernel)
  }
  d <- engine$spacing
  index <- round((values - engine$grid[1L]) / d) + 1L
  inside <- index >= 1L & index <= engine$L &
    abs(values - engine$grid[pmin(pmax(index, 1L), engine$L)]) <= d / 2 * (1 + 1e-9)
  out <- matrix(0, nrow = nrow(g), ncol = length(values))
  out[, inside] <- g[, index[inside], drop = FALSE] / d
  out
}

# P(theta <= value) under g for every draw: one row per draw and one column
# per value.
.bef_g_cdf_draws <- function(engine, g, values, reading) {
  finite <- is.finite(values)
  below <- matrix(0, nrow = engine$L, ncol = length(values))
  below[, !finite & values > 0] <- 1
  if (any(finite)) {
    below[, finite] <- if (engine$tau_c > 0) {
      vapply(values[finite], function(value) {
        stats::pnorm((value - engine$grid) / engine$tau_c)
      }, numeric(engine$L))
    } else {
      .bef_point_mass_below(engine, values[finite], reading)
    }
  }
  out <- g %*% below
  pmin(pmax(out, 0), 1)
}

# Quantiles of g for every draw: one row per draw and one column per
# probability.
.bef_g_quantile_draws <- function(engine, g, probs, reading) {
  n_draws <- nrow(g)
  if (engine$tau_c == 0) {
    # The draws take the place of the sites: each row of g is a set of
    # probabilities on the grid.
    return(.bef_site_quantile(engine, probs, weights = g, sites = seq_len(n_draws),
                              reading = reading))
  }
  vapply(probs, function(p) {
    lower <- rep(engine$grid[1L] - 10 * engine$tau_c, n_draws)
    upper <- rep(engine$grid[engine$L] + 10 * engine$tau_c, n_draws)
    for (step in seq_len(60L)) {
      middle <- (lower + upper) / 2
      cdf <- rowSums(g * stats::pnorm(outer(middle, engine$grid, "-") / engine$tau_c))
      below <- cdf < p
      lower[below] <- middle[below]
      upper[!below] <- middle[!below]
    }
    (lower + upper) / 2
  }, numeric(n_draws))
}

# The mean, variance and standard deviation of g: posterior means and
# intervals from the draws kept in the fit.
.bef_g_moments <- function(x, level, scale = 1) {
  probs <- .bef_interval_probs(level)
  parameters <- c("mean_g", "var_g", "sd_g")
  factor <- c(mean_g = scale, var_g = scale^2, sd_g = scale)
  rows <- lapply(parameters, function(parameter) {
    draws <- x$posterior[[parameter]] * factor[[parameter]]
    limits <- stats::quantile(draws, probs = probs, names = FALSE)
    data.frame(parameter = parameter, mean = mean(draws), lower = limits[1L],
               upper = limits[2L])
  })
  do.call(rbind, rows)
}

.bef_validate_scale <- function(scale) {
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    .bef_abort_invalid_args("`scale` must be a single positive number.", arg = "scale")
  }
  scale
}

.bef_validate_at <- function(at, grid) {
  if (is.null(at)) {
    return(grid)
  }
  if (!is.numeric(at) || length(at) == 0L || any(!is.finite(at))) {
    .bef_abort_invalid_args(
      "`at` must be NULL or a vector of finite numbers.",
      arg = "at"
    )
  }
  as.numeric(at)
}

.bef_validate_probs <- function(probs) {
  if (is.null(probs)) {
    return(c(0.05, 0.25, 0.5, 0.75, 0.95))
  }
  if (!is.numeric(probs) || length(probs) == 0L || anyNA(probs) ||
      any(probs <= 0) || any(probs >= 1)) {
    .bef_abort_invalid_args(
      "`probs` must be NULL or probabilities greater than 0 and less than 1.",
      arg = "probs"
    )
  }
  as.numeric(probs)
}

.bef_validate_prob_limits <- function(lower, upper) {
  ok <- function(value) is.numeric(value) && length(value) == 1L && !is.na(value)
  if (!ok(lower) || !ok(upper) || lower >= upper) {
    .bef_abort_invalid_args(
      "`lower` and `upper` must be single numbers with `lower` less than `upper`.",
      arg = "lower"
    )
  }
  c(lower, upper)
}
