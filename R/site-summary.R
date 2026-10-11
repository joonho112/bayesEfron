#' Summarize the posterior distributions of the site effects
#'
#' @description
#' `site_summary()` returns a data frame with one row per site: the site's
#' estimate and the mean, standard deviation, interval and mode of the
#' posterior distribution of its effect. Optionally, it also contains a
#' second, inner interval and the exceedance probability, i.e., the
#' posterior probability that the effect is greater than a threshold.
#' `exceedance_prob()` returns only that probability.
#'
#' @details
#' Both functions use the posterior distribution of each effect, which
#' averages over the posterior draws of \eqn{g}. They work for fits with
#' or without site-level draws (`site_draws` in [bayes_efron_fit()]).
#'
#' The posterior variance of an effect is the sum of two parts, and the
#' result contains the square root of each.
#'
#' * `sd_given_g` is the square root of the variance that the effect would
#'   have if \eqn{g} were known, averaged over the posterior draws of
#'   \eqn{g}.
#' * `sd_about_g` is the square root of the variance, over those draws, of
#'   the mean of the effect given \eqn{g}. This is the part that the
#'   uncertainty about \eqn{g} adds.
#'
#' The squares of the two columns sum to the square of `sd`. A plug-in
#' analysis, which treats one estimate of \eqn{g} as known, omits the
#' second part. Its standard deviation is still not equal to `sd_given_g`,
#' because the posterior mean of \eqn{g} is more spread out than a single
#' draw of \eqn{g}; for most sites it is larger than `sd_given_g`. See
#' [as_normalmix()] for how to compute it.
#'
#' For a fit with point masses, `reading` applies to the interval limits
#' and to the exceedance probabilities. With `reading = "step"`, a grid
#' point equal to `threshold` is counted as not above it, so the
#' probabilities for `direction = "above"` and `"below"` sum to 1.
#'
#' The result does not include ranks: the rank of a site depends on the
#' effects of all sites jointly, and the effects are independent only
#' given \eqn{g}, so ranks cannot be computed from the separate posterior
#' distributions summarized here. Joint draws of the effects are needed
#' for that (see `simulate()` in [bef_fit]).
#'
#' @param x A `bef_fit` object created by [bayes_efron_fit()].
#' @param level Posterior probability of the interval in
#'   the columns `lower` and `upper`, a number between 0 and 1.
#' @param inner_level Optional level of a second interval, a number
#'   between 0 and 1. If specified, the columns `inner_lower` and
#'   `inner_upper` contain the limits of that interval. This is intended
#'   for plots that show two interval widths.
#' @param threshold Threshold for the exceedance probability, a number in
#'   the units of the estimates, whatever the value of `scale`. In
#'   `site_summary()` it is optional: if specified, the column
#'   `prob_above` contains the posterior probability that the effect is
#'   greater than `threshold`.
#' @param sites Sites to include, as a numeric vector of positions in the
#'   data or a character vector of site labels. If `NULL` (the default),
#'   all sites are included.
#' @param scale Positive number by which the effect columns are multiplied
#'   in the result, e.g., `100` to express differences in proportions as
#'   percentage points.
#' @param reading How quantiles and probabilities are computed for a fit
#'   with point masses: `"interpolated"` (the default) spreads the
#'   probability of each grid point evenly over one grid spacing around
#'   it, and `"step"` leaves it on the grid point. See the Details of
#'   [g_summary()].
#' @param direction Either `"above"` (the default) for the probability
#'   that the effect is greater than `threshold`, or `"below"` for the
#'   probability that it is less than or equal to `threshold`.
#'
#' @return `site_summary()` returns a data frame with one row per site and
#'   the following columns.
#'   \describe{
#'     \item{`site`, `label`}{The site's position in the data and its
#'       label.}
#'     \item{`estimate`, `se`}{The estimate and its standard error, as
#'       given.}
#'     \item{`mean`, `sd`}{The posterior mean and the posterior standard
#'       deviation of the effect.}
#'     \item{`lower`, `upper`}{The limits of the posterior interval at
#'       `level`.}
#'     \item{`mode`}{The posterior mode.}
#'     \item{`sd_given_g`, `sd_about_g`}{The two parts of `sd` described in
#'       Details.}
#'     \item{`inner_lower`, `inner_upper`}{The limits of the interval at
#'       `inner_level`, only if `inner_level` is specified.}
#'     \item{`prob_above`}{The posterior probability that the effect is
#'       greater than `threshold`, only if `threshold` is specified.}
#'   }
#'   The values of `level`, `inner_level` and `threshold` are stored as
#'   attributes.
#'
#'   `exceedance_prob()` returns a numeric vector of probabilities, named
#'   by the site labels.
#'
#' @seealso [g_summary()], [plot_sites()], [confint.bef_fit_re()],
#'   [summary.bef_fit_re()]
#'
#' @examples
#' fit <- krw_fit
#'
#' # The three firms with the largest posterior mean gaps, in percentage
#' # points, with the probability that the gap exceeds 2 points.
#' sites <- site_summary(fit, threshold = 0.02, scale = 100)
#' top <- sites[order(-sites$mean)[1:3], ]
#' top[, c("label", "estimate", "se", "mean", "lower", "upper", "prob_above")]
#'
#' # How many firms have a positive gap with probability at least 0.95.
#' sum(exceedance_prob(fit) >= 0.95)
#'
#' # The part of the posterior standard deviation that comes from not
#' # knowing g is largest for the firms with the largest estimates.
#' largest <- order(-sites$sd_about_g)[1:3]
#' sites[largest, c("label", "estimate", "se", "sd", "sd_given_g", "sd_about_g")]
#'
#' @export
site_summary <- function(x,
                         level = 0.9,
                         inner_level = NULL,
                         threshold = NULL,
                         sites = NULL,
                         scale = 1,
                         reading = c("interpolated", "step")) {
  .bef_require_inherits(x, "bef_fit_re", "x", "bef_invalid_fit")
  level <- .bef_validate_summary_level(level)
  if (!is.null(inner_level)) {
    inner_level <- .bef_validate_summary_level(inner_level)
  }
  threshold <- .bef_validate_threshold(threshold, null_ok = TRUE)
  scale <- .bef_validate_scale(scale)
  reading <- .bef_validate_method_choice(
    reading, c("interpolated", "step"), arg = "reading"
  )
  x <- .bef_prepare_fit(x)

  table <- .bef_site_table(x)
  selected <- .bef_resolve_site_arg(sites, table$label)
  engine <- .bef_engine(x)
  weights <- .bef_site_weights(engine, sites = selected)
  moments <- .bef_site_moments(engine, sites = selected)
  limits <- .bef_site_quantile(engine, .bef_interval_probs(level), weights = weights,
                               sites = selected, reading = reading)
  out <- data.frame(
    site = table$site[selected],
    label = table$label[selected],
    estimate = engine$theta_hat[selected] * scale,
    se = engine$sigma[selected] * scale,
    mean = moments$mean * scale,
    sd = sqrt(moments$var) * scale,
    lower = limits[, 1L] * scale,
    upper = limits[, 2L] * scale,
    mode = .bef_site_mode(engine, weights = weights, sites = selected) * scale,
    sd_given_g = sqrt(moments$var_within) * scale,
    sd_about_g = sqrt(moments$var_between) * scale,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!is.null(inner_level)) {
    inner <- .bef_site_quantile(engine, .bef_interval_probs(inner_level),
                                weights = weights, sites = selected, reading = reading)
    out$inner_lower <- inner[, 1L] * scale
    out$inner_upper <- inner[, 2L] * scale
  }
  if (!is.null(threshold)) {
    out$prob_above <- .bef_site_exceedance(engine, threshold, weights = weights,
                                           sites = selected, reading = reading)
  }
  attr(out, "level") <- level
  attr(out, "inner_level") <- inner_level
  attr(out, "threshold") <- threshold
  out
}

#' @rdname site_summary
#' @export
exceedance_prob <- function(x,
                            threshold = 0,
                            direction = c("above", "below"),
                            sites = NULL,
                            reading = c("interpolated", "step")) {
  .bef_require_inherits(x, "bef_fit_re", "x", "bef_invalid_fit")
  threshold <- .bef_validate_threshold(threshold, null_ok = FALSE)
  direction <- .bef_validate_method_choice(
    direction, c("above", "below"), arg = "direction"
  )
  reading <- .bef_validate_method_choice(
    reading, c("interpolated", "step"), arg = "reading"
  )
  x <- .bef_prepare_fit(x)
  table <- .bef_site_table(x)
  selected <- .bef_resolve_site_arg(sites, table$label)
  engine <- .bef_engine(x)
  above <- .bef_site_exceedance(engine, threshold, sites = selected, reading = reading)
  out <- if (identical(direction, "above")) above else 1 - above
  names(out) <- table$label[selected]
  out
}

# The sites that `sites` selects, with an error that names that argument.
.bef_resolve_site_arg <- function(sites, labels) {
  tryCatch(
    .bef_resolve_site_parm(sites, labels),
    bef_invalid_args = function(err) {
      .bef_abort_invalid_args(
        gsub("`parm`", "`sites`", conditionMessage(err), fixed = TRUE),
        arg = "sites"
      )
    }
  )
}

.bef_validate_threshold <- function(threshold, null_ok) {
  if (is.null(threshold) && isTRUE(null_ok)) {
    return(NULL)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold)) {
    .bef_abort_invalid_args(
      "`threshold` must be a single finite number.",
      arg = "threshold"
    )
  }
  threshold
}
