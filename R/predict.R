# predict() and fitted() methods. Both use the draws stored in the fit and do
# not run Stan.

#' Posterior predictive distribution for a new site
#'
#' @description
#' Draws from the posterior predictive distribution for a site that was
#' not in the data: either its true effect, or the estimate that a new
#' study with a given standard error would produce. For a new site whose
#' estimate is already known, `predict()` instead returns the mean,
#' standard deviation, interval and mode of the posterior distribution of
#' its effect.
#'
#' @details
#' `predict()` has three uses, which are selected by `type` and `newdata`.
#'
#' * `type = "theta"` (the default) draws the true effect of a new site
#'   from \eqn{g}. Each draw uses a different posterior draw of \eqn{g},
#'   so the draws describe the distribution of effects and include the
#'   uncertainty about \eqn{g}.
#' * `type = "theta_hat"` draws the estimate that a new study would
#'   produce: sampling error with the standard error given in `newdata` is
#'   added to a draw of the effect. If `newdata` contains several standard
#'   errors, the result has one column per new site. Each row then uses
#'   one posterior draw of \eqn{g} for all sites and draws a separate
#'   effect and sampling error for each site, so the sites in a row are
#'   independent given that draw.
#' * If `newdata` is a data frame with the columns `estimate` and `se`,
#'   the estimates of the new sites are known. Nothing is drawn, and for
#'   each row `predict()` returns the posterior mean, the posterior
#'   standard deviation, the interval at `level` and the posterior mode of
#'   the effect. `type`, `n`, `jitter` and `seed` are
#'   not used.
#'
#' For a fit with point masses, \eqn{g} is a distribution on the grid, so
#' each draw of an effect is one of the grid points. If `jitter = TRUE`, a
#' uniform random number between minus and plus half the grid spacing is
#' added to the draw, which spreads the draws over the intervals between
#' the grid points. For a fit with `component = "normal"`, a draw comes
#' from the normal component at the grid point that was drawn, and
#' `jitter` is not used.
#'
#' A new site with a known estimate was not in the data, so before its
#' estimate is observed the distribution of its effect is the posterior
#' mean of \eqn{g}. Bayes' rule with the estimate then gives the posterior
#' distribution that `predict()` summarizes. The summaries are those that
#' the site would have in a fit of the data with that site added, on the
#' same grid, up to Monte Carlo error. They differ from the summaries that
#' [site_summary()] gives for a site in the fit with the same estimate and
#' standard error, because that site's estimate was already used to
#' estimate \eqn{g}. The difference is small when there are many sites and
#' is largest for an estimate in a tail of the distribution. Each row of
#' `newdata` is treated as the only new site, so the other rows do not
#' affect its result. For a fit with point masses, the interval limits
#' are computed as in [site_summary()] with its default `reading`.
#'
#' @param object A `bef_fit` object created by [bayes_efron_fit()].
#' @param newdata For `type = "theta_hat"`, a numeric vector of positive
#'   standard errors, one per new site. It is required for that type and
#'   not used for `type = "theta"`. Can also be a data frame with the
#'   columns `estimate` and `se`, for new sites whose estimates are known;
#'   see Details.
#' @param ... These dots are for future extensions and must be empty.
#' @param type Quantity to draw for a new site: `"theta"` (the default)
#'   for its true effect, or `"theta_hat"` for its estimate.
#' @param n Number of draws, a positive integer. If `NULL` (the default),
#'   the number of posterior draws in the fit is used. If `n` is larger
#'   than that, posterior draws of \eqn{g} are used more than once.
#' @param jitter If `TRUE`, the draws from a fit with point masses are
#'   spread between the grid points; see Details.
#' @param seed Integer seed for the draws. If `NULL` (the default), no
#'   seed is set. If a seed is specified, the state of R's random number
#'   generator is restored afterward.
#' @param level Posterior probability of the interval in
#'   the columns `lower` and `upper`, a number between 0 and 1. Only used
#'   if `newdata` is a data frame.
#'
#' @return For `type = "theta"`, and for `type = "theta_hat"` with one
#'   standard error in `newdata`, a numeric vector of `n` draws. For
#'   `type = "theta_hat"` with more than one standard error, a matrix with
#'   `n` rows and one column per standard error.
#'
#'   If `newdata` is a data frame, a data frame with one row per row of
#'   `newdata` and the following columns.
#'   \describe{
#'     \item{`estimate`, `se`}{The estimate and its standard error, as
#'       given.}
#'     \item{`mean`, `sd`}{The posterior mean and the posterior standard
#'       deviation of the effect.}
#'     \item{`lower`, `upper`}{The limits of the posterior interval at
#'       `level`.}
#'     \item{`mode`}{The posterior mode.}
#'   }
#'
#' @seealso [site_summary()], [bef_fit], [bayes_efron_fit()]
#'
#' @examples
#' # The callback gap of a firm that was not in the data.
#' theta_new <- predict(krw_fit, seed = 1)
#' quantile(theta_new, c(0.05, 0.5, 0.95))
#' mean(theta_new > 0)
#'
#' # An estimate with a standard error of 1 percentage point.
#' estimate_new <- predict(krw_fit, newdata = 0.01, type = "theta_hat",
#'                         seed = 1)
#' quantile(estimate_new, c(0.05, 0.5, 0.95))
#'
#' # Two new firms whose estimates are known: gaps of 0 and of 5 percentage
#' # points, with standard errors of 1 and 2 points.
#' predict(krw_fit, newdata = data.frame(estimate = c(0, 0.05),
#'                                       se = c(0.01, 0.02)))
#'
#' @export
predict.bef_fit_re <- function(object,
                               newdata = NULL,
                               ...,
                               type = c("theta", "theta_hat"),
                               n = NULL,
                               jitter = FALSE,
                               seed = NULL,
                               level = 0.9) {
  .bef_check_no_dots(list(...))
  .bef_require_inherits(object, "bef_fit_re", "object", "bef_invalid_fit")
  if (.bef_is_v01_fit(object)) {
    .bef_abort_v01_fit("predict")
  }
  if (is.data.frame(newdata)) {
    return(.bef_predict_new_sites(object, newdata, level = level))
  }
  type <- .bef_validate_method_choice(
    type, c("theta", "theta_hat"), arg = "type"
  )
  if (!is.logical(jitter) || length(jitter) != 1L || is.na(jitter)) {
    .bef_abort_invalid_args(
      "`jitter` must be a single non-missing logical.",
      arg = "jitter"
    )
  }

  grid <- object$metadata$data_list$grid
  tau_c <- object$metadata$data_list$tau_c
  tau_c <- if (is.null(tau_c)) 0 else as.numeric(tau_c)
  g <- .bef_grid_draws(object, "g")
  n_post <- nrow(g)

  if (is.null(n)) {
    n <- n_post
  }
  if (!.bef_is_whole_number(n) || n < 1L) {
    .bef_abort_invalid_args(
      "`n` must be a positive integer scalar.",
      arg = "n"
    )
  }
  n <- as.integer(n)

  if (!is.null(seed)) {
    if (!.bef_is_whole_number(seed)) {
      .bef_abort_invalid_args(
        "`seed` must be NULL or an integer scalar.",
        arg = "seed"
      )
    }
    # Leave R's random number generator as it was found: put the old state
    # back, or remove the state if there was none.
    had_state <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    old <- if (had_state) get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL
    set.seed(as.integer(seed))
    on.exit({
      if (had_state) {
        assign(".Random.seed", old, .GlobalEnv)
      } else {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
  }

  # one posterior draw of g per predictive draw
  rows <- if (n <= n_post) sample.int(n_post, n) else sample.int(n_post, n, replace = TRUE)
  # An effect is the grid point that was drawn, moved within its normal
  # component or, with `jitter`, within its cell.
  spread <- function(theta) {
    if (tau_c > 0) {
      theta + stats::rnorm(length(theta), 0, tau_c)
    } else if (jitter) {
      h <- if (length(grid) > 1L) diff(grid)[1L] else 0
      theta + stats::runif(length(theta), -h / 2, h / 2)
    } else {
      theta
    }
  }
  theta_new <- spread(vapply(
    rows,
    function(r) grid[sample.int(length(grid), 1L, prob = g[r, ])],
    numeric(1)
  ))

  if (identical(type, "theta")) {
    return(theta_new)
  }

  if (is.null(newdata)) {
    .bef_abort_invalid_args(
      "`type = \"theta_hat\"` needs `newdata`, the standard error of the new site.",
      arg = "newdata"
    )
  }
  sigma_labels <- names(newdata)
  sigma_new <- .bef_assert_numeric_vector(newdata, "newdata", min_len = 1L)
  if (any(sigma_new <= 0)) {
    .bef_abort_invalid_args(
      "`newdata` must contain strictly positive standard errors.",
      arg = "newdata"
    )
  }

  if (length(sigma_new) == 1L) {
    return(stats::rnorm(n, theta_new, sigma_new))
  }
  # The first site has the effect drawn above, and its estimates are drawn
  # first, so that they are those of a call with that site alone. The
  # other sites have effects of their own, drawn for each row from the same
  # draw of g in one call.
  first <- stats::rnorm(n, theta_new, sigma_new[1L])
  others <- length(sigma_new) - 1L
  index <- vapply(
    rows,
    function(r) sample.int(length(grid), others, replace = TRUE, prob = g[r, ]),
    integer(others)
  )
  theta_others <- spread(matrix(grid[index], nrow = others))
  out <- unname(cbind(
    first,
    t(theta_others + stats::rnorm(length(theta_others), 0, sigma_new[-1L]))
  ))
  colnames(out) <- sigma_labels
  out
}

#' @order 8
#' @describeIn bef_fit Returns the posterior means of the site effects,
#'   as `coef(object)` does.
#' @export
fitted.bef_fit_re <- function(object, ...) {
  stats::coef(object, type = "mean")
}

# The posterior distributions of the effects of sites that were not in the
# data: the posterior mean of g is the prior distribution of a new site, and
# each new site is conditioned on its own estimate. This is what averaging
# the site posteriors over the draws of g gives when a draw is weighted by
# the marginal likelihood that it gives the estimate, and it is what a fit
# with the site added gives.
.bef_predict_new_sites <- function(object, newdata, level = 0.9) {
  level <- .bef_validate_summary_level(level)
  missing <- setdiff(c("estimate", "se"), names(newdata))
  if (length(missing) > 0L) {
    .bef_abort_invalid_args(
      "A data frame given as `newdata` must have the columns `estimate` and `se`.",
      arg = "newdata",
      missing_fields = missing
    )
  }
  estimate <- newdata$estimate
  se <- newdata$se
  if (!is.numeric(estimate) || !is.numeric(se) || nrow(newdata) == 0L ||
      any(!is.finite(estimate)) || any(!is.finite(se)) || any(se <= 0)) {
    .bef_abort_invalid_args(
      "`newdata$estimate` must be finite and `newdata$se` finite and positive.",
      arg = "newdata"
    )
  }
  object <- .bef_prepare_fit(object)
  fitted <- .bef_engine(object)
  # Normalize each draw on the log scale, then take the posterior mean g.
  # Keeping small probabilities on the log scale avoids losing components
  # that can be relevant to a new estimate far from the original data.
  top <- apply(fitted$log_g, 1L, max)
  log_g <- fitted$log_g - (top + log(rowSums(exp(fitted$log_g - top))))
  log_mean_g <- apply(log_g, 2L, function(x) {
    largest <- max(x)
    if (!is.finite(largest)) return(largest)
    largest + log(mean(exp(x - largest)))
  })
  engine <- .bef_engine_from_parts(
    matrix(log_mean_g, nrow = 1L), fitted$grid, estimate, se, fitted$tau_c
  )
  weights <- .bef_site_weights(engine)
  moments <- .bef_site_moments(engine)
  limits <- .bef_site_quantile(engine, .bef_interval_probs(level), weights = weights)
  out <- data.frame(
    estimate = estimate,
    se = se,
    mean = moments$mean,
    sd = sqrt(moments$var),
    lower = limits[, 1L],
    upper = limits[, 2L],
    mode = .bef_site_mode(engine, weights = weights),
    row.names = rownames(newdata)
  )
  attr(out, "level") <- level
  out
}
