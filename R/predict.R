# predict() and fitted() methods. Both use the draws stored in the fit and do
# not run Stan.

#' Predict the effect or the estimate at a new site
#'
#' @description
#' `predict()` draws from the predictive distribution for a site that was
#' not in the data. With `type = "theta"` it draws the true effect of the
#' new site from \eqn{g}, using a different posterior draw of \eqn{g} each
#' time, so that the draws describe the distribution of effects with the
#' uncertainty about \eqn{g} included. With `type = "theta_hat"` it adds
#' sampling error with the standard error given in `newdata`, and so draws
#' the estimate that a new study of that precision would give.
#'
#' @details
#' \eqn{g} is a distribution on the grid, so a draw of the effect is one of
#' the grid points. With `jitter = TRUE` a uniform random number of up to
#' half the grid spacing is added or subtracted, which spreads the draws
#' over the intervals between the grid points.
#'
#' @param object A `bef_fit_re` object from [bayes_efron_fit()].
#' @param newdata Numeric vector of positive standard errors, one for each
#'   new site. It is needed for `type = "theta_hat"` and not used
#'   otherwise.
#' @param ... These dots are for future extensions and must be empty.
#' @param type Character string: `"theta"` (the default) for the true
#'   effect of a new site, or `"theta_hat"` for its estimate.
#' @param n Whole number, the number of draws. `NULL` (the default) means
#'   the number of posterior draws in the fit. If `n` is larger than that,
#'   posterior draws of \eqn{g} are used more than once.
#' @param jitter Logical. If `TRUE`, the draws are spread between the grid
#'   points, as described in Details. Defaults to `FALSE`.
#' @param seed Whole number, a seed for the draws, or `NULL` (the default).
#'   When a seed is given, R's random number generator is left afterward in
#'   the state in which it was found.
#'
#' @return A numeric vector of `n` draws. For `type = "theta_hat"` with
#'   more than one standard error in `newdata`, a matrix with `n` rows and
#'   one column for each standard error.
#'
#' @seealso [bef_fit], [bayes_efron_fit()]
#'
#' @examples
#' # The effect in a new study.
#' theta_new <- predict(raudenbush_fit, seed = 1)
#' quantile(theta_new, c(0.05, 0.5, 0.95))
#' mean(theta_new > 0)
#'
#' # The estimate that a new study with standard error 0.15 would give.
#' estimate_new <- predict(raudenbush_fit, newdata = 0.15, type = "theta_hat",
#'                         seed = 1)
#' quantile(estimate_new, c(0.05, 0.5, 0.95))
#'
#' @export
predict.bef_fit_re <- function(object,
                               newdata = NULL,
                               ...,
                               type = c("theta", "theta_hat"),
                               n = NULL,
                               jitter = FALSE,
                               seed = NULL) {
  .bef_check_no_dots(list(...))
  .bef_require_inherits(object, "bef_fit_re", "object", "bef_invalid_fit")
  if (.bef_is_v01_fit(object)) {
    .bef_abort_v01_fit("predict")
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
  theta_new <- vapply(
    rows,
    function(r) grid[sample.int(length(grid), 1L, prob = g[r, ])],
    numeric(1)
  )

  if (jitter) {
    h <- if (length(grid) > 1L) diff(grid)[1L] else 0
    theta_new <- theta_new + stats::runif(n, -h / 2, h / 2)
  }

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
  out <- matrix(
    vapply(
      sigma_new,
      function(s) stats::rnorm(n, theta_new, s),
      numeric(n)
    ),
    nrow = n,
    ncol = length(sigma_new)
  )
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
