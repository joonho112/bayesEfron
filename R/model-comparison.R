# loo() and waic() methods.

#' Leave-one-out cross-validation for a fitted model
#'
#' @description
#' `loo()` estimates how well the fitted model predicts the estimate of a
#' site that is left out. It applies Pareto smoothed importance sampling
#' (Vehtari, Gelman and Gabry, 2017) to the log likelihood of each site,
#' which is stored with the draws. The result serves to compare fits of the
#' same data that differ in the grid or in the number of spline functions;
#' [compare_efron_grids()] makes that comparison in one call.
#'
#' `waic()` computes the widely applicable information criterion from the
#' same log likelihood. We recommend `loo()`, because its Pareto k values
#' show for each site whether the approximation can be relied on, whereas
#' `waic()` gives only a general warning when some of its terms are large.
#'
#' @details
#' Both methods need the loo package. The relative efficiency that
#' [loo::loo()] uses for its Monte Carlo standard errors is computed from
#' the draws with [loo::relative_eff()].
#'
#' @param x A `bef_fit` object from [bayes_efron_fit()].
#' @param ... Passed to [loo::loo()] or [loo::waic()].
#'
#' @return The object that [loo::loo()] or [loo::waic()] returns.
#'
#' @references
#' Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. *Statistics
#' and Computing*, 27(5), 1413-1432. \doi{10.1007/s11222-016-9696-4}
#'
#' @seealso [compare_efron_grids()], [plot.bef_fit_re()] with
#'   `type = "loo"`
#'
#' @examplesIf requireNamespace("loo", quietly = TRUE)
#' fit_loo <- loo::loo(raudenbush_fit)
#' fit_loo
#' range(fit_loo$diagnostics$pareto_k)
#'
#' @exportS3Method loo::loo
loo.bef_fit <- function(x, ...) {
  .bef_require_loo()
  .bef_require_inherits(x, "bef_fit", "x", "bef_invalid_fit")
  if (.bef_is_v01_fit(x)) {
    .bef_abort_v01_fit("loo")
  }

  log_lik <- .bef_log_lik_draws(x)
  if (any(!is.finite(log_lik))) {
    .bef_abort_invalid_fit(
      paste0(
        "The log likelihood of some draws is not finite, so `loo()` cannot ",
        "be computed. Check the data and the grid, and fit the model again."
      )
    )
  }
  # Subtracting a constant from the log likelihood of a site changes neither
  # the relative efficiency nor the normalized importance weights. It keeps
  # exp(log_lik), which loo also forms for its Monte Carlo errors, away from
  # underflow. The constants are added back below.
  offsets <- unname(apply(log_lik, 2L, max))
  centered <- sweep(log_lik, 2L, offsets, "-")
  r_eff <- loo::relative_eff(
    exp(centered),
    chain_id = .bef_draw_chain_id(x$draws)
  )
  result <- loo::loo(centered, r_eff = r_eff, ...)
  .bef_restore_loo_scale(result, offsets)
}

# Put the scores and any saved importance weights back on the scale of the
# log likelihood. The standard errors over sites change with the constants,
# so they are computed again from the pointwise values.
.bef_restore_loo_scale <- function(result, offsets) {
  result$pointwise[, "elpd_loo"] <- result$pointwise[, "elpd_loo"] + offsets
  result$pointwise[, "looic"] <- -2 * result$pointwise[, "elpd_loo"]
  for (metric in c("elpd_loo", "p_loo", "looic")) {
    values <- result$pointwise[, metric]
    result$estimates[metric, ] <- c(
      sum(values), sqrt(length(values) * stats::var(values))
    )
    # loo keeps these copies of the estimates for older code.
    result[[metric]] <- result$estimates[metric, "Estimate"]
    result[[paste0("se_", metric)]] <- result$estimates[metric, "SE"]
  }

  # `is_method` comes through the dots, so the saved object can be any of
  # the three. Its log weights and their normalizing constant are shifted
  # together, which keeps weights() consistent.
  for (slot in c("psis_object", "tis_object", "sis_object")) {
    object <- result[[slot]]
    if (!is.null(object)) {
      object$log_weights <- sweep(object$log_weights, 2L, offsets, "-")
      attr(object, "norm_const_log") <- attr(object, "norm_const_log") - offsets
      result[[slot]] <- object
    }
  }
  result
}

#' @rdname loo.bef_fit
#' @exportS3Method loo::waic
waic.bef_fit <- function(x, ...) {
  .bef_require_loo()
  .bef_require_inherits(x, "bef_fit", "x", "bef_invalid_fit")
  if (.bef_is_v01_fit(x)) {
    .bef_abort_v01_fit("waic")
  }
  loo::waic(.bef_log_lik_draws(x), ...)
}

# The chain of each draw, in the order of the rows of the draws matrix.
.bef_draw_chain_id <- function(draws) {
  dims <- dim(draws)
  rep(seq_len(dims[[2L]]), each = dims[[1L]])
}

.bef_require_loo <- function() {
  if (!requireNamespace("loo", quietly = TRUE)) {
    .bef_abort_invalid_args(
      "The `loo` package is required for model comparison. Install it with install.packages(\"loo\").",
      arg = "loo"
    )
  }
  invisible(TRUE)
}
