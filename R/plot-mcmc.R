#' Trace plots and other MCMC plots with bayesplot
#'
#' @description
#' `plot_mcmc()` is a convenience wrapper around bayesplot functions for
#' the posterior draws of a fit: trace plots, rank plots and densities by
#' chain, pairs plots, and the energy plot of the sampler. For any other
#' plot, pass the draws (`x$draws`) and `bayesplot::nuts_params(x)` to
#' bayesplot directly.
#'
#' @details
#' Each type calls one bayesplot function: `"trace"` calls `mcmc_trace()`,
#' `"rank"` calls `mcmc_rank_overlay()`, `"dens"` calls
#' `mcmc_dens_overlay()`, `"pairs"` calls `mcmc_pairs()` (which also marks
#' the divergent transitions), and `"energy"` calls `mcmc_nuts_energy()`.
#' The bayesplot package must be installed, and the plots use its theme
#' and colors.
#'
#' Only the draws stored in the fit are plotted. If the fit was thinned
#' (`thin`), these are a subset of the iterations, whereas [diagnose()]
#' and [plot_diagnostics()] describe the whole run. The energy plot
#' compares the energy at each iteration with its change from the previous
#' iteration, so it needs every iteration. For a thinned fit, including
#' the fits stored in the package, `type = "energy"` therefore throws an
#' error. In that case [diagnose()] still reports the E-BFMI of the whole
#' run.
#'
#' In a trace plot, chains that have converged cover the same range and
#' none of them drifts. See `vignette("diagnostics")` and
#' `vignette("plots")` for how to read the trace plot and the rank plot.
#'
#' @inheritParams plot_g
#' @param type Type of plot: `"trace"` (the default), `"rank"`, `"dens"`,
#'   `"pairs"` or `"energy"`. See Details.
#' @param pars Character vector with the names of the parameters to plot.
#'   A name without an index selects all elements of that parameter, e.g.,
#'   `"alpha"` selects every spline coefficient. By default the spline
#'   coefficients, `lambda`, and the mean and standard deviation of \eqn{g}
#'   are plotted. Not used with `type = "energy"`.
#'
#' @return The object returned by the bayesplot function: a `ggplot`
#'   object or, for `type = "pairs"`, a grid of plots that is shown when it
#'   is printed.
#'
#' @seealso [plot_diagnostics()], [diagnose()], `vignette("diagnostics")`
#'
#' @examplesIf requireNamespace("bayesplot", quietly = TRUE)
#' plot_mcmc(krw_fit, pars = c("lambda", "mean_g", "sd_g"))
#' plot_mcmc(krw_fit, type = "rank", pars = "alpha")
#' plot_mcmc(krw_fit, type = "dens", pars = c("mean_g", "sd_g"))
#'
#' @export
plot_mcmc <- function(x,
                      type = c("trace", "rank", "dens", "pairs", "energy"),
                      pars = c("alpha", "lambda", "mean_g", "sd_g"),
                      ...) {
  .bef_check_dots_empty(list(...), "plot_mcmc")
  x <- .bef_plot_fit(x)
  type <- .bef_validate_method_choice(
    type, c("trace", "rank", "dens", "pairs", "energy"), arg = "type"
  )
  .bef_require_bayesplot()
  if (identical(type, "energy")) {
    # Between kept iterations that are `thin` apart the energy changes far
    # more than between neighbors, so the plot of a thinned fit would show a
    # sampler that moves freely whatever the run was like.
    thin <- .bef_metadata_attr(x$metadata, "storage")$thin
    if (isTRUE(thin > 1L)) {
      .bef_abort_invalid_args(
        paste0(
          "The energy plot compares each iteration with the one before it, and this fit ",
          "keeps one iteration in ", as.integer(thin), ". Fit the model without `thin` to ",
          "draw it; `diagnose()` gives the E-BFMI of the whole run."
        ),
        arg = "type"
      )
    }
    return(bayesplot::mcmc_nuts_energy(bayesplot::nuts_params(x)))
  }
  variables <- .bef_plot_mcmc_variables(pars, dimnames(x$draws)[[3L]])
  draws <- unclass(x$draws)[, , variables, drop = FALSE]
  names(dimnames(draws)) <- c("iterations", "chains", "parameters")
  dimnames(draws)[[2L]] <- sprintf("chain:%d", seq_len(dim(draws)[[2L]]))
  switch(
    type,
    trace = bayesplot::mcmc_trace(draws),
    rank = bayesplot::mcmc_rank_overlay(draws),
    dens = bayesplot::mcmc_dens_overlay(draws),
    pairs = bayesplot::mcmc_pairs(draws, np = bayesplot::nuts_params(x))
  )
}

# The variables that `pars` names: a name with an index is one variable, a
# name without stands for every element of that parameter.
.bef_plot_mcmc_variables <- function(pars, variables) {
  if (!is.character(pars) || length(pars) == 0L || anyNA(pars) || any(!nzchar(pars))) {
    .bef_abort_invalid_args("`pars` must be a character vector of parameter names.",
                            arg = "pars")
  }
  base <- sub("[[].*$", "", variables)
  found <- lapply(pars, function(name) variables[variables == name | base == name])
  missing <- pars[lengths(found) == 0L]
  if (length(missing) > 0L) {
    .bef_abort_invalid_args(
      sprintf("`pars` names no parameter of the fit: %s. The parameters are: %s.",
              paste(sprintf("\"%s\"", missing), collapse = ", "),
              paste(unique(base), collapse = ", ")),
      arg = "pars"
    )
  }
  unique(unlist(found))
}

.bef_require_bayesplot <- function() {
  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    .bef_abort_invalid_args(
      "The `bayesplot` package is required for `plot_mcmc()`. Install it with install.packages(\"bayesplot\").",
      arg = "type"
    )
  }
  invisible(TRUE)
}
