#' @details
#' To get started, read `vignette("bayesEfron")`, which walks through a
#' complete analysis of the callback gaps of 97 firms. The model, its prior
#' distributions and the way it is fitted are described in
#' [bayes_efron_fit()].
#'
#' @section Functions:
#' * [as_bef_data()] extracts the estimates and their standard errors from
#'   a list, a data frame or a `metafor::escalc()` object and checks them.
#' * [make_efron_grid()] builds the grid and the spline basis, and
#'   [plot()][plot.bef_grid] shows them together with the estimates.
#' * [bayes_efron_fit()] fits the model and returns a [bef_fit] object,
#'   which has `summary()`, `coef()`, `confint()` and other methods.
#' * [diagnose()] returns sampler and grid diagnostics.
#' * [g_summary()] summarizes the estimated distribution of effects, and
#'   [site_summary()] and [exceedance_prob()] summarize the site effects.
#' * [plot()][plot.bef_fit_re] plots a fit. The functions that it calls can
#'   also be called directly: [plot_g()], [plot_sites()],
#'   [plot_shrinkage()], [plot_exceedance()], [plot_site_posterior()],
#'   [plot_diagnostics()], [plot_loo()], [plot_ppc()], [plot_mcmc()] and,
#'   for simulated data, [plot_recovery()].
#' * [predict()][predict.bef_fit_re] returns draws of the effect or of the
#'   estimate of a new site.
#' * [loo()][loo.bef_fit] estimates how well a fit predicts the estimate of
#'   a site that is left out, and [compare_efron_grids()] compares fits on
#'   several grids by that measure.
#' * [as_normalmix()] returns the estimate of \eqn{g} as a `normalmix`
#'   object for the ebnm package.
#' * [bayes_efron_compile()] compiles the two Stan programs before the first
#'   fit, and [bayes_efron_clear_cache()] deletes the compiled models.
#'
#' @section Data sets and vignettes:
#' The following data sets come with a stored fit, which can be used
#' without CmdStan.
#'
#' * [krw_discrimination] and [krw_fit]: the callback gaps of 97 firms,
#'   analyzed in `vignette("bayesEfron")`.
#' * [raudenbush1985] and [raudenbush_fit]: a meta-analysis of 19 studies,
#'   analyzed in `vignette("teacher-expectancy")`.
#' * [twin_towers] and [twin_towers_fit]: a simulated data set with known
#'   true effects, analyzed in `vignette("twin-towers")`.
#'
#' The data set [star_schools] contains estimates of the effect of small
#' classes in kindergarten in 79 schools of Project STAR, a multisite
#' trial.
#'
#' The other vignettes cover data input, the choice of the grid,
#' diagnostics, plots, a comparison with empirical Bayes methods, the
#' model, and the implementation. `browseVignettes("bayesEfron")` lists
#' them all.
#'
#' @section Funding:
#' This research was supported by the Institute of Education Sciences,
#' U.S. Department of Education, through Grant R305D240078 to the
#' University of Alabama. The opinions expressed are those of the
#' authors and do not represent views of the Institute or the U.S.
#' Department of Education.
#'
#' @docType package
#' @name bayesEfron-package
#' @aliases bayesEfron
#' @keywords internal
#' @importFrom graphics plot
#' @importFrom stats AIC BIC coef confint logLik nobs quantile residuals simulate vcov
#' @importFrom rlang .data
"_PACKAGE"

# `density` is a variable that ggplot2 computes inside after_stat().
if (getRversion() >= "2.15.1") {
  utils::globalVariables("density")
}
