#' @details
#' For sites \eqn{i = 1, \ldots, K} the model is
#' \deqn{\hat\theta_i \mid \theta_i \sim N(\theta_i, \sigma_i^2),
#'   \qquad \theta_i \sim g,}
#' where the standard errors \eqn{\sigma_i} are taken as known and may
#' differ between sites. The distribution \eqn{g} is represented on a
#' grid, with log probabilities that are a combination of natural cubic
#' spline functions; the spline coefficients have a normal prior whose
#' precision has a half-Cauchy prior. [bayes_efron_fit()] gives the
#' details.
#'
#' @section Main functions:
#' * [as_bef_data()] takes the estimates and standard errors out of a list,
#'   a data frame or a `metafor::escalc()` object and checks them.
#' * [bayes_efron_fit()] fits the model and returns a [bef_fit] object,
#'   which has `summary()`, `coef()`, `confint()` and other methods.
#' * [diagnose()] returns the convergence diagnostics of a fit.
#' * [plot()][plot.bef_fit_re] draws the site effects, the estimated
#'   distribution of effects, the diagnostics or the Pareto k values.
#' * [predict()][predict.bef_fit_re] draws the effect or the estimate of a
#'   new site.
#' * [make_efron_grid()] builds the grid and the spline basis, and
#'   [compare_efron_grids()] compares several of them by leave-one-out
#'   cross-validation.
#' * [bayes_efron_compile()] compiles the Stan model before the first fit,
#'   and [bayes_efron_clear_cache()] removes the compiled model.
#'
#' The examples use the data set [raudenbush1985] and the fit
#' [raudenbush_fit]. `vignette("bayesEfron")` is the place to start.
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
#' @importFrom stats coef confint logLik nobs vcov
"_PACKAGE"

if (getRversion() >= "2.15.1") {
  utils::globalVariables(".data")
}
