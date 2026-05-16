#' bayesEfron: Fully Bayesian Inference for the Empirical-Bayes Deconvolution Problem
#'
#' @description
#' `bayesEfron` fits the fully Bayesian Efron log-spline prior for
#' univariate random-effects meta-analytic deconvolution with
#' heteroscedastic within-study standard errors. The package is for
#' the applied meta-analyst who has a vector of study-level effect
#' estimates and their standard errors and wants two related
#' quantities: posterior summaries of the latent site effects, and a
#' continuous, smoothly regularised estimate of the underlying mixing
#' distribution that generated those effects. Both objects are
#' returned from a single fit, with calibrated credible intervals
#' derived from the same posterior draws.
#'
#' @details
#' For sites \eqn{i = 1, \ldots, K},
#' \deqn{\hat\theta_i \mid \theta_i \sim \mathcal{N}(\theta_i, \sigma_i^2),
#'   \qquad \theta_i \mid g \sim g,}
#' where \eqn{g} is the unknown mixing distribution. The package
#' represents \eqn{g} as a discrete distribution on a fixed grid of
#' length \eqn{L} (default \eqn{L = 101}) whose log-density is a
#' linear combination of \eqn{M} natural-cubic-spline basis functions
#' (default \eqn{M = 6}). The spline coefficients receive a weakly
#' informative Gaussian prior whose precision \eqn{\lambda} is in turn
#' assigned a half-Cauchy hyperprior, so the smoothness of the
#' deconvolved density is itself part of what the sampler learns.
#'
#' Posterior computation is delegated to CmdStan via the `cmdstanr`
#' package, with two-tier compile caching so that repeated fits with
#' the same model code reuse the compiled binary. The within-study
#' standard errors enter exactly: each site's likelihood contribution
#' uses its own \eqn{\sigma_i}, with no homoscedasticity assumption.
#'
#' @section User-facing entry points:
#' \describe{
#'   \item{[bayes_efron_fit()]}{Primary fitting entry. Takes per-site
#'     effect estimates and standard errors, runs the eight-stage
#'     pipeline (validation, grid construction, Stan-data preparation,
#'     cache-backed model retrieval, sampling, draw extraction,
#'     postprocessing, assembly), and returns a fitted `bef_fit_re`
#'     object carrying posterior draws, the deconvolved density, and
#'     posterior summaries of the site effects.}
#'   \item{[make_efron_grid()]}{Constructs the discrete support of
#'     \eqn{g} and the natural-cubic-spline basis evaluated on it.
#'     Four recipes: paper real-data, paper simulation, paper
#'     sensitivity, and an experimental KL-target recipe.}
#'   \item{[as_bef_data()]}{Input adapter. Converts a plain list of
#'     `theta_hat`/`sigma`, an [metafor::escalc()] object, or an
#'     existing `bef_data` object into the canonical input class.}
#'   \item{[diagnose()]}{Diagnostic producer. Returns a
#'     `bef_diagnostic` object summarising R-hat, effective sample
#'     sizes, divergent transitions, max-treedepth saturations, and
#'     other sampler-health quantities.}
#' }
#'
#' @section Where to start:
#' New users should begin with the applied vignette \emph{A1 \enc{·}{.}
#' Getting started}, which walks through a minimal end-to-end fit on
#' a five-site toy dataset. The companion methodological vignette
#' \emph{M1 \enc{·}{.} The empirical-Bayes deconvolution problem}
#' develops the statistical background, motivates the log-spline
#' prior, and explains how the deconvolved density \eqn{\hat g}
#' relates to the posterior site effects \eqn{\theta_i}.
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
