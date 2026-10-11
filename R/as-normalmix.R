#' Export the estimated distribution of effects for the ebnm package
#'
#' @description
#' `as_normalmix()` returns the posterior mean of \eqn{g} as a `normalmix`
#' object, the class that the ashr and ebnm packages use for a prior
#' distribution that is a mixture of normal distributions. The means of
#' the mixture components are the grid points. The ebnm package can use
#' this object to compute the posterior distributions of the effects with
#' \eqn{g} fixed at that estimate, so the fit can be compared with
#' empirical Bayes fits of the same data.
#'
#' @details
#' The mixture has one component per grid point. The mixture proportions
#' are the posterior means of the probabilities of \eqn{g}. The component
#' standard deviations are zero for a fit with point masses. For a fit
#' with `component = "normal"`, they all equal the standard deviation of
#' its normal components.
#'
#' Pass the object to `ebnm::ebnm_npmle()` as `g_init` with `fix_g = TRUE`
#' to obtain the plug-in posterior distributions of the effects, in which
#' \eqn{g} is treated as known. They differ from the posterior
#' distributions of the fit, which average over the posterior draws of
#' \eqn{g}. For most sites the posterior means are close and the standard
#' deviation from the fit is slightly larger, but neither is true for
#' every site. See the Details of [site_summary()] for the relation
#' between the two standard deviations.
#'
#' `as_normalmix()` constructs the object itself, so neither ashr nor ebnm
#' needs to be installed.
#'
#' @param x A `bef_fit` object created by [bayes_efron_fit()].
#' @param ... These dots are for future extensions and must be empty.
#'
#' @return An object of class `normalmix`, as created by
#'   `ashr::normalmix()`. It is a list with the following components, each
#'   with one value per grid point.
#'   \describe{
#'     \item{`pi`}{The mixture proportions.}
#'     \item{`mean`}{The component means, which are the grid points.}
#'     \item{`sd`}{The component standard deviations.}
#'   }
#'
#' @seealso [g_summary()], [site_summary()], `vignette("empirical-bayes")`
#'
#' @examples
#' g <- as_normalmix(krw_fit)
#' class(g)
#' # The means of the components are the grid points.
#' range(g$mean)
#'
#' @examplesIf requireNamespace("ebnm", quietly = TRUE)
#' # The plug-in posterior distributions, with g fixed at this estimate.
#' plug_in <- ebnm::ebnm_npmle(krw_discrimination$yi, krw_discrimination$sei,
#'                             g_init = g, fix_g = TRUE)
#' sites <- site_summary(krw_fit)
#' # The posterior means are close for most firms: their absolute
#' # differences, in units of the posterior standard deviation.
#' summary(abs(coef(plug_in) - sites$mean) / sites$sd)
#' # The standard deviations from the fit are larger for most firms, but
#' # not for all.
#' summary(sites$sd / fitted(plug_in)$sd)
#'
#' @export
as_normalmix <- function(x, ...) {
  .bef_check_dots_empty(list(...), "as_normalmix")
  .bef_require_inherits(x, "bef_fit", "x", "bef_invalid_fit")
  x <- .bef_prepare_fit(x)
  engine <- .bef_engine(x)
  probability <- colMeans(exp(engine$log_g))
  # A mixture as ashr::normalmix() builds it: a data frame of the three
  # columns with the class in place of that of a data frame.
  structure(
    data.frame(pi = unname(probability / sum(probability)), mean = engine$grid,
               sd = rep(engine$tau_c, engine$L)),
    class = "normalmix"
  )
}
