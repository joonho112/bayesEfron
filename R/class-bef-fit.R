# The bef_data, bef_fit and bef_fit_re classes: constructors and validators.

#' Fitted Efron log-spline model objects
#'
#' @description
#' Objects of class `c("bef_fit_re", "bef_fit")`, as returned by
#' [bayes_efron_fit()], contain the posterior draws, the data, the settings
#' of the fit and summaries of the draws (see Details). There are methods
#' for printing and summarizing a fit, for posterior means, intervals and
#' quantiles of the site effects, and for extracting the posterior draws.
#' For plots, predictions, diagnostics and model comparison, see
#' [plot.bef_fit_re()], [predict.bef_fit_re()], [diagnose()] and
#' [loo.bef_fit()].
#'
#' @details
#' A fitted object is a list with the components `draws`, `posterior` and
#' `metadata`, and with `cmdstan_fit` if the model was fitted with
#' `keep_cmdstan_fit = TRUE`.
#'
#' `draws` is a `posterior::draws_array` with the posterior draws of the
#' following variables.
#'
#' * `lp__`: The log density of the posterior distribution up to a
#'   constant, as reported by CmdStan.
#' * `alpha`: The spline coefficients.
#' * `lambda`: The precision of the prior for the spline coefficients.
#' * `log_g`: The log probabilities of \eqn{g} at the grid points.
#' * `mean_g`, `var_g`, `sd_g`: The mean, variance and standard deviation
#'   of \eqn{g}.
#' * `effective_params`: The effective number of parameters defined below.
#' * `log_marginal_likelihood`: The marginal log likelihood of the
#'   estimates, with the site effects integrated out.
#'
#' A fit made with `site_draws = TRUE` stores five more variables per site.
#'
#' * `theta_mean`, `theta_sd`, `theta_map`: The mean, standard deviation
#'   and mode of the posterior distribution of \eqn{\theta_i}{theta_i}
#'   given \eqn{g}.
#' * `theta_rep`: A draw from that distribution.
#' * `log_lik`: The log likelihood of the site's estimate.
#'
#' Otherwise these variables are not stored, and the methods compute the
#' site summaries from the draws of \eqn{g}.
#'
#' `posterior` is a list with the draws of `mean_g`, `var_g`, `sd_g`,
#' `effective_params` and `log_marginal_likelihood` as numeric vectors.
#'
#' `metadata` is a list with the data, the settings of the fit and
#' summaries of the draws. It has the following components.
#'
#' * `model_family`: The model family, currently always `"RE"`.
#' * `grid_method`: The method that chose the limits of the grid, or
#'   `"support"` if the limits were specified directly.
#' * `seed`: The seed used by the sampler.
#' * `cmdstan_version`: The version of CmdStan.
#' * `stan_file_sha256`: The SHA-256 hash of the Stan file.
#' * `data_list`: The data passed to Stan, a list with the estimates
#'   `theta_hat`, the standard errors `sigma`, the grid points `grid`, the
#'   spline basis `B`, and the numbers of sites, grid points and spline
#'   functions (`K`, `L` and `M`).
#' * `runtime_seconds`: The time taken by sampling, in seconds.
#' * `theta_summary`: The table of site effects, a data frame with one row
#'   per site. Its columns are described below.
#' * `mean_g_summary`, `var_g_summary`, `effective_params_summary`,
#'   `log_marginal_likelihood_summary`: Lists with the mean, the standard
#'   deviation and the 5%, 50% and 95% quantiles of the draws.
#'
#' More information is stored in attributes of `metadata` and can be
#' extracted with `attr()`, e.g., `attr(fit$metadata, "grid")`.
#'
#' * `sampler_settings`: A list with the values of `chains`,
#'   `parallel_chains`, `iter_warmup`, `iter_sampling`, `adapt_delta`,
#'   `max_treedepth`, `seed`, `refresh` and `init` passed to CmdStan.
#' * `grid`: The grid of the fit, as returned by [make_efron_grid()]. It
#'   can be supplied as `grid` in another fit.
#' * `storage`: A list recording whether the site-level draws are stored,
#'   the number of chains, the number of iterations stored per chain, and
#'   the value of `thin`.
#' * `package_version`, `cmdstanr_version`: The versions of bayesEfron and
#'   cmdstanr used for the fit.
#'
#' `theta_summary`, the table of site effects, has the columns listed
#' below. The columns from `mean` to `mode` summarize the posterior
#' distribution of \eqn{\theta_i}{theta_i}, which averages over the
#' posterior draws of \eqn{g} (Lee and Sui, 2025, Section 4.3).
#'
#' * `site`: The site number, in the order of the data.
#' * `label`: The site label, or the site number if the data had no labels.
#' * `estimate`, `se`: The estimate and its standard error, as supplied.
#' * `mean`: The posterior mean.
#' * `sd`: The posterior standard deviation. The posterior variance is the
#'   variance of \eqn{\theta_i}{theta_i} given \eqn{g}, averaged over the
#'   draws of \eqn{g}, plus the variance over those draws of the mean of
#'   \eqn{\theta_i}{theta_i} given \eqn{g}, so `sd` includes the
#'   uncertainty about \eqn{g}.
#' * `lower`, `upper`: The 5% and 95% quantiles of the posterior
#'   distribution, i.e., the limits of an equal-tailed 90% interval. They
#'   are quantiles of the distribution itself and not sample quantiles of
#'   site-level draws. For a fit with point masses they are computed with
#'   the `"interpolated"` reading (described under `reading`), so they
#'   generally fall between grid points.
#' * `mode`: The posterior mode. For point masses this is the grid point
#'   with the largest posterior probability, and for normal components the
#'   value at which the posterior density is largest.
#' * `hpdi_lower`, `hpdi_upper`, `map`: Deprecated columns that repeat
#'   `lower`, `upper` and `mode`. They will be removed.
#'
#' Like all site summaries, the posterior summaries in this table are
#' computed from the stored draws of \eqn{g} and are subject to the Monte
#' Carlo error of those draws. [site_summary()] returns the same table for
#' any level, with additional columns.
#'
#' `effective_params` is the sum over sites of the ratio of the posterior
#' variance of \eqn{\theta_i}{theta_i} given \eqn{g} to
#' \eqn{\sigma_i^2}{sigma_i^2}. In the normal random-effects model this
#' ratio is the weight given to the site's own estimate, so the sum is
#' small when the estimates are shrunk strongly and close to the number of
#' sites when there is little shrinkage. Note that `p_loo` in
#' [loo.bef_fit()], which the loo package also calls the effective number
#' of parameters, is a different quantity.
#'
#' In the printed summary, the mean, variance and standard deviation of
#' \eqn{g} are labeled "Prior g", because \eqn{g} is the prior distribution
#' of each \eqn{\theta_i}{theta_i}.
#'
#' A fit saved by an earlier version of the package may contain a table of
#' site effects or a value of `effective_params` computed under other
#' definitions. The methods then recompute them from the stored draws of
#' \eqn{g} and issue a warning that names the recomputed quantities. The
#' saved object is not changed.
#'
#' @param x,object A `bef_fit` object created by [bayes_efron_fit()]. For
#'   the `summary.bef_fit` methods, the result of `summary()`.
#' @param ... Arguments passed on to `format()` by the `print()` methods.
#'   `summary()` throws an error if any are supplied; the other methods
#'   ignore them.
#' @param level Posterior probability of the intervals, a number between 0
#'   and 1.
#' @param type For `coef()`, the summary to return: `"mean"` (the default)
#'   for the posterior means or `"mode"` for the posterior modes (`"map"`,
#'   the earlier name of `"mode"`, is also accepted). For `confint()`, the
#'   quantities for which intervals are returned: `"theta"` (the default)
#'   for the site effects or `"g"` for the mean, variance and standard
#'   deviation of \eqn{g}.
#' @param parm The sites or parameters for which `confint()` returns
#'   intervals. For `type = "theta"`, a numeric vector of site positions in
#'   the data or a character vector of site labels. For `type = "g"`, one
#'   or more of `"mean_g"`, `"var_g"` and `"sd_g"`, or their positions. If
#'   `NULL` (the default), all are returned.
#' @param method How `confint()` computes the limits for `type = "theta"`.
#'   `"exact"` (the default) takes them from the posterior distribution of
#'   the site effect, which averages over the draws of \eqn{g}. `"draws"`
#'   takes them as sample quantiles of the `theta_rep` draws and requires a
#'   fit made with `site_draws = TRUE`.
#' @param interval Type of interval returned by `confint()`.
#'   `"equal-tailed"` (the default) has its limits at the quantiles
#'   `(1 - level) / 2` and `1 - (1 - level) / 2`. `"hpd"` is the shortest
#'   interval with posterior probability `level`; if the distribution has
#'   one mode, it is the highest posterior density interval. For
#'   `method = "draws"` and for `type = "g"`, `"hpd"` gives the shortest
#'   interval that contains at least that proportion of the draws. The two
#'   types differ whenever the posterior distribution is not symmetric,
#'   and most when it has two modes.
#' @param row.names Optional character vector of row names for
#'   `as.data.frame()`, one per site.
#' @param optional Not used; it is an argument of the generic.
#' @param use_cli If `TRUE` and the cli package is installed, the first
#'   line of the output is bold and, in a printed summary, the section
#'   titles are colored. If `FALSE`, the output is plain text. If `NULL`
#'   (the default), the
#'   value of `getOption("bayesEfron.use_cli", TRUE)` is used.
#'
#' @return `summary()` returns a list of class `summary.bef_fit`. `coef()`,
#'   `fitted()` and `residuals()` return numeric vectors with one element
#'   per site. `confint()` and `as.data.frame()` return data frames, and
#'   `vcov()`, `quantile()` and `simulate()` return matrices. The `print()`
#'   methods return `x` invisibly. The other methods are described in the
#'   two lists below: "Methods (by generic)" for class `bef_fit`, and
#'   "Functions" for class `bef_fit_re` and for a summary.
#'
#' @references
#' Lee, J. and Sui, D. (2025). Fully Bayesian inference for meta-analytic
#' deconvolution using Efron's log-spline prior. *Mathematics*, 13(16),
#' 2639. \doi{10.3390/math13162639}
#'
#' @seealso [bayes_efron_fit()], [diagnose()], [plot.bef_fit_re()],
#'   [predict.bef_fit_re()], [loo.bef_fit()]
#'
#' @examples
#' fit <- krw_fit
#' fit
#' summary(fit)$prior_summary
#'
#' # Posterior means, next to the estimates they are shrunk from.
#' head(cbind(estimate = krw_discrimination$yi, posterior_mean = coef(fit)))
#'
#' # Intervals for three sites, and for the mean and spread of g.
#' confint(fit, parm = c(29, 11, 65), level = 0.95)
#' confint(fit, type = "g")
#'
#' # Posterior standard deviations and posterior modes, with the standard
#' # errors of the estimates.
#' head(cbind(se = krw_discrimination$sei, as.data.frame(fit)[c("mean", "sd", "mode")]))
#'
#' # Posterior variances as a diagonal matrix.
#' diag(vcov(fit))[1:3]
#'
#' @name bef_fit
#' @aliases bef_fit_re
NULL

# The smallest number of sites of a data object and of a fit. With fewer the
# spline with its default six columns has more coefficients than there are
# estimates. make_efron_grid() alone accepts two, so that a grid can be
# looked at for any set of estimates.
.bef_min_sites <- function() {
  5L
}

new_bef_data <- function(theta_hat, sigma, names = NULL, source = "user") {
  structure(
    list(
      theta_hat = theta_hat,
      sigma = sigma,
      names = names,
      source = source
    ),
    class = "bef_data"
  )
}

validate_bef_data <- function(x) {
  class <- "bef_invalid_args"
  .bef_require_inherits(x, "bef_data", "x", class)
  .bef_require_fields(x, .bef_bef_data_fields(), "`bef_data`", class)

  if (!is.numeric(x$theta_hat) ||
      length(x$theta_hat) < .bef_min_sites() ||
      any(!is.finite(x$theta_hat))) {
    .bef_abort_validate(
      "`bef_data$theta_hat` must be a finite numeric vector of length at least 5.",
      class
    )
  }

  K <- length(x$theta_hat)
  if (!is.numeric(x$sigma) ||
      length(x$sigma) != K ||
      any(!is.finite(x$sigma)) ||
      any(x$sigma <= 0)) {
    .bef_abort_validate(
      "`bef_data$sigma` must be a strictly positive finite numeric vector with the same length as `theta_hat`.",
      class
    )
  }

  if (!is.null(x$names) &&
      (!is.character(x$names) || length(x$names) != K || anyNA(x$names))) {
    .bef_abort_validate(
      "`bef_data$names` must be NULL or a character vector with one non-missing label per site.",
      class
    )
  }
  if (!is.null(x$names) && anyDuplicated(x$names) > 0L) {
    repeated <- unique(x$names[duplicated(x$names)])
    .bef_abort_validate(
      sprintf(
        "The site labels must be distinct; repeated: %s.",
        paste(sprintf("\"%s\"", utils::head(repeated, 5L)), collapse = ", ")
      ),
      class,
      duplicated = repeated
    )
  }

  if (!.bef_is_string(x$source)) {
    .bef_abort_validate(
      "`bef_data$source` must be a non-empty character scalar.",
      class
    )
  }

  x
}

new_bef_fit <- function(draws, metadata, posterior = list(), cmdstan_fit = NULL) {
  .bef_fit_new(
    draws = draws,
    metadata = metadata,
    posterior = posterior,
    cmdstan_fit = cmdstan_fit,
    class = "bef_fit"
  )
}

validate_bef_fit <- function(x) {
  class <- "bef_invalid_fit"
  .bef_require_inherits(x, "bef_fit", "x", class)
  .bef_fit_summary_version(x)
  .bef_require_fields(x, .bef_fit_fields(), "`bef_fit`", class)

  .bef_validate_draws_array(x$draws, class)
  if (!is.list(x$metadata)) {
    .bef_abort_validate("`bef_fit$metadata` must be a list.", class)
  }
  if (!is.list(x$posterior)) {
    .bef_abort_validate("`bef_fit$posterior` must be a list.", class)
  }

  .bef_require_fields(
    x$metadata, .bef_universal_metadata_fields(), "`bef_fit$metadata`", class
  )
  .bef_validate_metadata_core(x$metadata, class)
  .bef_validate_stan_data_list(x$metadata$data_list, class)
  .bef_validate_summary_list(x$metadata$mean_g_summary, "mean_g_summary", class)
  .bef_validate_summary_list(x$metadata$var_g_summary, "var_g_summary", class)
  .bef_validate_summary_list(
    x$metadata$effective_params_summary,
    "effective_params_summary",
    class
  )
  .bef_validate_summary_list(
    x$metadata$log_marginal_likelihood_summary,
    "log_marginal_likelihood_summary",
    class
  )

  x
}

new_bef_fit_re <- function(draws,
                           metadata,
                           posterior = list(),
                           cmdstan_fit = NULL) {
  .bef_fit_new(
    draws = draws,
    metadata = metadata,
    posterior = posterior,
    cmdstan_fit = cmdstan_fit,
    class = c("bef_fit_re", "bef_fit")
  )
}

# The metadata fields of an object saved by bayesEfron 0.1. They are used only
# to recognize such an object: some methods still work on it, and the others
# say why they cannot.
.bef_v01_metadata_fields <- function() {
  c("model_family", "grid_method", "seed", "cmdstan_version",
    "stan_file_sha256", "data_list", "runtime_seconds", "mean_g_summary",
    "var_g_summary", "theta_summary", "theta_rep_draws",
    "effective_params_summary", "log_marginal_likelihood_summary")
}

# Was this object created by bayesEfron 0.1?
.bef_is_v01_fit <- function(x) {
  is.list(x) && is.list(x$metadata) &&
    setequal(names(x$metadata), .bef_v01_metadata_fields())
}

# The error for a method that cannot use an object created by bayesEfron 0.1.
.bef_abort_v01_fit <- function(what) {
  .bef_abort_invalid_fit(
    paste0(
      "This object was created by bayesEfron 0.1, and `", what, "()` needs ",
      "quantities that later versions store differently.\n",
      "These methods still work on it: summary(), coef(), confint(), ",
      "as.data.frame(), nobs(), print().\n",
      "Refit the model to use diagnose(), loo(), waic() and predict()."
    )
  )
}

validate_bef_fit_re <- function(x) {
  class <- "bef_invalid_fit"
  .bef_require_inherits(x, "bef_fit_re", "x", class)
  if (!identical(class(x), c("bef_fit_re", "bef_fit"))) {
    .bef_abort_validate(
      "`bef_fit_re` objects must have class vector c(\"bef_fit_re\", \"bef_fit\").",
      class
    )
  }

  validate_bef_fit(x)
  .bef_require_exact_fields(
    x$metadata, .bef_fit_re_metadata_fields(), "`bef_fit_re$metadata`", class
  )
  .bef_require_exact_fields(
    x$posterior, .bef_generated_quantity_fields(), "`bef_fit_re$posterior`",
    class
  )

  if (!identical(x$metadata$model_family, "RE")) {
    .bef_abort_validate(
      "`bef_fit_re$metadata$model_family` must be \"RE\".",
      class
    )
  }

  K <- as.integer(x$metadata$data_list$K)
  n_draws <- prod(dim(x$draws)[seq_len(2L)])
  .bef_validate_postprocess_metadata_attrs(x$metadata, class)
  .bef_validate_theta_summary(x$metadata$theta_summary, K, class)
  .bef_validate_generated_quantities(x$posterior, K, n_draws, class)
  .bef_validate_site_draws(x$draws, K, class)
  .bef_validate_theta_summary_consistency(x, class)

  x
}

.bef_fit_new <- function(draws,
                         metadata,
                         posterior = list(),
                         cmdstan_fit = NULL,
                         class = "bef_fit") {
  if (!is.null(metadata)) {
    attr(metadata, "summary_definition_version") <- .bef_summary_version()
  }
  out <- list(
    draws = draws,
    metadata = metadata,
    posterior = posterior
  )

  if (!is.null(cmdstan_fit)) {
    out$cmdstan_fit <- cmdstan_fit
  }

  structure(out, class = class)
}

.bef_abort_validate <- function(message, class, ...) {
  if (identical(class, "bef_invalid_args")) {
    .bef_abort_invalid_args(message, ..., validate = TRUE)
  } else if (identical(class, "bef_invalid_fit")) {
    .bef_abort_invalid_fit(message, ...)
  } else {
    .bef_abort(
      message,
      class,
      ...,
      extra_class = "bef_validate_error"
    )
  }
}

.bef_bef_data_fields <- function() {
  c("theta_hat", "sigma", "names", "source")
}

.bef_fit_fields <- function() {
  c("draws", "metadata", "posterior")
}

.bef_universal_metadata_fields <- function() {
  c(
    "model_family",
    "grid_method",
    "seed",
    "cmdstan_version",
    "stan_file_sha256",
    "data_list",
    "runtime_seconds",
    "mean_g_summary",
    "var_g_summary",
    "effective_params_summary",
    "log_marginal_likelihood_summary"
  )
}

.bef_fit_re_metadata_fields <- function() {
  c(
    "model_family",
    "grid_method",
    "seed",
    "cmdstan_version",
    "stan_file_sha256",
    "data_list",
    "runtime_seconds",
    "mean_g_summary",
    "var_g_summary",
    "theta_summary",
    "effective_params_summary",
    "log_marginal_likelihood_summary"
  )
}

# The scalar series kept in `fit$posterior`. The site-level matrices
# (theta_map, theta_mean, theta_sd, theta_rep) are in `fit$draws` and are read
# with .bef_site_draws().
.bef_generated_quantity_fields <- function() {
  c(
    "mean_g",
    "var_g",
    "sd_g",
    "effective_params",
    "log_marginal_likelihood"
  )
}

.bef_site_generated_quantity_fields <- function() {
  c("theta_map", "theta_mean", "theta_sd", "theta_rep")
}
