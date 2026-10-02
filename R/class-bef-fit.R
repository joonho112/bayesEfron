# The bef_data, bef_fit and bef_fit_re classes: constructors and validators.

#' Fitted Efron log-spline model objects
#'
#' @description
#' [bayes_efron_fit()] returns an object of class `bef_fit_re`, which
#' inherits from `bef_fit`. This page describes what the object holds and
#' the methods that summarize it. [plot.bef_fit_re()],
#' [predict.bef_fit_re()], [diagnose()] and [loo.bef_fit()] have pages of
#' their own.
#'
#' @details
#' The object is a list with three elements.
#'
#' `draws` is a `posterior::draws_array` with the posterior draws of every
#' quantity in the Stan model: the spline coefficients `alpha` and the
#' precision `lambda`; `log_g`, the log probabilities of \eqn{g} on the
#' grid; `mean_g`, `var_g` and `sd_g`, the mean, variance and standard
#' deviation of \eqn{g}; for each site the mean `theta_mean`, standard
#' deviation `theta_sd` and mode `theta_map` of the posterior distribution
#' of \eqn{\theta_i} given \eqn{g}, a draw `theta_rep` from that
#' distribution and the log likelihood `log_lik`; `effective_params`; and
#' `log_marginal_likelihood`, the sum of `log_lik` over the sites.
#'
#' `posterior` is a list with the draws of `mean_g`, `var_g`, `sd_g`,
#' `effective_params` and `log_marginal_likelihood` as numeric vectors.
#'
#' `metadata` is a list with the settings of the fit and summaries of the
#' draws:
#'
#' * `model_family` and `grid_method`, as given to [bayes_efron_fit()].
#' * `seed`, the seed that the sampler used; `cmdstan_version`; and
#'   `stan_file_sha256`, a hash of the Stan file.
#' * `data_list`, the data passed to Stan: the estimates `theta_hat`, the
#'   standard errors `sigma`, the grid points `grid`, the spline basis `B`,
#'   and the numbers of sites, grid points and spline functions, `K`, `L`
#'   and `M`.
#' * `runtime_seconds`, the time that sampling took.
#' * `theta_summary`, a data frame with one row for each site, described
#'   below.
#' * `mean_g_summary`, `var_g_summary`, `effective_params_summary` and
#'   `log_marginal_likelihood_summary`, each a list with the mean, the
#'   standard deviation and the 5%, 50% and 95% quantiles of the draws.
#'
#' The settings of the sampler are attributes of `metadata`.
#' `attr(fit$metadata, "sampler_settings")` is a list with the values of
#' `chains`, `parallel_chains`, `iter_warmup`, `iter_sampling`,
#' `adapt_delta`, `max_treedepth`, `seed`, `refresh` and `init` that were
#' passed to CmdStan, and the attributes `package_version` and
#' `cmdstanr_version` give the versions of bayesEfron and cmdstanr that
#' made the fit. A fit saved by an earlier version of the package does not
#' have these attributes.
#'
#' If the model was fitted with `keep_cmdstan_fit = TRUE` there is a fourth
#' element, `cmdstan_fit`.
#'
#' In `theta_summary`, `site` numbers the sites in the order of the data.
#' The other columns summarize the posterior distribution of
#' \eqn{\theta_i}, which averages over the posterior draws of \eqn{g}
#' (Lee and Sui, 2025, Section 4.3). `mean` is its mean and `sd` its
#' standard deviation. The posterior variance is the variance of
#' \eqn{\theta_i} given \eqn{g}, averaged over the draws of \eqn{g}, plus
#' the variance over those draws of the mean of \eqn{\theta_i} given
#' \eqn{g}, so that `sd` includes the uncertainty about \eqn{g}.
#' `hpdi_lower` and `hpdi_upper` are the 5% and 95% quantiles of the
#' posterior draws `theta_rep`; in spite of their names they are the limits
#' of an equal-tailed interval and not of a highest posterior density
#' interval. The draws take values on the grid, but a sample quantile is
#' interpolated between two adjacent values of the ordered draws, so a
#' limit can lie between two grid points. `map` is the posterior mode, the
#' grid point to which the posterior distribution gives the largest
#' probability.
#'
#' `effective_params` is the sum over the sites of the ratio of the
#' posterior variance of \eqn{\theta_i} given \eqn{g} to
#' \eqn{\sigma_i^2}. In the normal random-effects model that ratio is the
#' weight given to the site's own estimate, so the sum is small when the
#' estimates are shrunk strongly and near the number of sites when they
#' are hardly shrunk. It is a different quantity from the `p_loo` of
#' [loo.bef_fit()], which the loo package also calls an effective number
#' of parameters.
#'
#' \eqn{g} is the prior distribution of each \eqn{\theta_i}, and the
#' printed summary labels its mean, variance and standard deviation
#' "Prior g".
#'
#' The definitions of `sd`, `map` and `effective_params` have changed
#' between versions of the package, and a fit records which definitions its
#' summaries follow, in the attribute `summary_definition_version` of
#' `metadata`. For a fit that was saved by an earlier version, the methods
#' on this page compute `sd` and `map` again from the stored draws, and
#' `effective_params` too if it was stored under its first definition. They
#' give a warning that names the quantities that changed and leave the
#' saved object as it is. The result of `summary()` saved by an earlier
#' version holds no draws and cannot be brought up to date; `print()` stops
#' and asks for `summary()` to be called again on the fit. [diagnose()],
#' [predict.bef_fit_re()] and [loo.bef_fit()] need quantities that the
#' first released version did not store, and stop for a fit saved by it.
#'
#' @param x,object A `bef_fit` object from [bayes_efron_fit()], or for the
#'   `summary.bef_fit` methods the result of `summary()`.
#' @param ... For the `print()` methods, arguments passed on to `format()`.
#'   The other methods do not use them.
#' @param level Number between 0 and 1, the probability content of the
#'   intervals. Defaults to `0.9`.
#' @param type For `coef()`, `"mean"` (the default) for the posterior means
#'   or `"map"` for the posterior modes. For `confint()`, `"theta"` (the
#'   default) for intervals for the site effects or `"g"` for intervals for
#'   the mean, variance and standard deviation of \eqn{g}.
#' @param parm The rows to return: site numbers for `type = "theta"`, and
#'   for `type = "g"` one or more of `"mean_g"`, `"var_g"` and `"sd_g"` or
#'   their positions. `NULL` (the default) returns all.
#' @param row.names Character vector of row names, one for each site, or
#'   `NULL` (the default).
#' @param optional Not used; it is an argument of the generic.
#' @param use_cli Logical or `NULL`. With `TRUE` and the cli package
#'   installed, the first line is set in bold and section titles in
#'   color; `FALSE` gives plain text. `NULL` (the default) takes the value
#'   of `getOption("bayesEfron.use_cli", TRUE)`.
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
#' fit <- raudenbush_fit
#' fit
#' summary(fit)
#'
#' # Posterior means, next to the estimates they are shrunk from.
#' cbind(estimate = raudenbush1985$yi, posterior_mean = coef(fit))
#'
#' # Intervals for three sites, and for the mean and spread of g.
#' confint(fit, parm = c(4, 10, 18), level = 0.95)
#' confint(fit, type = "g")
#'
#' # Posterior standard deviations and posterior modes, with the standard
#' # errors of the estimates.
#' head(cbind(se = raudenbush1985$sei, as.data.frame(fit)[c("mean", "sd", "map")]))
#'
#' # Posterior variances as a diagonal matrix.
#' diag(vcov(fit))[1:3]
#'
#' @name bef_fit
#' @aliases bef_fit_re
NULL

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
      length(x$theta_hat) < 5L ||
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
  .bef_validate_site_draws_present(x$draws, K, class)
  .bef_validate_theta_summary_consistency(
    x$metadata$theta_summary, .bef_site_draws(x, "theta_mean"),
    .bef_site_draws(x, "theta_sd"), x$metadata$data_list$grid, class
  )

  x
}

.bef_fit_new <- function(draws,
                         metadata,
                         posterior = list(),
                         cmdstan_fit = NULL,
                         class = "bef_fit") {
  if (!is.null(metadata)) {
    attr(metadata, "summary_definition_version") <- 1L
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
