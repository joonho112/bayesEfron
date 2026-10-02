#' @order 1
#' @describeIn bef_fit Prints the number of sites, the grid method, the running time and
#'   the worst value of each convergence diagnostic.
#' @export
print.bef_fit <- function(x, ...) {
  x <- .bef_prepare_fit(x)
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @order 2
#' @describeIn bef_fit Returns a list of class `summary.bef_fit` with `prior_summary`,
#'   the posterior means of the mean, variance and standard deviation of
#'   \eqn{g}, and `diagnostics`, the convergence diagnostics.
#' @export
summary.bef_fit <- function(object, level = 0.9, ...) {
  level <- .bef_validate_summary_level(level)
  object <- .bef_prepare_fit(object)
  metadata <- object$metadata

  out <- list(
    prior_summary = list(
      mean = metadata$mean_g_summary$mean,
      var = metadata$var_g_summary$mean,
      sd = .bef_metadata_attr(metadata, "sd_g_summary")$mean
    ),
    diagnostics = list(
      rhat = .bef_metadata_attr(metadata, "diagnostics")$rhat,
      ess_bulk = .bef_metadata_attr(metadata, "diagnostics")$ess_bulk,
      ess_tail = .bef_metadata_attr(metadata, "diagnostics")$ess_tail,
      divergences = .bef_metadata_attr(metadata, "diagnostics")$divergences,
      max_treedepth = .bef_metadata_attr(metadata, "diagnostics")$max_treedepth,
      effective_params = metadata$effective_params_summary,
      log_marginal_likelihood = metadata$log_marginal_likelihood_summary,
      model_family = metadata$model_family,
      stan_file_sha256 = metadata$stan_file_sha256,
      runtime_seconds = metadata$runtime_seconds,
      diagnostic_skipped = .bef_metadata_attr(metadata, "diagnostic_skipped"),
      sampler_diagnostics_failed = .bef_metadata_attr(
        metadata,
        "sampler_diagnostics_failed"
      ),
      sampler_diagnostics_warned = {
        warned <- .bef_metadata_attr(metadata, "sampler_diagnostics_warned")
        if (is.null(warned)) character() else warned
      }
    )
  )
  attr(out, "level") <- level
  attr(out, "summary_definition_version") <- 1L
  class(out) <- "summary.bef_fit"
  out
}

#' @order 4
#' @describeIn bef_fit Prints the summary.
#' @export
print.summary.bef_fit <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @order 3
#' @describeIn bef_fit Adds `theta_summary` to the summary, with the columns
#'   `hpdi_lower` and `hpdi_upper` computed for `level`.
#' @export
summary.bef_fit_re <- function(object, level = 0.9, ...) {
  object <- .bef_prepare_fit(object)
  out <- NextMethod()
  out$theta_summary <- .bef_theta_summary_for_level(object, level = attr(out, "level"))
  class(out) <- c("summary.bef_fit_re", "summary.bef_fit")
  out
}

#' @order 5
#' @describeIn bef_fit Returns the posterior means or the posterior modes of the
#'   site effects, the column `mean` or `map` of `theta_summary`, as a
#'   numeric vector named by site number.
#' @export
coef.bef_fit_re <- function(object, type = c("mean", "map"), ...) {
  type <- .bef_validate_method_choice(
    type, c("mean", "map"), arg = "type"
  )
  object <- .bef_prepare_fit(object)
  theta_summary <- object$metadata$theta_summary
  out <- theta_summary[[type]]
  names(out) <- as.character(theta_summary$site)
  out
}

#' @order 7
#' @describeIn bef_fit Returns a diagonal matrix with the posterior variances of the
#'   site effects, the squares of the column `sd` of `theta_summary`. The
#'   site effects are not independent in the posterior distribution,
#'   because all of them depend on \eqn{g}, and their covariances are not
#'   computed. The matrix therefore serves for one site at a time and not
#'   for the variance of a sum or a difference of site effects.
#' @export
vcov.bef_fit_re <- function(object, ...) {
  object <- .bef_prepare_fit(object)
  theta_summary <- object$metadata$theta_summary
  out <- diag(theta_summary$sd^2, nrow = nrow(theta_summary))
  dimnames(out) <- list(
    as.character(theta_summary$site),
    as.character(theta_summary$site)
  )
  out
}

#' @order 6
#' @describeIn bef_fit Returns a data frame with the columns `site`, `lower`, `upper`
#'   and `point`. For `type = "theta"` the limits are equal-tailed quantiles
#'   of the `theta_rep` draws and `point` is the posterior mean. For
#'   `type = "g"` the rows are the mean, variance and standard deviation of
#'   \eqn{g}.
#' @export
confint.bef_fit_re <- function(object,
                               parm = NULL,
                               level = 0.9,
                               type = c("theta", "g"),
                               ...) {
  level <- .bef_validate_summary_level(level)
  type <- .bef_validate_method_choice(
    type, c("theta", "g"), arg = "type"
  )
  object <- .bef_prepare_fit(object)
  if (identical(type, "theta")) {
    return(.bef_confint_theta(object, parm = parm, level = level))
  }
  .bef_confint_g(object, parm = parm, level = level)
}

#' @order 9
#' @describeIn bef_fit Returns `theta_summary`.
#' @export
as.data.frame.bef_fit_re <- function(x,
                                     row.names = NULL,
                                     optional = FALSE,
                                     ...) {
  x <- .bef_prepare_fit(x)
  out <- x$metadata$theta_summary
  if (!is.null(row.names)) {
    if (!is.character(row.names) || length(row.names) != nrow(out) || anyNA(row.names)) {
      .bef_abort_invalid_args(
        "`row.names` must be NULL or a non-missing character vector with one value per site.",
        arg = "row.names"
      )
    }
    row.names(out) <- row.names
  }
  out
}

#' @order 6
#' @describeIn as_bef_data Prints the number of sites and the range of the estimates and
#'   of the standard errors.
#' @export
print.bef_data <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @order 7
#' @describeIn as_bef_data Returns a list with the number of sites `K` and summaries of
#'   `theta_hat` and `sigma`.
#' @export
summary.bef_data <- function(object, ...) {
  validate_bef_data(object)
  list(
    K = length(object$theta_hat),
    theta_hat = .bef_numeric_summary(object$theta_hat),
    sigma = .bef_numeric_summary(object$sigma),
    source = object$source,
    names = object$names
  )
}

#' @order 4
#' @describeIn diagnose Prints the worst value of each diagnostic.
#' @export
print.bef_diagnostic <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @order 5
#' @describeIn diagnose Returns a list in which `rhat`, `ess_bulk` and `ess_tail` are
#'   reduced to their worst value, with the name of the quantity that has
#'   it.
#' @export
summary.bef_diagnostic <- function(object, ...) {
  validate_bef_diagnostic(object)
  list(
    rhat = .bef_diagnostic_extreme(object$rhat, "max"),
    ess_bulk = .bef_diagnostic_extreme(object$ess_bulk, "min"),
    ess_tail = .bef_diagnostic_extreme(object$ess_tail, "min"),
    divergences = .bef_diagnostic_sum(object$divergences),
    max_treedepth = .bef_diagnostic_sum(object$max_treedepth),
    ebfmi = object$ebfmi,
    effective_params = object$effective_params_summary,
    model_family = object$model_family,
    stan_file_sha256 = object$stan_file_sha256,
    runtime_seconds = object$runtime_seconds,
    diagnostic_skipped = object$diagnostic_skipped,
    sampler_diagnostics_failed = object$sampler_diagnostics_failed,
    sampler_diagnostics_warned = object$sampler_diagnostics_warned
  )
}

#' @order 10
#' @describeIn bef_fit Returns the number of sites.
#' @export
nobs.bef_fit <- function(object, ...) {
  as.integer(object$metadata$data_list$K)
}

#' @order 11
#' @describeIn bef_fit Returns the posterior mean of `log_marginal_likelihood`,
#'   the log likelihood with the site effects summed out, as a `logLik`
#'   object whose degrees of freedom are the posterior mean of
#'   `effective_params`. It is not a maximized log likelihood, so `AIC()`
#'   and `BIC()` computed from it do not have their usual meaning; to
#'   compare fits, use [loo.bef_fit()].
#' @export
logLik.bef_fit <- function(object, ...) {
  object <- .bef_prepare_fit(object)
  value <- object$metadata$log_marginal_likelihood_summary$mean
  structure(
    value,
    df = object$metadata$effective_params_summary$mean,
    nobs = stats::nobs(object),
    class = "logLik"
  )
}

#' @order 12
#' @describeIn bef_fit Returns the draws as a `posterior::draws_array`.
#' @exportS3Method posterior::as_draws
as_draws.bef_fit <- function(x, ...) {
  x <- .bef_prepare_fit(x)
  posterior::as_draws_array(x$draws)
}

.bef_validate_summary_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1) {
    .bef_abort_invalid_args(
      "`level` must be a finite numeric scalar between 0 and 1.",
      arg = "level"
    )
  }
  level
}

.bef_metadata_attr <- function(metadata, name) {
  attr(metadata, name, exact = TRUE)
}

.bef_numeric_summary <- function(x) {
  list(
    min = min(x),
    median = stats::median(x),
    mean = mean(x),
    max = max(x),
    sd = stats::sd(x)
  )
}

.bef_diagnostic_extreme <- function(x, direction) {
  present <- which(!is.na(x))
  if (length(present) == 0L) {
    return(list(value = NA_real_, index = NA_integer_, variable = NA_character_))
  }
  local_index <- switch(
    direction,
    max = which.max(x[present]),
    min = which.min(x[present]),
    .bef_abort_invalid_args(
      "`direction` must be \"max\" or \"min\".",
      arg = "direction"
    )
  )
  index <- present[[local_index]]
  # The name of the variable says which parameter the extreme value belongs to.
  variable <- if (!is.null(names(x))) names(x)[[index]] else NA_character_
  list(value = x[[index]], index = index, variable = variable)
}

.bef_diagnostic_sum <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  sum(x, na.rm = TRUE)
}

.bef_validate_method_choice <- function(x, choices, arg) {
  if (identical(x, choices)) {
    return(choices[[1L]])
  }
  if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% choices) {
    .bef_abort_invalid_args(
      sprintf(
        "`%s` must be one of: %s.",
        arg,
        paste(sprintf("\"%s\"", choices), collapse = ", ")
      ),
      arg = arg
    )
  }
  x
}

.bef_theta_summary_for_level <- function(object, level) {
  theta_summary <- object$metadata$theta_summary
  if (isTRUE(all.equal(level, 0.9))) {
    return(theta_summary)
  }
  intervals <- confint(object, level = level, type = "theta")
  theta_summary$hpdi_lower <- intervals$lower
  theta_summary$hpdi_upper <- intervals$upper
  theta_summary
}

.bef_confint_theta <- function(object, parm, level) {
  theta_rep <- .bef_site_draws(object, "theta_rep")
  theta_summary <- object$metadata$theta_summary
  selected <- .bef_resolve_site_parm(parm, theta_summary$site)
  probs <- .bef_interval_probs(level)
  intervals <- vapply(
    selected,
    function(site) {
      posterior::quantile2(
        theta_rep[, site],
        probs = probs,
        names = FALSE
      )
    },
    numeric(2L)
  )

  data.frame(
    site = theta_summary$site[selected],
    lower = intervals[1L, ],
    upper = intervals[2L, ],
    point = theta_summary$mean[selected],
    row.names = NULL,
    check.names = FALSE
  )
}

.bef_confint_g <- function(object, parm, level) {
  parameters <- c("mean_g", "var_g", "sd_g")
  selected <- .bef_resolve_g_parm(parm, parameters)
  probs <- .bef_interval_probs(level)
  intervals <- vapply(
    parameters[selected],
    function(parameter) {
      posterior::quantile2(
        object$posterior[[parameter]],
        probs = probs,
        names = FALSE
      )
    },
    numeric(2L)
  )
  points <- vapply(
    parameters[selected],
    function(parameter) {
      switch(
        parameter,
        mean_g = object$metadata$mean_g_summary$mean,
        var_g = object$metadata$var_g_summary$mean,
        sd_g = .bef_metadata_attr(object$metadata, "sd_g_summary")$mean
      )
    },
    numeric(1L)
  )

  data.frame(
    site = parameters[selected],
    lower = intervals[1L, ],
    upper = intervals[2L, ],
    point = points,
    row.names = NULL,
    check.names = FALSE
  )
}

.bef_interval_probs <- function(level) {
  alpha <- (1 - level) / 2
  c(alpha, 1 - alpha)
}

.bef_resolve_site_parm <- function(parm, sites) {
  K <- length(sites)
  if (is.null(parm)) {
    return(seq_len(K))
  }
  if (is.numeric(parm)) {
    if (!.bef_is_whole_number_vector(parm) || any(parm < 1L) || any(parm > K)) {
      .bef_abort_invalid_args(
        "`parm` must contain valid site indices.",
        arg = "parm"
      )
    }
    return(as.integer(parm))
  }
  if (is.character(parm) && !anyNA(parm)) {
    matched <- match(parm, as.character(sites))
    if (anyNA(matched)) {
      .bef_abort_invalid_args(
        "`parm` must contain site labels present in the fit.",
        arg = "parm"
      )
    }
    return(matched)
  }
  .bef_abort_invalid_args(
    "`parm` must be NULL, numeric site indices, or character site labels.",
    arg = "parm"
  )
}

.bef_resolve_g_parm <- function(parm, parameters) {
  if (is.null(parm)) {
    return(seq_along(parameters))
  }
  if (is.numeric(parm)) {
    if (!.bef_is_whole_number_vector(parm) ||
        any(parm < 1L) ||
        any(parm > length(parameters))) {
      .bef_abort_invalid_args(
        "`parm` must contain valid prior-summary indices.",
        arg = "parm"
      )
    }
    return(as.integer(parm))
  }
  if (is.character(parm) && !anyNA(parm)) {
    matched <- match(parm, parameters)
    if (anyNA(matched)) {
      .bef_abort_invalid_args(
        "`parm` must contain one or more of mean_g, var_g, sd_g.",
        arg = "parm"
      )
    }
    return(matched)
  }
  .bef_abort_invalid_args(
    "`parm` must be NULL, numeric prior-summary indices, or character prior-summary labels.",
    arg = "parm"
  )
}

.bef_is_whole_number_vector <- function(x) {
  is.numeric(x) && all(is.finite(x)) && all(x == as.integer(x))
}
