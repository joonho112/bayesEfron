#' S3 methods for bayesEfron objects
#'
#' @description
#' The package exports a coordinated set of S3 methods for the
#' `bef_fit`, `bef_fit_re`, `bef_data`, and `bef_diagnostic` classes.
#' Together they cover printing, summarising, extracting credible
#' intervals, pulling point estimates, converting to a data frame,
#' counting observations, computing the marginal log-likelihood, and
#' coercing to the `posterior::draws_array` format.
#'
#' This help page documents the family-agnostic surface, which works
#' on any fitted `bef_fit` object regardless of the model family. The
#' RE-specific child class adds `coef()`, `vcov()`, `confint()`,
#' `as.data.frame()`, and a refined `summary()`. Plotting is
#' documented separately at [plot.bef_fit_re()]; the diagnostic
#' producer is at [diagnose()].
#'
#' @details
#' # Method index
#'
#' | Method | Class dispatched on | Returns |
#' |:-------|:-------------------|:--------|
#' | `print()` | `bef_fit`, `bef_data`, `bef_diagnostic`, `summary.bef_fit` | invisibly the input |
#' | `summary()` | `bef_fit`, `bef_fit_re`, `bef_data`, `bef_diagnostic` | a list (class `summary.bef_fit*` for fits) |
#' | `format()` | `bef_fit`, `summary.bef_fit`, `bef_data`, `bef_diagnostic` | character vector |
#' | `coef()` | `bef_fit_re` | named numeric vector of site point estimates |
#' | `vcov()` | `bef_fit_re` | diagonal matrix of site posterior variances |
#' | `confint()` | `bef_fit_re` | data frame of credible intervals |
#' | `as.data.frame()` | `bef_fit_re` | the per-site `theta_summary` table |
#' | `nobs()` | `bef_fit` | integer site count `K` |
#' | `logLik()` | `bef_fit` | `logLik` object with `df` = effective parameters |
#' | `posterior::as_draws()` | `bef_fit` | `draws_array` |
#'
#' # Metadata access
#'
#' `fit$metadata` is a closed list of 13 named fields; four
#' additional payloads are stored as attributes
#' (`diagnostics`, `diagnostic_skipped`,
#' `sampler_diagnostics_failed`, `sd_g_summary`). The methods on this
#' page are the recommended access path: `summary()` surfaces the
#' user-relevant fields; [diagnose()] returns a structured
#' `bef_diagnostic` object covering the diagnostic attributes; direct
#' `attr()` access is supported but unnecessary.
#'
#' # Sites and labels
#'
#' Site labels propagated from the input (named list or `escalc` row
#' names) appear in `coef()` names, `vcov()` dimnames, the `site`
#' column of `confint()` and `as.data.frame()`, and `summary()` print
#' output. When labels are absent the numeric site index is used.
#'
#' # Print and format backends
#'
#' `print()` and `format()` for `bef_fit`, `summary.bef_fit`,
#' `bef_data`, and `bef_diagnostic` prefer `cli`-styled output when
#' the `cli` package is installed and stdout supports it, and fall
#' back to a plain base-R representation otherwise. The `use_cli`
#' argument lets callers force one branch:
#'
#' * `use_cli = NULL` (default) — auto-detect.
#' * `use_cli = TRUE` — force `cli` styling (errors if `cli` is not
#'   installed).
#' * `use_cli = FALSE` — force plain base output.
#'
#' Setting the environment variable `BAYESEFRON_NO_CLI=1` is
#' equivalent to `use_cli = FALSE` for every call in the session and
#' is the recommended way to suppress styling in CI logs and
#' redirected stdout.
#'
#' # Site count requirement
#'
#' The `bef_data` validator (used by `as_bef_data()` and
#' `bayes_efron_fit()`) requires at least five sites. The constraint
#' reflects the methodological requirement of the log-spline
#' deconvolution prior, not a software limitation.
#'
#' @param x,object A bayesEfron S3 object (`bef_fit`, `bef_fit_re`,
#'   `bef_data`, or `bef_diagnostic`, depending on the method).
#' @param ... Additional arguments passed to methods. Most methods
#'   ignore `...`; `format()`-based printers accept it for
#'   compatibility with the generic.
#' @param level Numeric credible level in \eqn{(0, 1)}. Defaults to
#'   `0.9` (90 percent credible interval).
#' @param type Character option for methods with multiple estimands:
#'   `coef()` accepts `"mean"` (posterior mean, default) or `"map"`
#'   (posterior mode); `confint()` accepts `"theta"` (per-site
#'   latent-effect intervals, default) or `"g"` (mixing-distribution
#'   functionals).
#' @param parm Optional parameter or site subset. `NULL` (default)
#'   returns all sites or all `g`-functionals. Numeric values are
#'   1-based site indices; character values are matched against the
#'   site labels.
#' @param row.names Optional character vector of row names for
#'   `as.data.frame()`, with one value per site.
#' @param optional Included for `as.data.frame()` method
#'   compatibility; not used.
#' @param use_cli `NULL`, `TRUE`, or `FALSE`; controls optional
#'   `cli`-styled output for `format()` and `print()` methods. See
#'   the "Print and format backends" subsection in \strong{Details}.
#'
#' @return The return type depends on the method; see the index table
#'   in \strong{Details}.
#'
#' @seealso
#'   * [bayes_efron_fit()] for producing the fitted object.
#'   * [diagnose()] for the structured diagnostic producer.
#'   * [plot.bef_fit_re()] for visualization.
#'
#' @examples
#' # Load the cached five-site smoke fit shipped with the package.
#' fit <- readRDS(system.file(
#'   "examples", "cached_fit_re_smoke.rds",
#'   package = "bayesEfron"
#' ))
#'
#' summary(fit)
#' print(fit)
#' coef(fit)                       # posterior mean per site
#' coef(fit, type = "map")         # posterior MAP per site
#' confint(fit)                    # 90 percent credible intervals on theta
#' confint(fit, level = 0.95, type = "g")
#' nobs(fit)                       # site count K
#' logLik(fit)
#' as.data.frame(fit)              # per-site summary as a data frame
#' posterior::as_draws(fit)        # draws_array for downstream tools
#'
#' @name bayesEfron-methods
NULL

#' @rdname bayesEfron-methods
#' @export
print.bef_fit <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname bayesEfron-methods
#' @export
summary.bef_fit <- function(object, level = 0.9, ...) {
  level <- .bef_validate_summary_level(level)
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
      )
    )
  )
  attr(out, "level") <- level
  class(out) <- "summary.bef_fit"
  out
}

#' @rdname bayesEfron-methods
#' @export
print.summary.bef_fit <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname bayesEfron-methods
#' @export
summary.bef_fit_re <- function(object, level = 0.9, ...) {
  out <- NextMethod()
  out$theta_summary <- .bef_theta_summary_for_level(object, level = attr(out, "level"))
  class(out) <- c("summary.bef_fit_re", "summary.bef_fit")
  out
}

#' @rdname bayesEfron-methods
#' @export
coef.bef_fit_re <- function(object, type = c("mean", "map"), ...) {
  type <- .bef_validate_method_choice(
    type, c("mean", "map"), arg = "type", module = "coef.bef_fit_re"
  )
  theta_summary <- object$metadata$theta_summary
  out <- theta_summary[[type]]
  names(out) <- as.character(theta_summary$site)
  out
}

#' @rdname bayesEfron-methods
#' @export
vcov.bef_fit_re <- function(object, ...) {
  theta_summary <- object$metadata$theta_summary
  out <- diag(theta_summary$sd^2, nrow = nrow(theta_summary))
  dimnames(out) <- list(
    as.character(theta_summary$site),
    as.character(theta_summary$site)
  )
  out
}

#' @rdname bayesEfron-methods
#' @export
confint.bef_fit_re <- function(object,
                               parm = NULL,
                               level = 0.9,
                               type = c("theta", "g"),
                               ...) {
  level <- .bef_validate_summary_level(level)
  type <- .bef_validate_method_choice(
    type, c("theta", "g"), arg = "type", module = "confint.bef_fit_re"
  )
  if (identical(type, "theta")) {
    return(.bef_confint_theta(object, parm = parm, level = level))
  }
  .bef_confint_g(object, parm = parm, level = level)
}

#' @rdname bayesEfron-methods
#' @export
as.data.frame.bef_fit_re <- function(x,
                                     row.names = NULL,
                                     optional = FALSE,
                                     ...) {
  out <- x$metadata$theta_summary
  if (!is.null(row.names)) {
    if (!is.character(row.names) || length(row.names) != nrow(out) || anyNA(row.names)) {
      .bef_abort_invalid_args(
        "`row.names` must be NULL or a non-missing character vector with one value per site.",
        arg = "row.names",
        module = "as.data.frame.bef_fit_re"
      )
    }
    row.names(out) <- row.names
  }
  out
}

#' @rdname bayesEfron-methods
#' @export
print.bef_data <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname bayesEfron-methods
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

#' @rdname bayesEfron-methods
#' @export
print.bef_diagnostic <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname bayesEfron-methods
#' @export
summary.bef_diagnostic <- function(object, ...) {
  validate_bef_diagnostic(object)
  list(
    rhat = .bef_diagnostic_extreme(object$rhat, "max"),
    ess_bulk = .bef_diagnostic_extreme(object$ess_bulk, "min"),
    ess_tail = .bef_diagnostic_extreme(object$ess_tail, "min"),
    divergences = .bef_diagnostic_sum(object$divergences),
    max_treedepth = .bef_diagnostic_sum(object$max_treedepth),
    effective_params = object$effective_params_summary,
    model_family = object$model_family,
    stan_file_sha256 = object$stan_file_sha256,
    runtime_seconds = object$runtime_seconds,
    diagnostic_skipped = object$diagnostic_skipped,
    sampler_diagnostics_failed = object$sampler_diagnostics_failed
  )
}

#' @rdname bayesEfron-methods
#' @export
nobs.bef_fit <- function(object, ...) {
  as.integer(object$metadata$data_list$K)
}

#' @rdname bayesEfron-methods
#' @export
logLik.bef_fit <- function(object, ...) {
  value <- object$metadata$log_marginal_likelihood_summary$mean
  structure(
    value,
    df = object$metadata$effective_params_summary$mean,
    nobs = stats::nobs(object),
    class = "logLik"
  )
}

#' @rdname bayesEfron-methods
#' @exportS3Method posterior::as_draws
as_draws.bef_fit <- function(x, ...) {
  posterior::as_draws_array(x$draws)
}

.bef_validate_summary_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1) {
    .bef_abort_invalid_args(
      "`level` must be a finite numeric scalar between 0 and 1.",
      arg = "level",
      module = "summary.bef_fit"
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
    return(list(value = NA_real_, index = NA_integer_))
  }
  local_index <- switch(
    direction,
    max = which.max(x[present]),
    min = which.min(x[present]),
    .bef_abort_invalid_args(
      "`direction` must be \"max\" or \"min\".",
      arg = "direction",
      module = ".bef_diagnostic_extreme"
    )
  )
  index <- present[[local_index]]
  list(value = x[[index]], index = index)
}

.bef_diagnostic_sum <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  sum(x, na.rm = TRUE)
}

.bef_validate_method_choice <- function(x, choices, arg, module) {
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
      arg = arg,
      predicate = paste(choices, collapse = "|"),
      module = module
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
  theta_summary <- object$metadata$theta_summary
  selected <- .bef_resolve_site_parm(parm, theta_summary$site)
  probs <- .bef_interval_probs(level)
  intervals <- vapply(
    selected,
    function(site) {
      posterior::quantile2(
        object$metadata$theta_rep_draws[, site],
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
        arg = "parm",
        module = "confint.bef_fit_re"
      )
    }
    return(as.integer(parm))
  }
  if (is.character(parm) && !anyNA(parm)) {
    matched <- match(parm, as.character(sites))
    if (anyNA(matched)) {
      .bef_abort_invalid_args(
        "`parm` must contain site labels present in the fit.",
        arg = "parm",
        module = "confint.bef_fit_re"
      )
    }
    return(matched)
  }
  .bef_abort_invalid_args(
    "`parm` must be NULL, numeric site indices, or character site labels.",
    arg = "parm",
    module = "confint.bef_fit_re"
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
        arg = "parm",
        module = "confint.bef_fit_re"
      )
    }
    return(as.integer(parm))
  }
  if (is.character(parm) && !anyNA(parm)) {
    matched <- match(parm, parameters)
    if (anyNA(matched)) {
      .bef_abort_invalid_args(
        "`parm` must contain one or more of mean_g, var_g, sd_g.",
        arg = "parm",
        module = "confint.bef_fit_re"
      )
    }
    return(matched)
  }
  .bef_abort_invalid_args(
    "`parm` must be NULL, numeric prior-summary indices, or character prior-summary labels.",
    arg = "parm",
    module = "confint.bef_fit_re"
  )
}

.bef_is_whole_number_vector <- function(x) {
  is.numeric(x) && all(is.finite(x)) && all(x == as.integer(x))
}
