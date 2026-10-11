#' @order 1
#' @describeIn bef_fit Prints the number of sites, the number of draws,
#'   the worst value of each sampler diagnostic, the range of the standard
#'   errors, the grid and the method that chose its range, and the
#'   functions to call next.
#' @export
print.bef_fit <- function(x, ...) {
  x <- .bef_prepare_fit(x)
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @order 2
#' @describeIn bef_fit Returns a list of class `summary.bef_fit` with the
#'   components `prior_summary` (the posterior means of the mean, variance
#'   and standard deviation of \eqn{g}) and `diagnostics` (the sampler
#'   diagnostics; see [diagnose()]).
#' @export
summary.bef_fit <- function(object, level = 0.9, ...) {
  .bef_check_dots_empty(list(...), "summary")
  level <- .bef_validate_summary_level(level)
  object <- .bef_prepare_fit(object)
  metadata <- object$metadata
  flags <- .bef_fit_diagnostic_flags(object)

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
      diagnostic_skipped = flags$skipped,
      sampler_diagnostics_failed = flags$failed,
      sampler_diagnostics_warned = flags$warned
    )
  )
  attr(out, "level") <- level
  attr(out, "summary_definition_version") <- .bef_summary_version()
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
#' @describeIn bef_fit Adds the component `theta_summary` to the summary:
#'   the table of site effects described in Details, with equal-tailed
#'   intervals at `level`.
#' @export
summary.bef_fit_re <- function(object, level = 0.9, ...) {
  object <- .bef_prepare_fit(object)
  out <- NextMethod()
  out$theta_summary <- .bef_site_table_at_level(object, level = attr(out, "level"))
  .bef_note_site_table_names()
  class(out) <- c("summary.bef_fit_re", "summary.bef_fit")
  out
}

#' @order 5
#' @describeIn bef_fit Returns the posterior means or the posterior modes
#'   of the site effects (the column `mean` or `mode` of `theta_summary`)
#'   as a numeric vector. Its names are the site labels, or the site
#'   numbers if the data had no labels.
#' @export
coef.bef_fit_re <- function(object, type = c("mean", "mode", "map"), ...) {
  type <- .bef_validate_method_choice(
    type, c("mean", "mode", "map"), arg = "type"
  )
  object <- .bef_prepare_fit(object)
  table <- .bef_site_table(object)
  out <- table[[type]]
  names(out) <- table$label
  out
}

#' @order 7
#' @describeIn bef_fit Returns the posterior variances of the site effects
#'   (the squares of the column `sd` of `theta_summary`) as a diagonal
#'   matrix. If `full = TRUE`, the matrix also contains the posterior
#'   covariances between sites. The site effects are independent given
#'   \eqn{g} but not in the posterior distribution, because all of them
#'   depend on \eqn{g}. Their covariances are the covariances, over the
#'   draws of \eqn{g}, of the means of the effects given \eqn{g}. The
#'   diagonal matrix is sufficient when sites are considered one at a
#'   time. Use `full = TRUE` for the variance of a sum or a difference of
#'   site effects.
#' @param full If `TRUE`, `vcov()` also returns the covariances between
#'   sites. If `FALSE` (the default), the matrix is diagonal.
#' @export
vcov.bef_fit_re <- function(object, full = FALSE, ...) {
  if (!is.logical(full) || length(full) != 1L || is.na(full)) {
    .bef_abort_invalid_args("`full` must be TRUE or FALSE.", arg = "full")
  }
  object <- .bef_prepare_fit(object)
  table <- .bef_site_table(object)
  if (!full) {
    out <- diag(table$sd^2, nrow = nrow(table))
    dimnames(out) <- list(table$label, table$label)
    return(out)
  }
  K <- nrow(table)
  if (K >= 1000L) {
    message(sprintf(
      "The full covariance matrix has %d rows and columns and takes about %.0f MB.",
      K, 8 * K^2 / 2^20
    ))
  }
  conditional <- .bef_site_conditional(.bef_engine(object))
  centered <- sweep(conditional$mean, 2L, colMeans(conditional$mean))
  out <- crossprod(centered) / nrow(centered)
  diag(out) <- diag(out) + colMeans(conditional$var)
  dimnames(out) <- list(table$label, table$label)
  out
}

#' @order 13
#' @describeIn bef_fit Returns quantiles of the posterior distributions of
#'   the site effects as a matrix with one row per site and one column per
#'   probability in `probs`.
#' @param probs Numeric vector of probabilities for `quantile()`, each
#'   greater than 0 and less than 1.
#' @param reading How `quantile()` computes quantiles from point masses on
#'   the grid. `"interpolated"` (the default) spreads the probability of
#'   each grid point evenly over one grid spacing around it, and `"step"`
#'   leaves it on the grid point. See the Details of [g_summary()].
#' @export
quantile.bef_fit_re <- function(x,
                                probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                                reading = c("interpolated", "step"),
                                ...) {
  probs <- .bef_validate_probs(probs)
  reading <- .bef_validate_method_choice(
    reading, c("interpolated", "step"), arg = "reading"
  )
  x <- .bef_prepare_fit(x)
  table <- .bef_site_table(x)
  out <- .bef_site_quantile(.bef_engine(x), probs, reading = reading)
  dimnames(out) <- list(table$label, paste0(format(100 * probs, trim = TRUE), "%"))
  out
}

#' @order 14
#' @describeIn bef_fit Returns the estimates minus the posterior means of
#'   the site effects, as a numeric vector named by the site labels.
#' @export
residuals.bef_fit_re <- function(object, ...) {
  object <- .bef_prepare_fit(object)
  table <- .bef_site_table(object)
  stats::setNames(table$estimate - table$mean, table$label)
}

#' @order 15
#' @describeIn bef_fit Draws the site effects jointly from their posterior
#'   distribution. For each of the `nsim` draws, one posterior draw of
#'   \eqn{g} is taken and every site effect is drawn from its distribution
#'   given that \eqn{g}. Returns a matrix with one row per draw and one
#'   column per site. Each row is one joint draw, so functions of several
#'   sites (e.g., ranks or the largest effect) can be computed row by row.
#' @param nsim Number of joint draws for `simulate()`, a positive integer.
#' @param seed Integer seed for the draws of `simulate()`, or `NULL` (the
#'   default) for none. If a seed is given, the state of R's random number
#'   generator is restored afterward.
#' @export
simulate.bef_fit_re <- function(object, nsim = 1, seed = NULL, ...) {
  if (!.bef_is_whole_number(nsim) || nsim < 1L) {
    .bef_abort_invalid_args("`nsim` must be a positive whole number.", arg = "nsim")
  }
  object <- .bef_prepare_fit(object)
  if (!is.null(seed)) {
    if (!.bef_is_whole_number(seed)) {
      .bef_abort_invalid_args("`seed` must be NULL or a whole number.", arg = "seed")
    }
    # Leave R's random number generator as it was found.
    had_state <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    old <- if (had_state) get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL
    set.seed(as.integer(seed))
    on.exit({
      if (had_state) {
        assign(".Random.seed", old, .GlobalEnv)
      } else {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
  }
  table <- .bef_site_table(object)
  engine <- .bef_engine(object)
  out <- .bef_simulate_sites(engine, as.integer(nsim))
  colnames(out) <- table$label
  out
}

#' @order 6
#' @describeIn bef_fit Returns a data frame of posterior intervals with the
#'   attribute `level`. For `type = "theta"` it has one row per site and
#'   the columns `site`, `mean`, `lower` and `upper`. For `type = "g"` it
#'   has one row each for the mean, variance and standard deviation of
#'   \eqn{g} and the columns `parameter`, `mean`, `lower` and `upper`. The
#'   limits are computed as specified by `method` and `interval`. The
#'   column `point` and, for `type = "g"`, the column `site` repeat `mean`
#'   and `parameter`. They are deprecated and will be removed.
#' @export
confint.bef_fit_re <- function(object,
                               parm = NULL,
                               level = 0.9,
                               type = c("theta", "g"),
                               method = c("exact", "draws"),
                               interval = c("equal-tailed", "hpd"),
                               ...) {
  level <- .bef_validate_summary_level(level)
  type <- .bef_validate_method_choice(
    type, c("theta", "g"), arg = "type"
  )
  method <- .bef_validate_method_choice(
    method, c("exact", "draws"), arg = "method"
  )
  interval <- .bef_validate_method_choice(
    interval, c("equal-tailed", "hpd"), arg = "interval"
  )
  object <- .bef_prepare_fit(object)
  out <- if (identical(type, "theta")) {
    .bef_confint_theta(object, parm = parm, level = level, method = method,
                       interval = interval)
  } else {
    .bef_confint_g(object, parm = parm, level = level, interval = interval)
  }
  attr(out, "level") <- level
  attr(out, "interval") <- interval
  attr(out, "method") <- if (identical(type, "theta")) method else "draws"
  class(out) <- c("bef_intervals", "data.frame")
  out
}

# Printing the result of confint(): the level and the kind of interval,
# then the data frame.
#' @exportS3Method print bef_intervals
#' @noRd
print.bef_intervals <- function(x, ...) {
  level <- attr(x, "level", exact = TRUE)
  if (!is.null(level)) {
    cat(sprintf(
      "%s%% posterior intervals (%s)\n", format(100 * level),
      if (identical(attr(x, "interval", exact = TRUE), "hpd")) {
        "shortest"
      } else {
        "equal-tailed"
      }
    ))
  }
  print(as.data.frame(x), ...)
  invisible(x)
}

#' @order 9
#' @describeIn bef_fit Returns `theta_summary`, the table of site effects
#'   described in Details, with 90% intervals.
#' @export
as.data.frame.bef_fit_re <- function(x,
                                     row.names = NULL,
                                     optional = FALSE,
                                     ...) {
  x <- .bef_prepare_fit(x)
  out <- .bef_site_table(x)
  .bef_note_site_table_names()
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
#' @describeIn as_bef_data Prints the number of sites and the ranges of the
#'   estimates and of the standard errors.
#' @export
print.bef_data <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @order 7
#' @describeIn as_bef_data Returns a list with the number of sites `K` and
#'   summaries of `theta_hat` and `sigma`.
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

#' @order 8
#' @describeIn as_bef_data Returns a data frame with one row per site and
#'   the columns `site`, `label`, `estimate` and `se`. If there are no site
#'   labels, `label` is the site number.
#' @param row.names,optional Not used; they are arguments of the generic
#'   function.
#' @export
as.data.frame.bef_data <- function(x, row.names = NULL, optional = FALSE, ...) {
  validate_bef_data(x)
  K <- length(x$theta_hat)
  data.frame(
    site = seq_len(K),
    label = .bef_site_labels(x$names, K),
    estimate = x$theta_hat,
    se = x$sigma,
    stringsAsFactors = FALSE
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
#' @describeIn diagnose Returns a list in which `rhat`, `ess_bulk` and
#'   `ess_tail` are reduced to their worst value (the largest R-hat and
#'   the smallest effective sample sizes), with the name of the variable
#'   that has it.
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
    sampler_diagnostics_warned = object$sampler_diagnostics_warned,
    grid = object$grid
  )
}

#' @order 10
#' @describeIn bef_fit Returns the number of sites.
#' @export
nobs.bef_fit <- function(object, ...) {
  as.integer(object$metadata$data_list$K)
}

#' @order 11
#' @describeIn bef_fit Returns the posterior mean of
#'   `log_marginal_likelihood`, the marginal log likelihood with the site
#'   effects integrated out, as a `logLik` object. Its degrees of freedom
#'   are the posterior mean of `effective_params`. Note that this is not a
#'   maximized log likelihood, so `AIC()` and `BIC()` throw an error. Use
#'   [loo.bef_fit()] to compare fits.
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
#' @describeIn bef_fit Returns the posterior draws as a
#'   `posterior::draws_array`.
#' @exportS3Method posterior::as_draws
as_draws.bef_fit <- function(x, ...) {
  x <- .bef_prepare_fit(x)
  posterior::as_draws_array(x$draws)
}

#' @order 16
#' @describeIn bef_fit Returns the number of posterior draws stored in the
#'   fit, summed over the chains.
#' @exportS3Method posterior::ndraws
ndraws.bef_fit <- function(x) {
  as.integer(prod(dim(x$draws)[seq_len(2L)]))
}

#' @order 17
#' @describeIn bef_fit Returns the names of the variables in the draws.
#' @exportS3Method posterior::variables
variables.bef_fit <- function(x, ...) {
  dimnames(x$draws)[[3L]]
}

#' @order 18
#' @describeIn bef_fit `AIC()` and `BIC()` throw an error, for the reason
#'   given under `logLik()`.
#' @param k Not used.
#' @export
AIC.bef_fit <- function(object, ..., k = 2) {
  .bef_abort_no_information_criterion("AIC")
}

#' @order 19
#' @rdname bef_fit
#' @export
BIC.bef_fit <- function(object, ...) {
  .bef_abort_no_information_criterion("BIC")
}

.bef_abort_no_information_criterion <- function(name) {
  .bef_abort_invalid_args(
    sprintf(
      paste0(
        "`%s()` is not defined for this model: `logLik()` returns a posterior ",
        "mean and not a maximized log likelihood. Use `loo::loo()` to compare fits."
      ),
      name
    ),
    arg = "object"
  )
}

# An argument that a method does not have is an error and not ignored, so
# that a misspelled name does not pass unnoticed.
.bef_check_dots_empty <- function(dots, what) {
  if (length(dots) == 0L) {
    return(invisible(TRUE))
  }
  named <- names(dots)
  named <- named[nzchar(named)]
  .bef_abort_invalid_args(
    sprintf(
      "`%s()` has no argument%s %s.", what, if (length(dots) == 1L) "" else "s",
      if (length(named) > 0L) {
        paste(sprintf("`%s`", named), collapse = ", ")
      } else {
        "without a name"
      }
    ),
    arg = "..."
  )
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

.bef_confint_theta <- function(object, parm, level, method = "exact",
                               interval = "equal-tailed") {
  theta_summary <- .bef_site_table(object)
  selected <- .bef_resolve_site_parm(parm, theta_summary$label)
  probs <- .bef_interval_probs(level)
  if (identical(method, "draws")) {
    theta_rep <- .bef_theta_rep_draws(object)
    limits <- t(vapply(
      selected,
      function(site) {
        if (identical(interval, "hpd")) {
          .bef_draws_hpd(theta_rep[, site], level)
        } else {
          posterior::quantile2(theta_rep[, site], probs = probs, names = FALSE)
        }
      },
      numeric(2L)
    ))
  } else {
    engine <- .bef_engine(object)
    weights <- .bef_site_weights(engine, sites = selected)
    limits <- if (identical(interval, "hpd")) {
      .bef_site_hpd(engine, level, weights = weights, sites = selected)
    } else {
      .bef_site_quantile(engine, probs, weights = weights, sites = selected)
    }
  }

  data.frame(
    site = theta_summary$site[selected],
    mean = theta_summary$mean[selected],
    lower = limits[, 1L],
    upper = limits[, 2L],
    point = theta_summary$mean[selected],
    row.names = if (.bef_has_site_labels(theta_summary)) theta_summary$label[selected] else NULL,
    check.names = FALSE
  )
}

# The theta_rep draws of a fit that holds them, for `method = "draws"`.
.bef_theta_rep_draws <- function(object) {
  K <- as.integer(object$metadata$data_list$K)
  wanted <- sprintf("theta_rep[%d]", seq_len(K))
  if (!all(wanted %in% dimnames(object$draws)[[3L]])) {
    .bef_abort_invalid_args(
      paste0(
        "`method = \"draws\"` needs the `theta_rep` draws, which this fit does ",
        "not hold. Use `method = \"exact\"`, or fit the model again with ",
        "`site_draws = TRUE`."
      ),
      arg = "method"
    )
  }
  .bef_site_draws(object, "theta_rep")
}

.bef_confint_g <- function(object, parm, level, interval = "equal-tailed") {
  parameters <- c("mean_g", "var_g", "sd_g")
  selected <- .bef_resolve_g_parm(parm, parameters)
  # The same summaries as g_summary(type = "moments").
  moments <- .bef_g_moments(object, level = level)[selected, , drop = FALSE]
  intervals <- rbind(moments$lower, moments$upper)
  if (identical(interval, "hpd")) {
    intervals <- vapply(
      parameters[selected],
      function(parameter) .bef_draws_hpd(object$posterior[[parameter]], level),
      numeric(2L)
    )
  }
  points <- moments$mean

  data.frame(
    parameter = parameters[selected],
    mean = points,
    lower = intervals[1L, ],
    upper = intervals[2L, ],
    site = parameters[selected],
    point = points,
    row.names = NULL,
    check.names = FALSE
  )
}

.bef_interval_probs <- function(level) {
  alpha <- (1 - level) / 2
  c(alpha, 1 - alpha)
}

# The sites that `parm` selects: numbers are positions in the data, and
# character strings are site labels. Without labels the label of a site is
# its number as text, so the two agree.
.bef_resolve_site_parm <- function(parm, labels) {
  K <- length(labels)
  if (is.null(parm)) {
    return(seq_len(K))
  }
  if (is.numeric(parm)) {
    if (!.bef_is_whole_number_vector(parm) || any(parm < 1L) || any(parm > K)) {
      .bef_abort_invalid_args(
        sprintf("Numbers in `parm` are positions of sites and must lie between 1 and %d.", K),
        arg = "parm"
      )
    }
    return(as.integer(parm))
  }
  if (is.character(parm) && !anyNA(parm)) {
    matched <- match(parm, as.character(labels))
    if (anyNA(matched)) {
      .bef_abort_invalid_args(
        sprintf(
          "`parm` names site labels that the fit does not have: %s.",
          paste(sprintf("\"%s\"", parm[is.na(matched)]), collapse = ", ")
        ),
        arg = "parm"
      )
    }
    return(matched)
  }
  .bef_abort_invalid_args(
    "`parm` must be NULL, numbers (positions of sites) or character strings (site labels).",
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
