#' MCMC convergence and grid diagnostics for a fitted model
#'
#' @description
#' `diagnose()` returns the MCMC diagnostics of a fit (R-hat, effective
#' sample sizes, divergent transitions, tree depth and E-BFMI) and flags
#' those that indicate a problem. The diagnostics were computed by
#' [bayes_efron_fit()] when the model was fitted, so no sampling is done.
#' `diagnose()` also checks the grid, in particular whether the result for
#' any site depends on where the grid ends.
#'
#' @details
#' A diagnostic fails its check and is named in `sampler_diagnostics_failed`
#' in two cases.
#'
#' * An R-hat value exceeds 1.05.
#' * There is a divergent transition.
#'
#' A diagnostic leads to a warning and is named in
#' `sampler_diagnostics_warned` in the following cases.
#'
#' * An R-hat value is between 1.01 and 1.05.
#' * A bulk or tail effective sample size is below 400.
#' * An iteration reached the maximum tree depth.
#' * The energy Bayesian fraction of missing information (E-BFMI) of a
#'   chain is below 0.3.
#'
#' R-hat below 1.01 and effective sample sizes above 400 are the
#' recommendations of Vehtari et al. (2021), and 0.3 is the value below
#' which CmdStan warns about E-BFMI. The other thresholds are choices of
#' this package. See `vignette("diagnostics")` for what to do about each
#' kind of warning.
#'
#' The sampler diagnostics describe the whole run, even if the fit was
#' thinned. If the run has fewer than 400 draws in total, R-hat and the
#' effective sample sizes are too unreliable to act on. Their values are
#' still returned, but they are not compared with the thresholds, and
#' `diagnostic_skipped` contains `"rhat_check"`, `"ess_bulk_check"` and
#' `"ess_tail_check"`. Divergent transitions, the maximum tree depth and
#' E-BFMI are checked regardless of the length of the run.
#'
#' R-hat and the effective sample sizes are `NA` for a quantity that has
#' the same value in every draw, e.g., the posterior mode of a site that
#' never leaves one grid point in a fit with `site_draws = TRUE`. If a
#' quantity is constant within each chain but differs between chains, the
#' chains have not mixed and its R-hat is infinite. That R-hat is stored as
#' `NA`, and `"rhat"` is named in `sampler_diagnostics_failed` even for a
#' run of fewer than 400 draws. If R-hat or an effective sample size
#' cannot be computed at all, or an effective sample size is infinite, all
#' values of that diagnostic are `NA` and it is named in
#' `diagnostic_skipped`.
#'
#' The `grid` component contains the grid diagnostics. The main one is a
#' check of the ends of the grid, i.e., its lower and upper limits. With
#' point masses, the model places no effect beyond them. For each end,
#' `grid$ends` gives the largest posterior probability that any site has in
#' the outer 2% of the grid's range (`mass`). If this probability exceeds
#' 0.05, the result for that site depends on where the grid ends, and the
#' `status` of the end depends on how the end was set (`set_by`).
#' Otherwise the status is `"ok"`.
#'
#' * The status is `"warning"` if the end was specified in `support` or
#'   follows from the estimates alone (`set_by` is `"support"` or
#'   `"estimates"`). [bayes_efron_fit()] then issues a warning of class
#'   `bef_warn_grid_ends` that names the site, and `"grid_ends"` is added
#'   to `sampler_diagnostics_warned`, although it concerns the grid and not
#'   the sampler. A wider range may change the result, and
#'   [compare_efron_grids()] shows whether it does.
#' * The status is `"note"` if the `"pilot"` method moved the end in
#'   (`set_by` is `"pilot"`). The pilot fit, which used a wider range,
#'   placed at most 0.05% of \eqn{g} beyond that end. A site with a large
#'   probability there has an estimate beyond the end of \eqn{g} and is
#'   shrunk to it, so the probability is reported without a warning. The
#'   same applies to a range around the true effects of a simulation
#'   (`set_by` is `"true effects"`).
#'
#' See the examples of [plot.bef_grid_comparison()] for a comparison in
#' which the narrowest range gives a probability above 0.05.
#'
#' No single number shows whether a range is wider than necessary. If
#' there are few sites, the data carry little information about the parts
#' of the range that contain no estimates, and \eqn{g} is determined by the
#' prior there. The `grid` component therefore gives the posterior mean and
#' 90% interval of the standard deviation of \eqn{g} (`sd_g`) next to the
#' 5%, 50% and 95% quantiles of its prior distribution on the same grid
#' (`sd_g_prior`). If the two are close, the data have hardly determined
#' the tails of \eqn{g}, and it is worth comparing fits with other ranges.
#' The proportion of the grid's range that contains the central 99.9% of
#' the estimate of \eqn{g} (`occupied`) is informative only if there are
#' many sites.
#'
#' @param fit A `bef_fit` object created by [bayes_efron_fit()].
#' @param ... For `diagnose()`, these dots are for future extensions and
#'   must be empty. The `print()` method passes them on to `format()`;
#'   `format()` and `summary()` ignore them.
#' @param x,object A `bef_diagnostic` object returned by `diagnose()`. For
#'   the bayesplot methods (`nuts_params()`, `log_posterior()`, `rhat()`
#'   and `neff_ratio()`), `object` is a `bef_fit` object.
#' @param use_cli If `TRUE` and the cli package is installed, the first
#'   line is formatted in bold; `FALSE` gives plain text. If `NULL` (the
#'   default), the value of `getOption("bayesEfron.use_cli", TRUE)` is
#'   used.
#'
#' @return An object of class `bef_diagnostic`, which is a list with the
#'   following components.
#'   \describe{
#'     \item{`rhat`}{R-hat for each variable in the draws, a named numeric
#'       vector.}
#'     \item{`ess_bulk`, `ess_tail`}{The bulk and tail effective sample
#'       sizes, in the same form.}
#'     \item{`divergences`, `max_treedepth`}{The number of divergent
#'       transitions and the number of iterations that reached the maximum
#'       tree depth, each summed over the chains.}
#'     \item{`ebfmi`}{The E-BFMI of each chain. Like the two counts, it is
#'       `NA` if CmdStan did not report it.}
#'     \item{`sampler_diagnostics_failed`}{Names of the diagnostics that
#'       fail a check.}
#'     \item{`sampler_diagnostics_warned`}{Names of the diagnostics that
#'       lead to a warning. `"grid_ends"` is among them if an end of the
#'       grid has the status `"warning"`.}
#'     \item{`diagnostic_skipped`}{Names of the diagnostics that could not
#'       be computed or were not checked. These three character vectors
#'       are empty if the run has at least 400 draws and there are no
#'       problems.}
#'     \item{`grid`}{A list with the grid diagnostics (see below), or
#'       `NULL` if they cannot be computed for the fit.}
#'     \item{`model_family`, `stan_file_sha256`, `runtime_seconds`,
#'       `effective_params_summary`}{Values copied from `fit$metadata` (see
#'       [bef_fit]).}
#'   }
#'   The list `grid` has the following components.
#'   \describe{
#'     \item{`ends`}{A data frame with one row per end of the grid (see
#'       below).}
#'     \item{`share`}{The outer fraction of the grid's range that is
#'       examined at each end (0.02).}
#'     \item{`threshold`}{The posterior probability above which an end is
#'       flagged (0.05).}
#'     \item{`sd_g`}{The posterior mean and the 5% and 95% quantiles of the
#'       standard deviation of \eqn{g}.}
#'     \item{`sd_g_prior`}{The 5%, 50% and 95% quantiles of the prior
#'       distribution of that standard deviation on this grid.}
#'     \item{`occupied`}{The proportion of the grid's range that contains
#'       the central 99.9% of the estimate of \eqn{g}.}
#'     \item{`tau_c_over_sd_g`}{The standard deviation of a normal
#'       component divided by the standard deviation of \eqn{g} (`NA` for
#'       point masses).}
#'   }
#'   The data frame `grid$ends` has the following columns.
#'   \describe{
#'     \item{`end`}{`"lower"` or `"upper"`.}
#'     \item{`value`}{The grid point at that end.}
#'     \item{`set_by`}{How the end was set: `"support"`, `"estimates"`,
#'       `"pilot"` or `"true effects"`.}
#'     \item{`mass`}{The largest posterior probability that any site has
#'       in the outer 2% of the range at that end.}
#'     \item{`site`, `label`}{The position in the data and the label of
#'       that site.}
#'     \item{`estimate_beyond`}{How far the estimate of that site lies
#'       beyond the end, in units of its standard error (negative if it
#'       lies inside the grid).}
#'     \item{`sites_over`}{The number of sites for which this probability
#'       exceeds 0.05.}
#'     \item{`limits_at_end`}{The number of sites whose 90% interval has
#'       a limit at the grid point at that end when `reading = "step"` is
#'       used, because that grid point alone has at least the tail
#'       probability of the interval.}
#'     \item{`status`}{`"ok"`, `"warning"` or `"note"`.}
#'   }
#'
#' @references
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Bürkner, P.-C.
#' (2021). Rank-normalization, folding, and localization: An improved
#' R-hat for assessing convergence of MCMC. *Bayesian Analysis*, 16(2).
#' \doi{10.1214/20-BA1221}
#'
#' @seealso [bayes_efron_fit()], [plot.bef_fit_re()], `vignette("diagnostics")`
#'
#' @examples
#' diag <- diagnose(krw_fit)
#' diag
#'
#' # The largest R-hat and the variable that has it.
#' summary(diag)$rhat
#'
#' # The three variables with the smallest bulk effective sample size.
#' head(sort(diag$ess_bulk), 3)
#'
#' @order 1
#' @export
diagnose <- function(fit, ...) {
  UseMethod("diagnose")
}

#' @order 2
#' @rdname diagnose
#' @export
diagnose.default <- function(fit, ...) {
  .bef_abort_invalid_fit(
    "`fit` must inherit from class \"bef_fit\".",
    arg = "fit"
  )
}

#' @order 3
#' @rdname diagnose
#' @export
diagnose.bef_fit <- function(fit, ...) {
  dots <- list(...)
  if (length(dots) > 0L) {
    .bef_abort_invalid_args(
      "`diagnose()` takes no arguments other than `fit`.",
      arg = "..."
    )
  }
  .bef_require_inherits(fit, "bef_fit", "fit", "bef_invalid_fit")
  .bef_fit_summary_version(fit)
  if (.bef_is_v01_fit(fit)) {
    .bef_abort_v01_fit("diagnose")
  }

  metadata <- fit$metadata
  diagnostics <- .bef_metadata_attr(metadata, "diagnostics")
  flags <- .bef_fit_diagnostic_flags(fit)
  grid <- .bef_fit_grid_diagnostics(fit)
  if (!is.null(grid) && any(grid$ends$status == "warning")) {
    flags$warned <- union(flags$warned, "grid_ends")
  }
  out <- new_bef_diagnostic(
    rhat = diagnostics$rhat,
    ess_bulk = diagnostics$ess_bulk,
    ess_tail = diagnostics$ess_tail,
    divergences = diagnostics$divergences,
    max_treedepth = diagnostics$max_treedepth,
    ebfmi = if (is.null(diagnostics$ebfmi)) NA_real_ else diagnostics$ebfmi,
    model_family = metadata$model_family,
    stan_file_sha256 = metadata$stan_file_sha256,
    effective_params_summary = metadata$effective_params_summary,
    runtime_seconds = metadata$runtime_seconds,
    diagnostic_skipped = flags$skipped,
    sampler_diagnostics_failed = flags$failed,
    sampler_diagnostics_warned = flags$warned,
    grid = grid
  )
  validate_bef_diagnostic(out)
}

#' @order 6
#' @describeIn diagnose Returns a data frame with one row per diagnostic
#'   and the columns `diagnostic`, `value` (the worst value), `threshold`,
#'   `rule` (whether `value` should be `"at most"` or `"at least"` the
#'   threshold) and `status` (`"ok"`, `"warning"`, `"failed"` or
#'   `"not checked"`). The last row, `grid_ends`, is for the end of the
#'   grid at which a site has the largest posterior probability in the
#'   outer 2% of the range. Its status is `"note"` if that probability
#'   exceeds the threshold at an end set by the `"pilot"` method or from
#'   the true effects.
#' @param row.names,optional Ignored. They are arguments of the
#'   `as.data.frame()` generic.
#' @export
as.data.frame.bef_diagnostic <- function(x, row.names = NULL, optional = FALSE, ...) {
  validate_bef_diagnostic(x)
  thresholds <- .bef_diagnostic_thresholds()
  worst <- c(
    rhat = .bef_diagnostic_worst(x$rhat, max),
    ess_bulk = .bef_diagnostic_worst(x$ess_bulk, min),
    ess_tail = .bef_diagnostic_worst(x$ess_tail, min),
    divergences = .bef_diagnostic_worst(x$divergences, sum),
    max_treedepth = .bef_diagnostic_worst(x$max_treedepth, sum),
    ebfmi = .bef_diagnostic_worst(x$ebfmi, min)
  )
  status <- rep("ok", length(worst))
  names(status) <- names(worst)
  status[names(worst) %in% x$sampler_diagnostics_warned] <- "warning"
  status[names(worst) %in% x$sampler_diagnostics_failed] <- "failed"
  unchecked <- sub("_check$", "", x$diagnostic_skipped)
  status[names(worst) %in% unchecked & status == "ok"] <- "not checked"
  if ("sampler_diagnostics" %in% x$diagnostic_skipped) {
    status[c("divergences", "max_treedepth", "ebfmi")] <- "not checked"
  }
  out <- data.frame(
    diagnostic = names(worst),
    value = unname(worst),
    threshold = c(thresholds$rhat_warn, thresholds$ess_bulk_min, thresholds$ess_tail_min,
                  thresholds$divergences_max, thresholds$max_treedepth_max,
                  thresholds$ebfmi_min),
    rule = c("at most", "at least", "at least", "at most", "at most", "at least"),
    status = unname(status),
    stringsAsFactors = FALSE
  )
  if (is.null(x$grid)) {
    return(out)
  }
  # The end of the grid with the larger share; of two, the one that is a
  # warning comes first.
  ends <- x$grid$ends
  order <- order(match(ends$status, c("warning", "note", "ok")), -ends$mass)
  worst_end <- ends[order[[1L]], ]
  rbind(out, data.frame(
    diagnostic = "grid_ends", value = worst_end$mass, threshold = x$grid$threshold,
    rule = "at most", status = worst_end$status, stringsAsFactors = FALSE
  ))
}

#' @order 7
#' @describeIn diagnose Methods for bayesplot generics, so that the
#'   bayesplot plots of sampler diagnostics can be made from a fit.
#'   `nuts_params()` returns the divergence indicator, the tree depth, the
#'   energy, the acceptance statistic, the step size and the number of
#'   leapfrog steps for every stored iteration, and `log_posterior()`
#'   returns the log density for every stored iteration. `rhat()` and
#'   `neff_ratio()` return, for each variable, R-hat and the bulk effective
#'   sample size divided by the number of iterations of the whole run. In a
#'   thinned fit (see `thin` in [bayes_efron_fit()]) the stored iterations
#'   are not consecutive, so a plot that compares each iteration with the
#'   previous one, such as `bayesplot::mcmc_nuts_energy()`, should not be
#'   made from it. `nuts_params()` throws an error if the fit does not
#'   contain the sampler quantities.
#' @param pars An optional character vector. For `nuts_params()`, the
#'   names of the sampler quantities to return; for `rhat()` and
#'   `neff_ratio()`, the names of the variables. If `NULL` (the default),
#'   all are returned.
#' @exportS3Method bayesplot::nuts_params
nuts_params.bef_fit <- function(object, pars = NULL, ...) {
  nuts <- .bef_metadata_attr(object$metadata, "nuts")
  if (is.null(nuts) || length(dim(nuts)) != 3L) {
    .bef_abort_invalid_fit(
      paste0(
        "This fit does not hold the sampler's quantities of each iteration, ",
        "which fits made by bayesEfron 0.4 or later keep. Fit the model again."
      )
    )
  }
  dims <- dim(nuts)
  parameters <- dimnames(nuts)[[3L]]
  if (!is.null(pars)) {
    parameters <- intersect(parameters, pars)
  }
  rows <- lapply(parameters, function(parameter) {
    data.frame(
      Chain = rep(seq_len(dims[[2L]]), each = dims[[1L]]),
      Iteration = rep(seq_len(dims[[1L]]), times = dims[[2L]]),
      Parameter = parameter,
      Value = as.numeric(nuts[, , parameter]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$Parameter <- factor(out$Parameter, levels = parameters)
  out
}

#' @rdname diagnose
#' @exportS3Method bayesplot::log_posterior
log_posterior.bef_fit <- function(object, ...) {
  dims <- dim(object$draws)
  if (!"lp__" %in% dimnames(object$draws)[[3L]]) {
    .bef_abort_invalid_fit("This fit does not hold `lp__`, the log density of each iteration.")
  }
  data.frame(
    Chain = rep(seq_len(dims[[2L]]), each = dims[[1L]]),
    Iteration = rep(seq_len(dims[[1L]]), times = dims[[2L]]),
    Value = as.numeric(object$draws[, , "lp__"])
  )
}

#' @rdname diagnose
#' @exportS3Method bayesplot::rhat
rhat.bef_fit <- function(object, pars = NULL, ...) {
  .bef_diagnostic_by_variable(object, "rhat", pars)
}

#' @rdname diagnose
#' @exportS3Method bayesplot::neff_ratio
neff_ratio.bef_fit <- function(object, pars = NULL, ...) {
  .bef_diagnostic_by_variable(object, "ess_bulk", pars) / .bef_run_draws(object)
}

# The number of sampling iterations of the whole run, over all chains. The
# stored diagnostics are those of the whole run, and a fit that was thinned
# keeps fewer draws than that. A fit made before version 0.4 does not record
# the length of its run, and was not thinned.
.bef_run_draws <- function(object) {
  settings <- .bef_metadata_attr(object$metadata, "sampler_settings")
  run <- suppressWarnings(
    as.numeric(settings$iter_sampling) * as.numeric(settings$chains)
  )
  if (length(run) == 1L && is.finite(run) && run > 0) {
    return(run)
  }
  prod(dim(object$draws)[seq_len(2L)])
}

.bef_diagnostic_by_variable <- function(object, field, pars) {
  values <- .bef_metadata_attr(object$metadata, "diagnostics")[[field]]
  if (!is.null(pars)) {
    values <- values[intersect(names(values), pars)]
  }
  values
}
