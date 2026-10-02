#' Convergence diagnostics of a fitted model
#'
#' @description
#' `diagnose()` collects the diagnostics of the sampler that
#' [bayes_efron_fit()] computed when the model was fitted and names those
#' that fail a check. It computes nothing new from the draws.
#'
#' @details
#' A diagnostic is named in `sampler_diagnostics_failed` if an R-hat value
#' exceeds 1.05 or if there is a divergent transition. It is named in
#' `sampler_diagnostics_warned` if an R-hat value lies between 1.01 and
#' 1.05, a bulk or tail effective sample size is below 400, an iteration
#' reached the maximum tree depth, or the E-BFMI of a chain is below 0.2.
#' R-hat below 1.01 and effective sample sizes above 400 are the
#' recommendations of Vehtari et al. (2021); the other thresholds are
#' choices of this package.
#'
#' With fewer than 400 draws in all, R-hat and the effective sample sizes
#' are too unreliable to act on. Their values are still returned, but they
#' are not compared with the thresholds, and `diagnostic_skipped` holds
#' `"rhat_check"`, `"ess_bulk_check"` and `"ess_tail_check"` to say so.
#' Divergent transitions, the maximum tree depth and E-BFMI are checked
#' however short the run.
#'
#' R-hat and the effective sample sizes are undefined for a quantity that
#' takes the same value in every draw, such as the posterior mode of a site
#' whose mode never moves from one grid point, and are `NA` for it. If one
#' of the three diagnostics cannot be computed at all, or comes out
#' infinite for some quantity, all its values are `NA` and its name is in
#' `diagnostic_skipped`. `vignette("diagnostics")` discusses what to do
#' about each kind of warning.
#'
#' @param fit A `bef_fit` object from [bayes_efron_fit()].
#' @param ... For `diagnose()`, these dots are for future extensions and
#'   must be empty. The `print()` method passes them on to `format()`;
#'   `format()` and `summary()` do not use them.
#' @param x,object A `bef_diagnostic` object.
#' @inheritParams bef_fit
#'
#' @return A list of class `bef_diagnostic`. `rhat`, `ess_bulk` and
#'   `ess_tail` are numeric vectors with one value for each quantity in the
#'   draws. `divergences` and `max_treedepth` are the numbers of divergent
#'   transitions and of iterations that reached the maximum tree depth,
#'   totaled over the chains, and `ebfmi` has one value for each chain; the
#'   three are `NA` if CmdStan did not report them. `model_family`,
#'   `stan_file_sha256`, `runtime_seconds` and `effective_params_summary`
#'   repeat values from `fit$metadata`. `sampler_diagnostics_failed`,
#'   `sampler_diagnostics_warned` and `diagnostic_skipped` are character
#'   vectors that name the diagnostics that fail a check, those that
#'   deserve a warning, and those that could not be computed or were not
#'   checked (see Details). All three are empty for a fit of 400 or more
#'   draws without problems.
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
#' diag <- diagnose(raudenbush_fit)
#' diag
#'
#' # The worst value of R-hat and the quantity that has it.
#' summary(diag)$rhat
#'
#' # The three quantities with the smallest bulk effective sample size.
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
  validate_bef_diagnostic(out)
}
