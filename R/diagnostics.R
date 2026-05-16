#' Extract sampler diagnostics from a bayesEfron fit
#'
#' @description
#' Return a structured `bef_diagnostic` object that bundles the
#' sampler-health and model-quality summaries already computed by
#' [bayes_efron_fit()]. The function does not re-process draws or
#' run new computations; it surfaces information that the fit
#' pipeline stored as attributes of `fit$metadata`, behind the
#' supported access path.
#'
#' Use `diagnose()` together with `summary(fit)` to triage a fit:
#' `summary()` gives the user-facing posterior summaries; `diagnose()`
#' gives the structured sampler-health view that pairs with
#' `plot(fit, type = "diagnostic")`.
#'
#' @details
#' # The `bef_diagnostic` schema
#'
#' A `bef_diagnostic` object carries:
#'
#' | Field | Type | Meaning |
#' |:------|:-----|:--------|
#' | `rhat` | numeric vector | Per-parameter R-hat (Gelman–Rubin) statistic. Values close to 1 (typically `< 1.01`) indicate convergence. |
#' | `ess_bulk` | numeric vector | Per-parameter bulk effective sample size; large values support reliable posterior means. |
#' | `ess_tail` | numeric vector | Per-parameter tail ESS; large values support reliable tail-quantile / interval estimates. |
#' | `divergences` | numeric vector or NA | Number of divergent transitions per chain. |
#' | `max_treedepth` | numeric vector or NA | Number of max-treedepth saturations per chain. |
#' | `effective_params_summary` | named list | Posterior summary of effective parameters. |
#' | `model_family` | character | `"RE"` for v0.1. |
#' | `stan_file_sha256` | character | SHA-256 of the locked Stan source used for the fit. |
#' | `runtime_seconds` | numeric | Sampler wall-clock for the fit. |
#' | `diagnostic_skipped` | character vector | Diagnostics intentionally not computed (rare). |
#' | `sampler_diagnostics_failed` | character vector | Diagnostics requested but failed at extract time (rare). |
#'
#' The returned object is validated internally before return.
#'
#' # Reading the result
#'
#' Call `summary()` on a `bef_diagnostic` for the most-extreme value
#' per field (max R-hat, min ESS, total divergences, total
#' max-treedepth saturations); call `print()` for a one-screen
#' textual view. The `plot(fit, type = "diagnostic")` view consumes
#' the same attribute payload through this generic.
#'
#' @param fit A `bef_fit` object returned by [bayes_efron_fit()].
#' @param ... Reserved for future expansion; must be empty in v0.1.
#'
#' @return A validated `bef_diagnostic` object with the schema
#'   tabulated in \strong{Details}.
#'
#' @seealso
#'   * [bayes_efron_fit()] for the upstream pipeline that populates
#'     these diagnostics.
#'   * [summary.bef_diagnostic()][bayesEfron-methods] for the
#'     extreme-value summary.
#'   * [plot.bef_fit_re()] for the graphical companion view.
#'   * The methodological vignette M6 ("Verification and
#'     calibration") for what each diagnostic means in the context of
#'     the package's verification ledger.
#'
#' @examples
#' # Load the cached five-site smoke fit shipped with the package.
#' fit <- readRDS(system.file(
#'   "examples", "cached_fit_re_smoke.rds",
#'   package = "bayesEfron"
#' ))
#'
#' diag <- diagnose(fit)
#' print(diag)
#' summary(diag)
#'
#' @export
diagnose <- function(fit, ...) {
  UseMethod("diagnose")
}

#' @rdname diagnose
#' @export
diagnose.default <- function(fit, ...) {
  .bef_abort_invalid_fit(
    "`fit` must inherit from class \"bef_fit\".",
    arg = "fit",
    module = "diagnose"
  )
}

#' @rdname diagnose
#' @export
diagnose.bef_fit <- function(fit, ...) {
  dots <- list(...)
  if (length(dots) > 0L) {
    .bef_abort_invalid_args(
      "`diagnose()` does not accept additional arguments in bayesEfron v0.1.",
      arg = "...",
      module = "diagnose"
    )
  }
  .bef_require_inherits(fit, "bef_fit", "fit", "bef_invalid_fit")

  metadata <- fit$metadata
  diagnostics <- .bef_metadata_attr(metadata, "diagnostics")
  out <- new_bef_diagnostic(
    rhat = diagnostics$rhat,
    ess_bulk = diagnostics$ess_bulk,
    ess_tail = diagnostics$ess_tail,
    divergences = diagnostics$divergences,
    max_treedepth = diagnostics$max_treedepth,
    model_family = metadata$model_family,
    stan_file_sha256 = metadata$stan_file_sha256,
    effective_params_summary = metadata$effective_params_summary,
    runtime_seconds = metadata$runtime_seconds,
    diagnostic_skipped = .bef_metadata_attr(metadata, "diagnostic_skipped"),
    sampler_diagnostics_failed = .bef_metadata_attr(
      metadata,
      "sampler_diagnostics_failed"
    )
  )
  validate_bef_diagnostic(out)
}
