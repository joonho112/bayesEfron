# The bef_diagnostic class: constructor and validator.

new_bef_diagnostic <- function(rhat,
                               ess_bulk,
                               ess_tail,
                               divergences,
                               max_treedepth,
                               ...,
                               model_family,
                               stan_file_sha256,
                               effective_params_summary = NULL,
                               runtime_seconds = NULL,
                               diagnostic_skipped = character(),
                               sampler_diagnostics_failed = character(),
                               sampler_diagnostics_warned = character(),
                               ebfmi = NA_real_,
                               grid = NULL) {
  structure(
    list(
      rhat = rhat,
      ess_bulk = ess_bulk,
      ess_tail = ess_tail,
      divergences = divergences,
      max_treedepth = max_treedepth,
      ebfmi = ebfmi,
      effective_params_summary = effective_params_summary,
      model_family = model_family,
      stan_file_sha256 = stan_file_sha256,
      runtime_seconds = runtime_seconds,
      diagnostic_skipped = diagnostic_skipped,
      sampler_diagnostics_failed = sampler_diagnostics_failed,
      sampler_diagnostics_warned = sampler_diagnostics_warned,
      # What the fit says about its grid; `NULL` where it could not be
      # computed.
      grid = grid
    ),
    class = "bef_diagnostic"
  )
}

validate_bef_diagnostic <- function(x) {
  class <- "bef_invalid_fit"
  .bef_require_inherits(x, "bef_diagnostic", "x", class)
  .bef_require_fields(
    x, .bef_diagnostic_fields(), "`bef_diagnostic`", class
  )

  .bef_validate_diagnostic_numeric(
    x$rhat, "rhat", class, lower = 0, open_lower = TRUE
  )
  .bef_validate_diagnostic_numeric(
    x$ess_bulk, "ess_bulk", class, lower = 0
  )
  .bef_validate_diagnostic_numeric(
    x$ess_tail, "ess_tail", class, lower = 0
  )
  .bef_validate_diagnostic_integerish(
    x$divergences, "divergences", class, lower = 0
  )
  .bef_validate_diagnostic_integerish(
    x$max_treedepth, "max_treedepth", class, lower = 0
  )
  .bef_validate_diagnostic_numeric(
    x$ebfmi, "ebfmi", class, lower = 0
  )

  if (!identical(x$model_family, "RE")) {
    .bef_abort_validate(
      "`bef_diagnostic$model_family` must be \"RE\"; other families are not implemented.",
      class
    )
  }
  .bef_validate_sha256(x$stan_file_sha256, "`bef_diagnostic$stan_file_sha256`", class)

  if (!is.null(x$runtime_seconds) &&
      (!.bef_is_number(x$runtime_seconds) || x$runtime_seconds < 0)) {
    .bef_abort_validate(
      "`bef_diagnostic$runtime_seconds` must be NULL or a non-negative finite numeric scalar.",
      class
    )
  }
  if (!is.character(x$diagnostic_skipped) || anyNA(x$diagnostic_skipped)) {
    .bef_abort_validate(
      "`bef_diagnostic$diagnostic_skipped` must be a character vector without missing values.",
      class
    )
  }
  if (!is.character(x$sampler_diagnostics_warned) ||
      anyNA(x$sampler_diagnostics_warned)) {
    .bef_abort_validate(
      "`bef_diagnostic$sampler_diagnostics_warned` must be a character vector without missing values.",
      class
    )
  }
  if (!is.character(x$sampler_diagnostics_failed) ||
      anyNA(x$sampler_diagnostics_failed)) {
    .bef_abort_validate(
      "`bef_diagnostic$sampler_diagnostics_failed` must be a character vector without missing values.",
      class
    )
  }
  if (!is.null(x$effective_params_summary)) {
    .bef_validate_summary_list(
      x$effective_params_summary, "effective_params_summary", class
    )
  }
  if (!is.null(x$grid) && !(is.list(x$grid) && is.data.frame(x$grid$ends) &&
                              nrow(x$grid$ends) == 2L)) {
    .bef_abort_validate(
      "`bef_diagnostic$grid` must be NULL or a list whose element `ends` has one row for each end of the grid.",
      class
    )
  }

  x
}

.bef_diagnostic_fields <- function() {
  c(
    "rhat",
    "ess_bulk",
    "ess_tail",
    "divergences",
    "max_treedepth",
    "ebfmi",
    "effective_params_summary",
    "model_family",
    "stan_file_sha256",
    "runtime_seconds",
    "diagnostic_skipped",
    "sampler_diagnostics_failed",
    "sampler_diagnostics_warned",
    "grid"
  )
}

# The lines that print() shows about the grid of a fit.
.bef_format_grid_diagnostics <- function(grid, K = NA_integer_) {
  if (is.null(grid)) {
    return(character())
  }
  ends <- grid$ends
  percent <- function(x) format(round(100 * x, 1), trim = TRUE)
  end_line <- function(i) {
    note <- switch(
      ends$status[[i]],
      warning = sprintf("; the result for that site depends on where the grid ends"),
      note = sprintf("; that end was set %s", if (identical(ends$set_by[[i]], "pilot")) {
        "by the pilot fit on a wider range"
      } else {
        "from the true effects"
      }),
      ""
    )
    sprintf("  %s end %s: at most %s%% of a site's posterior in the outer %s%% (site %s)%s",
            ends$end[[i]], .bef_format_number(ends$value[[i]]), percent(ends$mass[[i]]),
            format(100 * grid$share, trim = TRUE), .bef_quote_label(ends$label[[i]]), note)
  }
  lines <- c("Grid:", end_line(1L), end_line(2L))
  lines <- c(lines, sprintf(
    "  sd of g: posterior %s (%s to %s)%s",
    .bef_format_number(grid$sd_g[["mean"]]), .bef_format_number(grid$sd_g[["q5"]]),
    .bef_format_number(grid$sd_g[["q95"]]),
    if (is.finite(grid$sd_g_prior[["q50"]])) {
      sprintf("; prior %s (%s to %s)", .bef_format_number(grid$sd_g_prior[["q50"]]),
              .bef_format_number(grid$sd_g_prior[["q5"]]),
              .bef_format_number(grid$sd_g_prior[["q95"]]))
    } else {
      ""
    }
  ))
  lines <- c(lines, sprintf("  the central 99.9%% of g spans %s%% of the grid",
                            percent(grid$occupied)))
  if (is.finite(grid$tau_c_over_sd_g)) {
    lines <- c(lines, sprintf("  a normal component is %s times as wide as the sd of g",
                              .bef_format_number(grid$tau_c_over_sd_g)))
  }
  lines
}

.bef_sampler_diagnostic_warning_fields <- function() {
  c("rhat", "ess_bulk", "ess_tail", "max_treedepth", "ebfmi", "grid_ends")
}

.bef_diagnostic_skipped_fields <- function() {
  c(
    "rhat", "ess_bulk", "ess_tail", "sampler_diagnostics",
    "rhat_check", "ess_bulk_check", "ess_tail_check"
  )
}

.bef_sampler_diagnostic_failure_fields <- function() {
  c("rhat", "ess_bulk", "ess_tail", "divergences")
}

.bef_validate_diagnostic_skip_consistency <- function(diagnostics,
                                                      diagnostic_skipped,
                                                      class) {
  # Rhat and the effective sample sizes are undefined for a variable that is
  # constant across the draws (a site whose posterior mode never moves, for
  # example), so some NA values are normal. A diagnostic counts as skipped
  # only when all of its values are NA. Names ending in _check instead
  # identify thresholds that were not applied to a short run; their
  # numerical values can still be shown.
  scalar_skipped <- intersect(diagnostic_skipped, c("rhat", "ess_bulk", "ess_tail"))
  for (field in scalar_skipped) {
    if (!all(is.na(diagnostics[[field]]))) {
      .bef_abort_validate(
        sprintf("Skipped diagnostic `%s` must have an `NA` diagnostic value.", field),
        class
      )
    }
  }

  if ("sampler_diagnostics" %in% diagnostic_skipped &&
      (!any(is.na(diagnostics$divergences)) ||
       !any(is.na(diagnostics$max_treedepth)))) {
    .bef_abort_validate(
      "Skipped sampler diagnostics must set divergences and max_treedepth to `NA`.",
      class
    )
  }

  unexpected_na <- c(
    setdiff(c("rhat", "ess_bulk", "ess_tail"), diagnostic_skipped)[
      vapply(
        setdiff(c("rhat", "ess_bulk", "ess_tail"), diagnostic_skipped),
        function(field) all(is.na(diagnostics[[field]])),
        logical(1)
      )
    ],
    if (!"sampler_diagnostics" %in% diagnostic_skipped &&
        (any(is.na(diagnostics$divergences)) ||
         any(is.na(diagnostics$max_treedepth)))) {
      "sampler_diagnostics"
    }
  )
  if (length(unexpected_na) > 0L) {
    .bef_abort_validate(
      "Diagnostic `NA` values must be recorded in attr(metadata, \"diagnostic_skipped\").",
      class,
      diagnostics = unexpected_na
    )
  }
}

.bef_validate_diagnostic_numeric <- function(x,
                                             field,
                                             class,
                                             lower = -Inf,
                                             open_lower = FALSE) {
  if (!is.numeric(x) || length(x) < 1L ||
      any(is.infinite(x)) ||
      any(!is.na(x) & if (open_lower) x <= lower else x < lower)) {
    .bef_abort_validate(
      sprintf("`bef_diagnostic$%s` must be numeric and respect its lower bound.", field),
      class
    )
  }
}

.bef_validate_diagnostic_integerish <- function(x, field, class, lower = 0) {
  if (!is.numeric(x) || length(x) < 1L || any(is.infinite(x)) ||
      any(!is.na(x) & x != as.integer(x)) ||
      any(!is.na(x) & x < lower)) {
    .bef_abort_validate(
      sprintf("`bef_diagnostic$%s` must be an integerish non-negative value.", field),
      class
    )
  }
}
