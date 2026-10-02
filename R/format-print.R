#' @order 13
#' @describeIn bef_fit Returns the lines that `print()` shows, as a character vector.
#' @export
format.bef_fit <- function(x, ..., use_cli = NULL) {
  x <- .bef_prepare_fit(x)
  if (.bef_format_use_cli(use_cli)) {
    return(format_bef_fit_cli(x))
  }
  format_bef_fit_base(x)
}

#' @order 14
#' @describeIn bef_fit Returns the lines of the printed summary.
#' @export
format.summary.bef_fit <- function(x, ..., use_cli = NULL) {
  .bef_check_saved_summary(x)
  if (.bef_format_use_cli(use_cli)) {
    return(format_summary_bef_fit_cli(x))
  }
  format_summary_bef_fit_base(x)
}

#' @order 8
#' @describeIn as_bef_data Returns the lines that `print()` shows, as a character vector.
#' @export
format.bef_data <- function(x, ..., use_cli = NULL) {
  if (.bef_format_use_cli(use_cli)) {
    return(format_bef_data_cli(x))
  }
  format_bef_data_base(x)
}

#' @order 6
#' @describeIn diagnose Returns the lines that `print()` shows, as a character vector.
#' @export
format.bef_diagnostic <- function(x, ..., use_cli = NULL) {
  if (.bef_format_use_cli(use_cli)) {
    return(format_bef_diagnostic_cli(x))
  }
  format_bef_diagnostic_base(x)
}

format_bef_fit_base <- function(x) {
  metadata <- x$metadata
  diagnostics <- .bef_metadata_attr(metadata, "diagnostics")
  skipped <- .bef_metadata_attr(metadata, "diagnostic_skipped")
  failed <- .bef_metadata_attr(metadata, "sampler_diagnostics_failed")
  warned <- .bef_metadata_attr(metadata, "sampler_diagnostics_warned")

  lines <- c(
    "<bayesEfron fit>",
    sprintf("Model family: %s", metadata$model_family),
    sprintf("Sites: %d", as.integer(metadata$data_list$K)),
    sprintf("Grid method: %s", metadata$grid_method),
    sprintf("Runtime: %s sec", .bef_format_number(metadata$runtime_seconds)),
    sprintf(
      "Diagnostics: Rhat %s; ESS bulk %s; ESS tail %s; divergences %s; treedepth hits %s",
      .bef_format_number(.bef_diagnostic_worst(diagnostics$rhat, max)),
      .bef_format_number(.bef_diagnostic_worst(diagnostics$ess_bulk, min)),
      .bef_format_number(.bef_diagnostic_worst(diagnostics$ess_tail, min)),
      .bef_format_integerish(.bef_diagnostic_worst(diagnostics$divergences, sum)),
      .bef_format_integerish(.bef_diagnostic_worst(diagnostics$max_treedepth, sum))
    )
  )

  if (length(skipped) > 0L) {
    lines <- c(lines, sprintf("Skipped diagnostics: %s", paste(skipped, collapse = ", ")))
  }
  if (length(failed) > 0L) {
    lines <- c(lines, sprintf("Failed checks: %s", paste(failed, collapse = ", ")))
  }
  if (length(warned) > 0L) {
    lines <- c(lines, sprintf("Warnings: %s", paste(warned, collapse = ", ")))
  }

  c(lines, "Use summary() for posterior summaries.")
}

format_summary_bef_fit_base <- function(x) {
  prior <- x$prior_summary
  diagnostics <- x$diagnostics
  lines <- c(
    "<summary.bef_fit>",
    "",
    "Prior g:",
    sprintf("  mean: %s", .bef_format_number(prior$mean)),
    sprintf("  var:  %s", .bef_format_number(prior$var)),
    sprintf("  sd:   %s", .bef_format_number(prior$sd)),
    "",
    "Diagnostics:",
    sprintf("  Rhat (max):        %s",
            .bef_format_number(.bef_diagnostic_worst(diagnostics$rhat, max))),
    sprintf("  ESS bulk (min):    %s",
            .bef_format_number(.bef_diagnostic_worst(diagnostics$ess_bulk, min))),
    sprintf("  ESS tail (min):    %s",
            .bef_format_number(.bef_diagnostic_worst(diagnostics$ess_tail, min))),
    sprintf("  Divergences:       %s",
            .bef_format_integerish(.bef_diagnostic_worst(diagnostics$divergences, sum))),
    sprintf("  Treedepth hits:    %s",
            .bef_format_integerish(.bef_diagnostic_worst(diagnostics$max_treedepth, sum))),
    sprintf(
      "  Effective params:  mean %s, sd %s",
      .bef_format_number(diagnostics$effective_params$mean),
      .bef_format_number(diagnostics$effective_params$sd)
    ),
    sprintf(
      "  Log marginal lik.: mean %s, sd %s",
      .bef_format_number(diagnostics$log_marginal_likelihood$mean),
      .bef_format_number(diagnostics$log_marginal_likelihood$sd)
    ),
    sprintf("  Runtime:           %s sec", .bef_format_number(diagnostics$runtime_seconds))
  )

  if (length(diagnostics$diagnostic_skipped) > 0L) {
    lines <- c(
      lines,
      sprintf(
        "  Skipped:           %s",
        paste(diagnostics$diagnostic_skipped, collapse = ", ")
      )
    )
  }
  if (length(diagnostics$sampler_diagnostics_failed) > 0L) {
    lines <- c(
      lines,
      sprintf(
        "  Failed checks:     %s",
        paste(diagnostics$sampler_diagnostics_failed, collapse = ", ")
      )
    )
  }
  if (length(diagnostics$sampler_diagnostics_warned) > 0L) {
    lines <- c(
      lines,
      sprintf(
        "  Warnings:          %s",
        paste(diagnostics$sampler_diagnostics_warned, collapse = ", ")
      )
    )
  }

  if (!is.null(x$theta_summary)) {
    lines <- c(lines, "", "Theta summary:", .bef_format_theta_summary(x$theta_summary))
  }

  lines
}

format_bef_data_base <- function(x) {
  x_summary <- summary(x)
  lines <- c(
    "<bef_data>",
    sprintf("Sites: %d", x_summary$K),
    sprintf("Source: %s", x_summary$source),
    sprintf(
      "theta_hat: min %s; median %s; max %s",
      .bef_format_number(x_summary$theta_hat$min),
      .bef_format_number(x_summary$theta_hat$median),
      .bef_format_number(x_summary$theta_hat$max)
    ),
    sprintf(
      "sigma: min %s; median %s; max %s",
      .bef_format_number(x_summary$sigma$min),
      .bef_format_number(x_summary$sigma$median),
      .bef_format_number(x_summary$sigma$max)
    )
  )
  if (!is.null(x_summary$names)) {
    shown <- utils::head(x_summary$names, 3L)
    suffix <- if (length(x_summary$names) > length(shown)) ", ..." else ""
    lines <- c(
      lines,
      sprintf("Names: %s%s", paste(shown, collapse = ", "), suffix)
    )
  }
  lines
}

format_bef_diagnostic_base <- function(x) {
  x_summary <- summary(x)
  lines <- c(
    "<bef_diagnostic>",
    sprintf("Model family: %s", x_summary$model_family),
    sprintf("Rhat max: %s%s", .bef_format_number(x_summary$rhat$value),
            .bef_format_at_variable(x_summary$rhat$variable)),
    sprintf("ESS bulk min: %s%s", .bef_format_number(x_summary$ess_bulk$value),
            .bef_format_at_variable(x_summary$ess_bulk$variable)),
    sprintf("ESS tail min: %s%s", .bef_format_number(x_summary$ess_tail$value),
            .bef_format_at_variable(x_summary$ess_tail$variable)),
    sprintf(
      "Divergences: %s",
      .bef_format_integerish(x_summary$divergences)
    ),
    sprintf(
      "Max treedepth hits: %s",
      .bef_format_integerish(x_summary$max_treedepth)
    ),
    sprintf("E-BFMI min: %s",
            .bef_format_number(.bef_diagnostic_worst(x_summary$ebfmi, min))),
    sprintf("Runtime: %s sec", .bef_format_number(x_summary$runtime_seconds))
  )

  if (length(x_summary$sampler_diagnostics_failed) > 0L) {
    lines <- c(lines, sprintf(
      "Failed checks: %s",
      paste(x_summary$sampler_diagnostics_failed, collapse = ", ")
    ))
  }
  if (length(x_summary$sampler_diagnostics_warned) > 0L) {
    lines <- c(lines, sprintf(
      "Warnings: %s",
      paste(x_summary$sampler_diagnostics_warned, collapse = ", ")
    ))
  }

  if (!is.null(x_summary$effective_params)) {
    lines <- c(
      lines,
      sprintf(
        "Effective params: mean %s; sd %s",
        .bef_format_number(x_summary$effective_params$mean),
        .bef_format_number(x_summary$effective_params$sd)
      )
    )
  }
  if (length(x_summary$diagnostic_skipped) > 0L) {
    lines <- c(
      lines,
      sprintf(
        "Skipped diagnostics: %s",
        paste(x_summary$diagnostic_skipped, collapse = ", ")
      )
    )
  }
  lines
}

format_bef_fit_cli <- function(x) {
  .bef_format_cli_lines(format_bef_fit_base(x), heading = 1L)
}

format_summary_bef_fit_cli <- function(x) {
  lines <- format_summary_bef_fit_base(x)
  .bef_format_cli_lines(
    lines,
    heading = 1L,
    sections = which(lines %in% c("Prior g:", "Diagnostics:", "Theta summary:"))
  )
}

format_bef_data_cli <- function(x) {
  .bef_format_cli_lines(format_bef_data_base(x), heading = 1L)
}

format_bef_diagnostic_cli <- function(x) {
  .bef_format_cli_lines(format_bef_diagnostic_base(x), heading = 1L)
}

.bef_format_use_cli <- function(use_cli) {
  if (!is.null(use_cli) &&
      (!is.logical(use_cli) || length(use_cli) != 1L || is.na(use_cli))) {
    .bef_abort_invalid_args(
      "`use_cli` must be NULL, TRUE, or FALSE.",
      arg = "use_cli"
    )
  }
  if (is.null(use_cli)) {
    use_cli <- getOption("bayesEfron.use_cli", TRUE)
  }
  isTRUE(use_cli) && requireNamespace("cli", quietly = TRUE)
}

.bef_format_cli_lines <- function(lines, heading = integer(), sections = integer()) {
  if (!requireNamespace("cli", quietly = TRUE)) {
    return(lines)
  }
  out <- lines
  if (length(heading) > 0L) {
    out[heading] <- cli::style_bold(out[heading])
  }
  if (length(sections) > 0L) {
    out[sections] <- cli::col_blue(out[sections])
  }
  out
}

.bef_format_theta_summary <- function(theta_summary, max_rows = 12L) {
  n <- nrow(theta_summary)
  rows <- seq_len(n)
  omitted <- 0L
  if (n > max_rows) {
    rows <- c(seq_len(max_rows / 2L), seq.int(n - max_rows / 2L + 1L, n))
    omitted <- n - length(rows)
  }

  out <- c("  site      mean        sd     lower     upper       map")
  body <- vapply(
    rows,
    function(i) {
      sprintf(
        "  %-4s %9s %9s %9s %9s %9s",
        as.character(theta_summary$site[[i]]),
        .bef_format_number(theta_summary$mean[[i]]),
        .bef_format_number(theta_summary$sd[[i]]),
        .bef_format_number(theta_summary$hpdi_lower[[i]]),
        .bef_format_number(theta_summary$hpdi_upper[[i]]),
        .bef_format_number(theta_summary$map[[i]])
      )
    },
    character(1L)
  )
  out <- c(out, body)
  if (omitted > 0L) {
    out <- c(out, sprintf("  ... %d sites omitted ...", omitted))
  }
  out
}

.bef_format_number <- function(x, digits = 4L) {
  if (length(x) != 1L || is.na(x)) {
    return("NA")
  }
  format(signif(as.numeric(x), digits = digits), trim = TRUE)
}

.bef_format_integerish <- function(x) {
  if (length(x) != 1L || is.na(x)) {
    return("NA")
  }
  as.character(as.integer(round(x)))
}

.bef_format_at_variable <- function(variable) {
  if (is.null(variable) || length(variable) != 1L || is.na(variable)) {
    return("")
  }
  sprintf(" (%s)", variable)
}
