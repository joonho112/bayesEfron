postprocess_stan_draws <- function(cmdstan_fit, stan_data, model_family = "RE") {
  if (!identical(model_family, "RE")) {
    .bef_abort_extraction_failed(
      "`model_family` must be \"RE\" for bayesEfron v0.1 postprocessing.",
      stage = 8L,
      model_family = model_family
    )
  }

  K <- as.integer(stan_data$K)
  draws <- .bef_extract_draws_array(cmdstan_fit)
  draw_matrix <- .bef_draws_matrix(draws)
  posterior <- .bef_extract_re_generated_quantities(draw_matrix, K = K)
  diagnostics <- .bef_postprocess_diagnostics(
    cmdstan_fit = cmdstan_fit,
    draws = draws
  )

  list(
    draws = draws,
    posterior = posterior,
    diagnostics = diagnostics$values,
    diagnostic_skipped = diagnostics$diagnostic_skipped,
    sampler_diagnostics_failed = diagnostics$sampler_diagnostics_failed,
    mean_g_summary = .bef_summary_vector(posterior$mean_g),
    var_g_summary = .bef_summary_vector(posterior$var_g),
    sd_g_summary = .bef_summary_vector(posterior$sd_g),
    effective_params_summary = .bef_summary_vector(posterior$effective_params),
    log_marginal_likelihood_summary = .bef_summary_vector(
      posterior$log_marginal_likelihood
    ),
    theta_summary = .bef_theta_summary(posterior),
    theta_rep_draws = posterior$theta_rep
  )
}

.bef_postprocess_diagnostics <- function(cmdstan_fit, draws) {
  values <- list(
    rhat = .bef_draw_diagnostic(draws, "rhat", posterior::rhat),
    ess_bulk = .bef_draw_diagnostic(draws, "ess_bulk", posterior::ess_bulk),
    ess_tail = .bef_draw_diagnostic(draws, "ess_tail", posterior::ess_tail)
  )

  sampler <- .bef_sampler_diagnostic_counts(cmdstan_fit)
  values$divergences <- sampler$divergences
  values$max_treedepth <- sampler$max_treedepth

  skipped <- c(
    attr(values$rhat, "bef_skipped_diagnostic", exact = TRUE),
    attr(values$ess_bulk, "bef_skipped_diagnostic", exact = TRUE),
    attr(values$ess_tail, "bef_skipped_diagnostic", exact = TRUE),
    sampler$diagnostic_skipped
  )
  values <- lapply(values, .bef_drop_diagnostic_attrs)

  failed <- .bef_sampler_diagnostic_failures(values, draws)
  if (length(failed) > 0L) {
    .bef_warn_sampler_diagnostics_failed(
      "Sampler diagnostics exceeded bayesEfron warning thresholds.",
      diagnostics = failed
    )
  }

  list(
    values = values,
    diagnostic_skipped = unique(skipped),
    sampler_diagnostics_failed = failed
  )
}

.bef_draw_diagnostic <- function(draws, diagnostic, fun) {
  value <- tryCatch(
    fun(draws),
    error = function(err) {
      .bef_warn_diagnostic_skipped(
        sprintf("Could not compute sampler diagnostic `%s`.", diagnostic),
        diagnostic = diagnostic,
        parent = err
      )
      out <- NA_real_
      attr(out, "bef_skipped_diagnostic") <- diagnostic
      out
    }
  )

  skipped <- attr(value, "bef_skipped_diagnostic", exact = TRUE)
  value <- suppressWarnings(as.numeric(value))
  if (!is.null(skipped)) {
    attr(value, "bef_skipped_diagnostic") <- skipped
    return(value)
  }
  if (length(value) != 1L || is.infinite(value)) {
    .bef_warn_diagnostic_skipped(
      sprintf("Sampler diagnostic `%s` did not return a finite scalar.", diagnostic),
      diagnostic = diagnostic
    )
    out <- NA_real_
    attr(out, "bef_skipped_diagnostic") <- diagnostic
    return(out)
  }
  if (is.na(value)) {
    .bef_warn_diagnostic_skipped(
      sprintf("Sampler diagnostic `%s` returned `NA`.", diagnostic),
      diagnostic = diagnostic
    )
    attr(value, "bef_skipped_diagnostic") <- diagnostic
  }
  value
}

.bef_sampler_diagnostic_counts <- function(cmdstan_fit) {
  diagnostic_fun <- tryCatch(
    cmdstan_fit$diagnostic_summary,
    error = function(err) NULL
  )
  if (!is.function(diagnostic_fun)) {
    .bef_warn_diagnostic_skipped(
      "CmdStan fit does not expose `diagnostic_summary()`; divergence and treedepth diagnostics are unavailable.",
      diagnostic = "sampler_diagnostics"
    )
    return(list(
      divergences = NA_real_,
      max_treedepth = NA_real_,
      diagnostic_skipped = "sampler_diagnostics"
    ))
  }

  summary <- tryCatch(
    diagnostic_fun(
      diagnostics = c("divergences", "treedepth"),
      quiet = TRUE
    ),
    error = function(err) {
      .bef_warn_diagnostic_skipped(
        "Could not extract CmdStan sampler diagnostics.",
        diagnostic = "sampler_diagnostics",
        parent = err
      )
      NULL
    }
  )
  if (is.null(summary)) {
    return(list(
      divergences = NA_real_,
      max_treedepth = NA_real_,
      diagnostic_skipped = "sampler_diagnostics"
    ))
  }

  divergences <- .bef_integerish_diagnostic(summary$num_divergent)
  max_treedepth <- .bef_integerish_diagnostic(summary$num_max_treedepth)
  if (is.na(divergences) || is.na(max_treedepth)) {
    .bef_warn_diagnostic_skipped(
      "CmdStan sampler diagnostics did not return finite non-negative counts.",
      diagnostic = "sampler_diagnostics"
    )
    return(list(
      divergences = NA_real_,
      max_treedepth = NA_real_,
      diagnostic_skipped = "sampler_diagnostics"
    ))
  }

  list(
    divergences = divergences,
    max_treedepth = max_treedepth,
    diagnostic_skipped = character()
  )
}

.bef_integerish_diagnostic <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0L || any(is.infinite(x)) || any(!is.na(x) & x < 0)) {
    return(NA_real_)
  }
  if (any(is.na(x))) {
    return(NA_real_)
  }
  sum(x)
}

.bef_drop_diagnostic_attrs <- function(x) {
  attributes(x) <- NULL
  x
}

.bef_sampler_diagnostic_failures <- function(values, draws) {
  n_draws <- prod(dim(draws)[seq_len(2L)])
  if (!is.finite(n_draws) || n_draws < 400L) {
    return(character())
  }

  failed <- character()
  if (is.finite(values$rhat) && values$rhat > 1.05) {
    failed <- c(failed, "rhat")
  }
  if (is.finite(values$ess_bulk) && values$ess_bulk < 400) {
    failed <- c(failed, "ess_bulk")
  }
  if (is.finite(values$ess_tail) && values$ess_tail < 400) {
    failed <- c(failed, "ess_tail")
  }
  if (is.finite(values$divergences) && values$divergences / n_draws > 0.01) {
    failed <- c(failed, "divergences")
  }
  failed
}

.bef_extract_draws_array <- function(cmdstan_fit) {
  draws_fun <- tryCatch(cmdstan_fit$draws, error = function(err) NULL)
  if (!is.function(draws_fun)) {
    .bef_abort_extraction_failed(
      "`cmdstan_fit` does not expose a callable `draws()` method.",
      stage = 8L
    )
  }

  draws <- tryCatch(
    draws_fun(format = "draws_array"),
    error = function(err) {
      .bef_abort_extraction_failed(
        "Failed to extract CmdStan draws as a draws_array.",
        stage = 8L,
        parent = err
      )
    }
  )

  if (!is.array(draws) || !is.numeric(draws) || length(dim(draws)) != 3L ||
      any(!is.finite(draws))) {
    .bef_abort_extraction_failed(
      "`cmdstan_fit$draws(format = \"draws_array\")` must return a finite numeric 3D array.",
      stage = 8L
    )
  }

  draws
}

.bef_draws_matrix <- function(draws) {
  variables <- dimnames(draws)[[3L]]
  if (!is.character(variables) || anyNA(variables) || any(!nzchar(variables))) {
    .bef_abort_extraction_failed(
      "Draws array must carry variable names in its third dimension.",
      stage = 8L
    )
  }

  posterior::as_draws_matrix(draws)
}

.bef_extract_re_generated_quantities <- function(draw_matrix, K) {
  list(
    mean_g = .bef_draws_scalar(draw_matrix, "mean_g"),
    var_g = .bef_draws_scalar(draw_matrix, "var_g"),
    sd_g = .bef_draws_scalar(draw_matrix, "sd_g"),
    theta_map = .bef_draws_vector(draw_matrix, "theta_map", K),
    theta_mean = .bef_draws_vector(draw_matrix, "theta_mean", K),
    theta_sd = .bef_draws_vector(draw_matrix, "theta_sd", K),
    theta_rep = .bef_draws_vector(draw_matrix, "theta_rep", K),
    effective_params = .bef_draws_scalar(draw_matrix, "effective_params"),
    log_marginal_likelihood = .bef_draws_scalar(
      draw_matrix,
      "log_marginal_likelihood"
    )
  )
}

.bef_draws_scalar <- function(draw_matrix, field) {
  if (!field %in% colnames(draw_matrix)) {
    .bef_abort_extraction_failed(
      sprintf("Draws array is missing generated quantity `%s`.", field),
      stage = 8L,
      field = field
    )
  }
  as.numeric(draw_matrix[, field])
}

.bef_draws_vector <- function(draw_matrix, field, K) {
  columns <- sprintf("%s[%d]", field, seq_len(K))
  missing <- setdiff(columns, colnames(draw_matrix))
  if (length(missing) > 0L) {
    .bef_abort_extraction_failed(
      sprintf("Draws array is missing generated quantity columns for `%s`.", field),
      stage = 8L,
      field = field,
      missing_fields = missing
    )
  }
  out <- draw_matrix[, columns, drop = FALSE]
  .bef_plain_draw_matrix(out)
}

.bef_plain_draw_matrix <- function(x) {
  matrix(
    as.numeric(x),
    nrow = nrow(x),
    ncol = ncol(x),
    dimnames = dimnames(x)
  )
}

.bef_summary_vector <- function(x) {
  x <- as.numeric(x)
  q <- stats::quantile(x, probs = c(0.05, 0.5, 0.95), names = FALSE)
  list(
    mean = mean(x),
    sd = .bef_sd0(x),
    q5 = q[[1L]],
    q50 = q[[2L]],
    q95 = q[[3L]]
  )
}

.bef_sd0 <- function(x) {
  if (length(x) <= 1L) {
    return(0)
  }
  stats::sd(x)
}

.bef_theta_summary <- function(posterior) {
  K <- ncol(posterior$theta_mean)
  hpdi <- vapply(
    seq_len(K),
    function(site) {
      posterior::quantile2(
        posterior$theta_rep[, site],
        probs = c(0.05, 0.95),
        names = FALSE
      )
    },
    numeric(2L)
  )

  data.frame(
    site = seq_len(K),
    mean = colMeans(posterior$theta_mean),
    sd = colMeans(posterior$theta_sd),
    hpdi_lower = hpdi[1L, ],
    hpdi_upper = hpdi[2L, ],
    map = colMeans(posterior$theta_map),
    row.names = seq_len(K),
    check.names = FALSE
  )
}

.bef_matrix_quantiles <- function(x, probs) {
  vapply(
    seq_len(ncol(x)),
    function(col) {
      posterior::quantile2(x[, col], probs = probs, names = FALSE)
    },
    numeric(length(probs))
  )
}
