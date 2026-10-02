postprocess_stan_draws <- function(cmdstan_fit, stan_data, model_family = "RE") {
  if (!identical(model_family, "RE")) {
    .bef_abort_extraction_failed(
      "`model_family` must be \"RE\"; other families are not implemented.",
      model_family = model_family
    )
  }

  K <- as.integer(stan_data$K)
  draws <- .bef_extract_draws_array(cmdstan_fit)
  draw_matrix <- .bef_draws_matrix(draws)
  # All generated quantities are needed for the summaries. Only the scalar
  # series are kept in `posterior`; the site-level matrices stay in `draws`.
  posterior_full <- .bef_extract_re_generated_quantities(draw_matrix, K = K)
  posterior <- posterior_full[.bef_generated_quantity_fields()]
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
    sampler_diagnostics_warned = diagnostics$sampler_diagnostics_warned,
    mean_g_summary = .bef_summary_vector(posterior_full$mean_g),
    var_g_summary = .bef_summary_vector(posterior_full$var_g),
    sd_g_summary = .bef_summary_vector(posterior_full$sd_g),
    effective_params_summary = .bef_summary_vector(posterior_full$effective_params),
    log_marginal_likelihood_summary = .bef_summary_vector(
      posterior_full$log_marginal_likelihood
    ),
    theta_summary = .bef_theta_summary(
      posterior_full,
      log_g = .bef_draws_vector(draw_matrix, "log_g", as.integer(stan_data$L)),
      stan_data = stan_data
    )
  )
}

# The worst value of a diagnostic over parameters or chains: the largest Rhat,
# the smallest effective sample size or E-BFMI, the total count of divergent
# transitions or of iterations that reached the maximum tree depth. `NA` when
# no finite value is available, which is how a skipped diagnostic appears.
.bef_diagnostic_worst <- function(x, fun) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else fun(x)
}

# The thresholds used by diagnose(), kept in one place. Rhat above `rhat_fail`
# fails the check and Rhat between `rhat_warn` and `rhat_fail` is a warning;
# the other values are the smallest acceptable effective sample size and
# E-BFMI and the largest acceptable counts.
.bef_diagnostic_thresholds <- function() {
  list(
    rhat_warn = 1.01,
    rhat_fail = 1.05,
    ess_bulk_min = 400,
    ess_tail_min = 400,
    divergences_max = 0L,
    max_treedepth_max = 0L,
    ebfmi_min = 0.2,
    min_draws_for_check = 400L
  )
}

# Rhat and the bulk and tail effective sample sizes, one value per monitored
# variable. They must be computed variable by variable: a single Rhat over the
# whole draws array is not a convergence diagnostic.
.bef_per_parameter_diagnostics <- function(draws) {
  summary <- tryCatch(
    posterior::summarize_draws(draws, "rhat", "ess_bulk", "ess_tail"),
    error = function(err) {
      .bef_warn_diagnostic_skipped(
        "Could not compute per-parameter convergence diagnostics.",
        diagnostic = "per_parameter",
        parent = err
      )
      NULL
    }
  )
  if (is.null(summary)) {
    na <- stats::setNames(NA_real_, "all")
    return(list(
      rhat = na, ess_bulk = na, ess_tail = na,
      diagnostic_skipped = c("rhat", "ess_bulk", "ess_tail")
    ))
  }

  # With fewer draws than `min_draws_for_check`, Rhat and the effective sample
  # sizes are not meaningful, so a missing value is not reported as a problem.
  n_draws <- prod(dim(draws)[seq_len(2L)])
  quiet <- !is.finite(n_draws) ||
    n_draws < .bef_diagnostic_thresholds()$min_draws_for_check

  summary <- as.data.frame(summary)
  nm <- summary$variable
  pull <- function(field) {
    value <- suppressWarnings(as.numeric(summary[[field]]))
    stats::setNames(value, nm)
  }
  rhat <- pull("rhat")
  ess_bulk <- pull("ess_bulk")
  ess_tail <- pull("ess_tail")

  skipped <- character()
  for (field in c("rhat", "ess_bulk", "ess_tail")) {
    value <- get(field)
    if (all(is.na(value)) || any(is.infinite(value))) {
      # A nonfinite result cannot satisfy the stored diagnostic contract.
      # Mark the entire field unavailable, just as for a computation error;
      # a skipped field must contain only NA values.
      if (any(is.infinite(value))) {
        value[] <- NA_real_
        assign(field, value)
      }
      if (!quiet) {
        .bef_warn_diagnostic_skipped(
          sprintf("Convergence diagnostic `%s` is unavailable.", field),
          diagnostic = field
        )
      }
      skipped <- c(skipped, field)
    }
  }

  list(
    rhat = rhat, ess_bulk = ess_bulk, ess_tail = ess_tail,
    diagnostic_skipped = skipped
  )
}

.bef_postprocess_diagnostics <- function(cmdstan_fit, draws) {
  per_parameter <- .bef_per_parameter_diagnostics(draws)

  values <- list(
    rhat = per_parameter$rhat,
    ess_bulk = per_parameter$ess_bulk,
    ess_tail = per_parameter$ess_tail
  )

  sampler <- .bef_sampler_diagnostic_counts(cmdstan_fit)
  values$divergences <- sampler$divergences
  values$max_treedepth <- sampler$max_treedepth
  values$ebfmi <- sampler$ebfmi

  checks <- .bef_diagnostic_check(values, draws)
  skipped <- unique(c(
    per_parameter$diagnostic_skipped, sampler$diagnostic_skipped,
    checks$skipped
  ))
  if (length(checks$failed) > 0L) {
    .bef_warn_sampler_diagnostics_failed(
      "Some sampler diagnostics fail the convergence checks; see diagnose().",
      diagnostics = checks$failed
    )
  }

  list(
    values = values,
    diagnostic_skipped = skipped,
    sampler_diagnostics_failed = checks$failed,
    sampler_diagnostics_warned = checks$warned
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
      ebfmi = NA_real_,
      diagnostic_skipped = "sampler_diagnostics"
    ))
  }

  summary <- tryCatch(
    diagnostic_fun(
      diagnostics = c("divergences", "treedepth", "ebfmi"),
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
      ebfmi = NA_real_,
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
      ebfmi = NA_real_,
      diagnostic_skipped = "sampler_diagnostics"
    ))
  }

  ebfmi <- suppressWarnings(as.numeric(summary$ebfmi))
  if (length(ebfmi) == 0L || all(is.na(ebfmi))) {
    ebfmi <- NA_real_
  }

  list(
    divergences = divergences,
    max_treedepth = max_treedepth,
    ebfmi = ebfmi,
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

# Compare the diagnostics with the thresholds. Returns the names of the
# diagnostics that fail, deserve a warning or have unchecked thresholds.
.bef_diagnostic_check <- function(values, draws) {
  thresholds <- .bef_diagnostic_thresholds()
  n_draws <- prod(dim(draws)[seq_len(2L)])
  check_series <- is.finite(n_draws) &&
    n_draws >= thresholds$min_draws_for_check

  failed <- character()
  warned <- character()
  skipped <- if (check_series) character() else
    c("rhat_check", "ess_bulk_check", "ess_tail_check")

  finite_of <- function(x) x[is.finite(x)]

  if (check_series) {
    rhat <- finite_of(values$rhat)
    if (length(rhat) > 0L) {
      if (any(rhat > thresholds$rhat_fail)) {
        failed <- c(failed, "rhat")
      } else if (any(rhat > thresholds$rhat_warn)) {
        warned <- c(warned, "rhat")
      }
    }
    ess_bulk <- finite_of(values$ess_bulk)
    if (length(ess_bulk) > 0L && any(ess_bulk < thresholds$ess_bulk_min)) {
      warned <- c(warned, "ess_bulk")
    }
    ess_tail <- finite_of(values$ess_tail)
    if (length(ess_tail) > 0L && any(ess_tail < thresholds$ess_tail_min)) {
      warned <- c(warned, "ess_tail")
    }
  }

  divergences <- finite_of(values$divergences)
  if (length(divergences) > 0L &&
      sum(divergences) > thresholds$divergences_max) {
    failed <- c(failed, "divergences")
  }

  max_treedepth <- finite_of(values$max_treedepth)
  if (length(max_treedepth) > 0L &&
      sum(max_treedepth) > thresholds$max_treedepth_max) {
    warned <- c(warned, "max_treedepth")
  }

  ebfmi <- finite_of(values$ebfmi)
  if (length(ebfmi) > 0L && any(ebfmi < thresholds$ebfmi_min)) {
    warned <- c(warned, "ebfmi")
  }

  list(failed = failed, warned = warned, skipped = skipped)
}

.bef_extract_draws_array <- function(cmdstan_fit) {
  draws_fun <- tryCatch(cmdstan_fit$draws, error = function(err) NULL)
  if (!is.function(draws_fun)) {
    .bef_abort_extraction_failed(
      "`cmdstan_fit` does not expose a callable `draws()` method."
    )
  }

  draws <- tryCatch(
    draws_fun(format = "draws_array"),
    error = function(err) {
      .bef_abort_extraction_failed(
        "Failed to extract CmdStan draws as a draws_array.",
        parent = err
      )
    }
  )

  if (!is.array(draws) || !is.numeric(draws) || length(dim(draws)) != 3L ||
      any(!is.finite(draws))) {
    .bef_abort_extraction_failed(
      "`cmdstan_fit$draws(format = \"draws_array\")` must return a finite numeric 3D array."
    )
  }

  draws
}

.bef_draws_matrix <- function(draws) {
  variables <- dimnames(draws)[[3L]]
  if (!is.character(variables) || anyNA(variables) || any(!nzchar(variables))) {
    .bef_abort_extraction_failed(
      "Draws array must carry variable names in its third dimension."
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

# The table of site effects. Each row summarizes the posterior distribution of
# one theta[i], which averages over the posterior draws of g (Lee and Sui
# 2025, Section 4.3 and Appendix A): its mean, its standard deviation, the 5%
# and 95% quantiles and its mode on the grid.
.bef_theta_summary <- function(posterior, log_g, stan_data) {
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
    sd = .bef_posterior_sd(posterior$theta_mean, posterior$theta_sd),
    hpdi_lower = hpdi[1L, ],
    hpdi_upper = hpdi[2L, ],
    map = .bef_posterior_mode(
      log_g,
      theta_hat = as.numeric(stan_data$theta_hat),
      sigma = as.numeric(stan_data$sigma),
      grid = as.numeric(stan_data$grid)
    ),
    row.names = seq_len(K),
    check.names = FALSE
  )
}

# The posterior standard deviation of each site effect. The posterior variance
# has two parts: the variance of theta[i] given g, averaged over the draws of
# g, and the variance over the draws of the mean of theta[i] given g.
# `theta_mean` and `theta_sd` have one row per draw and one column per site.
.bef_posterior_sd <- function(theta_mean, theta_sd) {
  within <- colMeans(theta_sd^2)
  n_draws <- nrow(theta_mean)
  if (n_draws <= 1L) {
    return(unname(sqrt(within)))
  }
  centered <- sweep(theta_mean, 2L, colMeans(theta_mean))
  between <- colSums(centered^2) / (n_draws - 1L)
  unname(sqrt(within + between))
}

# The posterior mode of each site effect: the grid point to which the
# posterior distribution of theta[i] gives the largest probability. With
# equal probabilities the smallest grid point is taken, as in the Stan
# program.
.bef_posterior_mode <- function(log_g, theta_hat, sigma, grid) {
  probability <- .bef_posterior_probabilities(log_g, theta_hat, sigma, grid)
  grid[apply(probability, 1L, which.max)]
}

# P(theta[i] = grid[j] | data) as a matrix with one row per site and one
# column per grid point: the posterior probability given g, averaged over the
# draws of g. `log_g` has one row per draw and one column per grid point.
# The sites are handled in blocks, so that no intermediate matrix has more
# than `max_cells` elements.
.bef_posterior_probabilities <- function(log_g, theta_hat, sigma, grid,
                                         max_cells = 2e6) {
  n_draws <- nrow(log_g)
  K <- length(theta_hat)
  L <- length(grid)
  g <- exp(log_g)

  # Log likelihood of each estimate at each grid point. Before it is
  # exponentiated, the largest value of each site is subtracted; that factor
  # cancels in the posterior probabilities.
  log_lik <- matrix(
    stats::dnorm(
      rep(theta_hat, times = L),
      mean = rep(grid, each = K),
      sd = rep(sigma, times = L),
      log = TRUE
    ),
    nrow = K,
    ncol = L
  )
  lik <- exp(log_lik - apply(log_lik, 1L, max))

  out <- matrix(NA_real_, nrow = K, ncol = L)
  per_block <- max(1L, as.integer(max_cells %/% n_draws))
  for (first in seq(1L, K, by = per_block)) {
    sites <- first:min(K, first + per_block - 1L)
    block <- lik[sites, , drop = FALSE]
    # marginal[s, i] = sum over j of g[s, j] * lik[i, j]
    marginal <- g %*% t(block)
    # mean over s of g[s, j] * lik[i, j] / marginal[s, i]
    out[sites, ] <- block * t(crossprod(g, 1 / marginal)) / n_draws
  }

  # If a product underflowed for some draw, redo that site on the log scale.
  for (site in which(rowSums(!is.finite(out)) > 0L)) {
    log_post <- sweep(log_g, 2L, log_lik[site, ], "+")
    weights <- exp(log_post - apply(log_post, 1L, max))
    out[site, ] <- colMeans(weights / rowSums(weights))
  }
  out
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

# A site-level generated quantity (theta_map, theta_mean, theta_sd or
# theta_rep) as a matrix with one row per draw and one column per site.
.bef_site_draws <- function(fit, name) {
  K <- as.integer(fit$metadata$data_list$K)
  .bef_draws_vector(.bef_draws_matrix(fit$draws), name, K)
}

# g or log_w as a matrix with one row per draw and one column per grid point.
# They are in the draws only when the model was fitted with
# store_grid_quantities = TRUE; otherwise they are computed from log_g and
# alpha.
.bef_grid_draws <- function(fit, which = c("g", "log_w")) {
  which <- match.arg(which)
  dm <- .bef_draws_matrix(fit$draws)
  L <- as.integer(fit$metadata$data_list$L)

  if (identical(which, "g")) {
    stored <- paste0("g[", seq_len(L), "]")
    if (all(stored %in% colnames(dm))) {
      return(dm[, stored, drop = FALSE])
    }
    return(exp(dm[, paste0("log_g[", seq_len(L), "]"), drop = FALSE]))
  }

  stored <- paste0("log_w[", seq_len(L), "]")
  if (all(stored %in% colnames(dm))) {
    return(dm[, stored, drop = FALSE])
  }
  # log_w = B %*% alpha, up to the additive constant that log_softmax removes.
  M <- as.integer(fit$metadata$data_list$M)
  B <- matrix(fit$metadata$data_list$B, nrow = L, ncol = M)
  alpha <- dm[, paste0("alpha[", seq_len(M), "]"), drop = FALSE]
  alpha %*% t(B)
}

# The log likelihood of each site as a matrix with one row per draw, the
# shape that the loo package expects.
.bef_log_lik_draws <- function(fit) {
  K <- as.integer(fit$metadata$data_list$K)
  dm <- .bef_draws_matrix(fit$draws)
  wanted <- paste0("log_lik[", seq_len(K), "]")
  if (!all(wanted %in% colnames(dm))) {
    .bef_abort_invalid_fit(
      paste0(
        "This fit has no `log_lik` draws, which loo() and waic() need. It was ",
        "created by bayesEfron 0.1; refit the model."
      )
    )
  }
  dm[, wanted, drop = FALSE]
}
