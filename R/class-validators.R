# Checks of fields and values shared by the class validators.

.bef_require_inherits <- function(x, class_name, arg, class) {
  if (!inherits(x, class_name)) {
    .bef_abort_validate(
      sprintf("`%s` must inherit from class \"%s\".", arg, class_name),
      class
    )
  }
}

.bef_require_fields <- function(x, required, what, class) {
  missing <- setdiff(required, names(x))
  if (!is.list(x) || length(missing) > 0L) {
    .bef_abort_validate(
      sprintf(
        "%s must contain required fields: %s.",
        what,
        paste(required, collapse = ", ")
      ),
      class,
      missing_fields = missing
    )
  }
}

.bef_require_exact_fields <- function(x, required, what, class) {
  missing <- setdiff(required, names(x))
  extra <- setdiff(names(x), required)
  if (!is.list(x) || length(missing) > 0L || length(extra) > 0L ||
      !identical(names(x), required)) {
    .bef_abort_validate(
      sprintf(
        "%s must contain exactly these fields: %s.",
        what,
        paste(required, collapse = ", ")
      ),
      class,
      missing_fields = missing,
      extra_fields = extra
    )
  }
}

.bef_validate_draws_array <- function(x, class) {
  if (!is.array(x) || !is.numeric(x) || length(dim(x)) != 3L ||
      any(dim(x) <= 0L) ||
      any(!is.finite(x))) {
    .bef_abort_validate(
      "`bef_fit$draws` must be a finite numeric 3D draws array with positive dimensions.",
      class
    )
  }
}

.bef_validate_metadata_core <- function(metadata, class) {
  if (!identical(metadata$model_family, "RE")) {
    .bef_abort_validate(
      "`metadata$model_family` must be \"RE\"; other families are not implemented.",
      class
    )
  }
  # "support" is recorded when the end points of the grid were given.
  if (!metadata$grid_method %in% c(.bef_grid_methods(), "support")) {
    .bef_abort_validate(
      "`metadata$grid_method` must be one of the supported grid methods.",
      class
    )
  }
  if (!.bef_is_whole_number(metadata$seed)) {
    .bef_abort_validate("`metadata$seed` must be a finite integer scalar.", class)
  }
  if (!.bef_is_string(metadata$cmdstan_version)) {
    .bef_abort_validate("`metadata$cmdstan_version` must be a non-empty string.", class)
  }
  .bef_validate_sha256(metadata$stan_file_sha256, "`metadata$stan_file_sha256`", class)
  if (!.bef_is_number(metadata$runtime_seconds) ||
      metadata$runtime_seconds < 0) {
    .bef_abort_validate(
      "`metadata$runtime_seconds` must be a non-negative finite numeric scalar.",
      class
    )
  }
}

.bef_validate_stan_data_list <- function(data_list, class) {
  .bef_require_fields(
    data_list,
    c("K", "theta_hat", "sigma", "L", "grid", "M", "B", "store_grid_quantities"),
    "`metadata$data_list`", class
  )
  if (!.bef_is_whole_number(data_list$K) || data_list$K < .bef_min_sites()) {
    .bef_abort_validate("`metadata$data_list$K` must be an integer scalar >= 5.", class)
  }
  if (!.bef_is_whole_number(data_list$L) || data_list$L < 1L) {
    .bef_abort_validate("`metadata$data_list$L` must be a positive integer scalar.", class)
  }
  if (!.bef_is_whole_number(data_list$M) || data_list$M < 1L) {
    .bef_abort_validate("`metadata$data_list$M` must be a positive integer scalar.", class)
  }

  K <- as.integer(data_list$K)
  L <- as.integer(data_list$L)
  M <- as.integer(data_list$M)

  if (!is.numeric(data_list$theta_hat) ||
      length(data_list$theta_hat) != K ||
      any(!is.finite(data_list$theta_hat))) {
    .bef_abort_validate(
      "`metadata$data_list$theta_hat` must be a finite numeric vector of length K.",
      class
    )
  }
  if (!is.numeric(data_list$sigma) ||
      length(data_list$sigma) != K ||
      any(!is.finite(data_list$sigma)) ||
      any(data_list$sigma <= 0)) {
    .bef_abort_validate(
      "`metadata$data_list$sigma` must be a strictly positive finite numeric vector of length K.",
      class
    )
  }
  if (!is.numeric(data_list$grid) ||
      length(data_list$grid) != L ||
      any(!is.finite(data_list$grid)) ||
      any(diff(data_list$grid) <= 0)) {
    .bef_abort_validate(
      "`metadata$data_list$grid` must be a finite strictly increasing numeric vector of length L.",
      class
    )
  }
  if (!is.matrix(data_list$B) ||
      !is.numeric(data_list$B) ||
      nrow(data_list$B) != L ||
      ncol(data_list$B) != M ||
      any(!is.finite(data_list$B))) {
    .bef_abort_validate(
      "`metadata$data_list$B` must be a finite numeric matrix with dimensions L by M.",
      class
    )
  }
}

.bef_validate_summary_list <- function(x, field, class) {
  required <- c("mean", "sd", "q5", "q50", "q95")
  if (!is.list(x) || !identical(names(x), required)) {
    .bef_abort_validate(
      sprintf(
        "`%s` must be a summary list with exactly these fields: %s.",
        field,
        paste(required, collapse = ", ")
      ),
      class
    )
  }
  valid <- vapply(
    x,
    function(value) is.numeric(value) && length(value) == 1L && is.finite(value),
    logical(1)
  )
  if (!all(valid)) {
    .bef_abort_validate(
      sprintf("All fields in `%s` must be finite numeric scalars.", field),
      class
    )
  }
}

.bef_validate_postprocess_metadata_attrs <- function(metadata, class) {
  .bef_validate_summary_list(
    attr(metadata, "sd_g_summary", exact = TRUE),
    "attr(metadata, \"sd_g_summary\")",
    class
  )

  diagnostics <- attr(metadata, "diagnostics", exact = TRUE)
  .bef_require_exact_fields(
    diagnostics,
    c("rhat", "ess_bulk", "ess_tail", "divergences", "max_treedepth", "ebfmi"),
    "attr(metadata, \"diagnostics\")",
    class
  )
  .bef_validate_diagnostic_numeric(
    diagnostics$rhat, "rhat", class, lower = 0, open_lower = TRUE
  )
  .bef_validate_diagnostic_numeric(
    diagnostics$ess_bulk, "ess_bulk", class, lower = 0
  )
  .bef_validate_diagnostic_numeric(
    diagnostics$ess_tail, "ess_tail", class, lower = 0
  )
  .bef_validate_diagnostic_integerish(
    diagnostics$divergences, "divergences", class, lower = 0
  )
  .bef_validate_diagnostic_integerish(
    diagnostics$max_treedepth, "max_treedepth", class, lower = 0
  )
  .bef_validate_diagnostic_numeric(
    diagnostics$ebfmi, "ebfmi", class, lower = 0
  )

  diagnostic_skipped <- attr(metadata, "diagnostic_skipped", exact = TRUE)
  if (!is.character(diagnostic_skipped) || anyNA(diagnostic_skipped)) {
    .bef_abort_validate(
      "attr(metadata, \"diagnostic_skipped\") must be a character vector without missing values.",
      class
    )
  }
  if (any(duplicated(diagnostic_skipped)) ||
      length(setdiff(diagnostic_skipped, .bef_diagnostic_skipped_fields())) > 0L) {
    .bef_abort_validate(
      "attr(metadata, \"diagnostic_skipped\") contains unsupported diagnostic names.",
      class
    )
  }
  .bef_validate_diagnostic_skip_consistency(diagnostics, diagnostic_skipped, class)

  sampler_failed <- attr(metadata, "sampler_diagnostics_failed", exact = TRUE)
  if (!is.character(sampler_failed) || anyNA(sampler_failed)) {
    .bef_abort_validate(
      "attr(metadata, \"sampler_diagnostics_failed\") must be a character vector without missing values.",
      class
    )
  }
  if (any(duplicated(sampler_failed)) ||
      length(setdiff(sampler_failed, .bef_sampler_diagnostic_failure_fields())) > 0L ||
      length(intersect(sampler_failed, diagnostic_skipped)) > 0L) {
    .bef_abort_validate(
      "attr(metadata, \"sampler_diagnostics_failed\") contains unsupported or skipped diagnostic names.",
      class
    )
  }

  sampler_warned <- attr(metadata, "sampler_diagnostics_warned", exact = TRUE)
  if (!is.character(sampler_warned) || anyNA(sampler_warned)) {
    .bef_abort_validate(
      "attr(metadata, \"sampler_diagnostics_warned\") must be a character vector without missing values.",
      class
    )
  }
  if (any(duplicated(sampler_warned)) ||
      length(setdiff(sampler_warned, .bef_sampler_diagnostic_warning_fields())) > 0L) {
    .bef_abort_validate(
      "attr(metadata, \"sampler_diagnostics_warned\") contains unsupported diagnostic names.",
      class
    )
  }
}

.bef_validate_theta_summary <- function(x, K, class) {
  if (!is.data.frame(x) || nrow(x) != K) {
    .bef_abort_validate(
      "`metadata$theta_summary` must be a data frame with one row per site.",
      class
    )
  }
  # The limits and the mode may carry the names of this version, those of
  # version 0.3.0, or both.
  new_names <- c("lower", "upper", "mode")
  old_names <- c("hpdi_lower", "hpdi_upper", "map")
  if (!all(new_names %in% names(x)) && !all(old_names %in% names(x))) {
    .bef_abort_validate(
      paste0(
        "`metadata$theta_summary` must have the columns `site`, `mean`, `sd`, ",
        "`lower`, `upper` and `mode`."
      ),
      class,
      missing_fields = setdiff(new_names, names(x))
    )
  }
  .bef_require_fields(x, c("site", "mean", "sd"), "`metadata$theta_summary`", class)
  unknown <- setdiff(names(x), .bef_site_table_columns())
  if (length(unknown) > 0L) {
    .bef_abort_validate(
      sprintf(
        "`metadata$theta_summary` has unknown columns: %s.",
        paste(sprintf("`%s`", unknown), collapse = ", ")
      ),
      class,
      unknown_fields = unknown
    )
  }
  if (!is.integer(x$site) && !is.numeric(x$site)) {
    .bef_abort_validate("`metadata$theta_summary$site` must be numeric.", class)
  }
  if (!identical(as.integer(x$site), seq_len(K))) {
    .bef_abort_validate(
      "`metadata$theta_summary$site` must enumerate sites from 1 to K.",
      class
    )
  }
  if (!is.null(x$label) &&
      (!is.character(x$label) || anyNA(x$label) || anyDuplicated(x$label) > 0L)) {
    .bef_abort_validate(
      "`metadata$theta_summary$label` must hold one distinct label per site.",
      class
    )
  }
  .bef_validate_numeric_column(x$mean, "`metadata$theta_summary$mean`", class)
  .bef_validate_numeric_column(x$sd, "`metadata$theta_summary$sd`", class, lower = 0)
  lower <- if (is.null(x$lower)) x$hpdi_lower else x$lower
  upper <- if (is.null(x$upper)) x$hpdi_upper else x$upper
  mode <- if (is.null(x$mode)) x$map else x$mode
  .bef_validate_numeric_column(lower, "`metadata$theta_summary$lower`", class)
  .bef_validate_numeric_column(upper, "`metadata$theta_summary$upper`", class)
  if (any(lower > upper)) {
    .bef_abort_validate(
      "`metadata$theta_summary$lower` must be <= `upper`.",
      class
    )
  }
  .bef_validate_numeric_column(mode, "`metadata$theta_summary$mode`", class)
  for (pair in list(c("lower", "hpdi_lower"), c("upper", "hpdi_upper"), c("mode", "map"))) {
    if (all(pair %in% names(x)) && !isTRUE(all.equal(x[[pair[1L]]], x[[pair[2L]]]))) {
      .bef_abort_validate(
        sprintf("`metadata$theta_summary$%s` and `$%s` must hold the same values.",
                pair[1L], pair[2L]),
        class
      )
    }
  }
}

# The table of site effects must agree with the draws that it summarizes:
# `mean` and `sd` with the mean and the standard deviation of the posterior
# distribution computed from `log_g`, and, with point masses, `mode` must be
# a grid point. The moments are computed again for a few sites only, so that
# the check stays cheap for a fit with many sites.
.bef_validate_theta_summary_consistency <- function(fit, class) {
  table <- .bef_site_table(fit)
  data <- fit$metadata$data_list
  grid <- as.numeric(data$grid)
  scale <- max(1, abs(grid))
  tau_c <- if (is.null(data$tau_c)) 0 else data$tau_c
  if (tau_c == 0) {
    distance <- vapply(
      table$mode,
      function(value) min(abs(value - grid)),
      numeric(1L)
    )
    if (any(distance > sqrt(.Machine$double.eps) * scale)) {
      .bef_abort_validate(
        "`metadata$theta_summary$mode` must hold points of the grid.",
        class
      )
    }
  }
  # The moments can be checked only against the draws of g, and only for a
  # table with the definitions of this version: one saved by an earlier
  # version is computed again when a method uses the fit.
  variables <- dimnames(fit$draws)[[3L]]
  current <- identical(
    attr(fit$metadata, "summary_definition_version", exact = TRUE),
    .bef_summary_version()
  )
  if (!current || !any(c("log_g[1]", "g[1]") %in% variables)) {
    return(invisible(TRUE))
  }
  engine <- .bef_engine(fit)
  K <- engine$K
  check <- unique(round(seq(1, K, length.out = min(K, 5L))))
  moments <- .bef_site_moments(engine, sites = check)
  .bef_expect_close(
    table$mean[check] / scale,
    moments$mean / scale,
    "`metadata$theta_summary$mean` must match the posterior means computed from `log_g` in `draws`.",
    class
  )
  .bef_expect_close(
    table$sd[check] / scale,
    sqrt(moments$var) / scale,
    paste0(
      "`metadata$theta_summary$sd` must be the posterior standard deviation ",
      "computed from `log_g` in `draws`."
    ),
    class
  )
  invisible(TRUE)
}

.bef_expect_close <- function(x, y, message, class) {
  # Looser than the precision of the arithmetic: a table computed from the
  # values that CmdStan wrote with eight significant digits must pass.
  tolerance <- 1e-6
  if (length(x) != length(y) || any(abs(x - y) > tolerance)) {
    .bef_abort_validate(message, class)
  }
}

.bef_validate_generated_quantities <- function(posterior, K, n_draws, class) {
  for (field in .bef_generated_quantity_fields()) {
    value <- posterior[[field]]
    if (!is.numeric(value) || any(!is.finite(value))) {
      .bef_abort_validate(
        sprintf("`bef_fit_re$posterior$%s` must be finite and numeric.", field),
        class
      )
    }
    if (length(value) != n_draws) {
      .bef_abort_validate(
        sprintf("`bef_fit_re$posterior$%s` must have one value per draw.", field),
        class
      )
    }
  }
}

# The site-level generated quantities must be in the draws array.
# The site-level draws are optional. A fit that keeps any of them must keep
# each kept block whole, one variable for each site.
.bef_validate_site_draws <- function(draws, K, class) {
  available <- dimnames(draws)[[3L]]
  blocks <- sub("\\[.*$", "", available)
  for (field in c(.bef_site_generated_quantity_fields(), "log_lik")) {
    if (!field %in% blocks) next
    expected <- paste0(field, "[", seq_len(K), "]")
    if (!all(expected %in% available)) {
      .bef_abort_validate(
        sprintf("`bef_fit_re$draws` must contain %s[1:K] or none of it.", field),
        class
      )
    }
  }
}

# Does the fit keep the site-level draws of version 0.3.0?
.bef_has_site_draws <- function(fit) {
  K <- as.integer(fit$metadata$data_list$K)
  wanted <- unlist(lapply(
    .bef_site_generated_quantity_fields(),
    function(field) paste0(field, "[", c(1L, K), "]")
  ))
  all(wanted %in% dimnames(fit$draws)[[3L]])
}

.bef_validate_numeric_column <- function(x, field, class, lower = -Inf) {
  if (!is.numeric(x) || any(!is.finite(x)) || any(x < lower)) {
    .bef_abort_validate(
      sprintf("%s must be finite and numeric.", field),
      class
    )
  }
}

.bef_validate_sha256 <- function(x, field, class) {
  if (!.bef_is_string(x) || !grepl("^[0-9a-fA-F]{64}$", x)) {
    .bef_abort_validate(
      sprintf("%s must be a 64-character hexadecimal SHA-256 string.", field),
      class
    )
  }
}

.bef_is_string <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

.bef_is_number <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x)
}

.bef_is_whole_number <- function(x) {
  .bef_is_number(x) && x == as.integer(x)
}
