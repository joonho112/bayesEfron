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
  if (!metadata$grid_method %in% .bef_grid_methods()) {
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
  if (!.bef_is_whole_number(data_list$K) || data_list$K < 5L) {
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
  .bef_require_exact_fields(
    x,
    c("site", "mean", "sd", "hpdi_lower", "hpdi_upper", "map"),
    "`metadata$theta_summary`",
    class
  )
  if (!is.integer(x$site) && !is.numeric(x$site)) {
    .bef_abort_validate("`metadata$theta_summary$site` must be numeric.", class)
  }
  if (!identical(as.integer(x$site), seq_len(K))) {
    .bef_abort_validate(
      "`metadata$theta_summary$site` must enumerate sites from 1 to K.",
      class
    )
  }
  .bef_validate_numeric_column(x$mean, "`metadata$theta_summary$mean`", class)
  .bef_validate_numeric_column(x$sd, "`metadata$theta_summary$sd`", class, lower = 0)
  .bef_validate_numeric_column(x$hpdi_lower, "`metadata$theta_summary$hpdi_lower`", class)
  .bef_validate_numeric_column(x$hpdi_upper, "`metadata$theta_summary$hpdi_upper`", class)
  if (any(x$hpdi_lower > x$hpdi_upper)) {
    .bef_abort_validate(
      "`metadata$theta_summary$hpdi_lower` must be <= `hpdi_upper`.",
      class
    )
  }
  .bef_validate_numeric_column(x$map, "`metadata$theta_summary$map`", class)
}

# The table of site effects must agree with the draws that it summarizes:
# `mean` with the mean of `theta_mean`, `sd` with the posterior standard
# deviation computed from `theta_mean` and `theta_sd`, and `map` must be a
# grid point.
.bef_validate_theta_summary_consistency <- function(theta_summary,
                                                    theta_mean_draws,
                                                    theta_sd_draws,
                                                    grid,
                                                    class) {
  .bef_expect_close(
    theta_summary$mean,
    colMeans(theta_mean_draws),
    "`metadata$theta_summary$mean` must match `theta_mean` column means in `draws`.",
    class
  )
  .bef_expect_close(
    theta_summary$sd,
    .bef_posterior_sd(theta_mean_draws, theta_sd_draws),
    paste0(
      "`metadata$theta_summary$sd` must be the posterior standard deviation ",
      "computed from `theta_mean` and `theta_sd` in `draws`."
    ),
    class
  )
  grid <- as.numeric(grid)
  distance <- vapply(
    theta_summary$map,
    function(value) min(abs(value - grid)),
    numeric(1L)
  )
  if (any(distance > sqrt(.Machine$double.eps) * max(1, abs(grid)))) {
    .bef_abort_validate(
      "`metadata$theta_summary$map` must hold points of the grid.",
      class
    )
  }
}

.bef_expect_close <- function(x, y, message, class) {
  tolerance <- sqrt(.Machine$double.eps)
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
.bef_validate_site_draws_present <- function(draws, K, class) {
  available <- dimnames(draws)[[3L]]
  for (field in .bef_site_generated_quantity_fields()) {
    expected <- paste0(field, "[", seq_len(K), "]")
    if (!all(expected %in% available)) {
      .bef_abort_validate(
        sprintf("`bef_fit_re$draws` must contain %s[1:K].", field),
        class
      )
    }
  }
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
