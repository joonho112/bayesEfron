new_bef_data <- function(theta_hat, sigma, names = NULL, source = "user") {
  structure(
    list(
      theta_hat = theta_hat,
      sigma = sigma,
      names = names,
      source = source
    ),
    class = "bef_data"
  )
}

validate_bef_data <- function(x) {
  class <- "bef_invalid_args"
  .bef_require_inherits(x, "bef_data", "x", class)
  .bef_require_fields(x, .bef_bef_data_fields(), "`bef_data`", class)

  if (!is.numeric(x$theta_hat) ||
      length(x$theta_hat) < 5L ||
      any(!is.finite(x$theta_hat))) {
    .bef_abort_validate(
      "`bef_data$theta_hat` must be a finite numeric vector of length at least 5.",
      class
    )
  }

  K <- length(x$theta_hat)
  if (!is.numeric(x$sigma) ||
      length(x$sigma) != K ||
      any(!is.finite(x$sigma)) ||
      any(x$sigma <= 0)) {
    .bef_abort_validate(
      "`bef_data$sigma` must be a strictly positive finite numeric vector with the same length as `theta_hat`.",
      class
    )
  }

  if (!is.null(x$names) &&
      (!is.character(x$names) || length(x$names) != K || anyNA(x$names))) {
    .bef_abort_validate(
      "`bef_data$names` must be NULL or a character vector with one non-missing label per site.",
      class
    )
  }

  if (!.bef_is_string(x$source)) {
    .bef_abort_validate(
      "`bef_data$source` must be a non-empty character scalar.",
      class
    )
  }

  x
}

new_bef_fit <- function(draws, metadata, posterior = list(), cmdstan_fit = NULL) {
  .bef_fit_new(
    draws = draws,
    metadata = metadata,
    posterior = posterior,
    cmdstan_fit = cmdstan_fit,
    class = "bef_fit"
  )
}

validate_bef_fit <- function(x) {
  class <- "bef_invalid_fit"
  .bef_require_inherits(x, "bef_fit", "x", class)
  .bef_require_fields(x, .bef_fit_fields(), "`bef_fit`", class)

  .bef_validate_draws_array(x$draws, class)
  if (!is.list(x$metadata)) {
    .bef_abort_validate("`bef_fit$metadata` must be a list.", class)
  }
  if (!is.list(x$posterior)) {
    .bef_abort_validate("`bef_fit$posterior` must be a list.", class)
  }

  .bef_require_fields(
    x$metadata, .bef_universal_metadata_fields(), "`bef_fit$metadata`", class
  )
  .bef_validate_metadata_core(x$metadata, class)
  .bef_validate_stan_data_list(x$metadata$data_list, class)
  .bef_validate_summary_list(x$metadata$mean_g_summary, "mean_g_summary", class)
  .bef_validate_summary_list(x$metadata$var_g_summary, "var_g_summary", class)
  .bef_validate_summary_list(
    x$metadata$effective_params_summary,
    "effective_params_summary",
    class
  )
  .bef_validate_summary_list(
    x$metadata$log_marginal_likelihood_summary,
    "log_marginal_likelihood_summary",
    class
  )

  x
}

new_bef_fit_re <- function(draws,
                           metadata,
                           posterior = list(),
                           cmdstan_fit = NULL) {
  .bef_fit_new(
    draws = draws,
    metadata = metadata,
    posterior = posterior,
    cmdstan_fit = cmdstan_fit,
    class = c("bef_fit_re", "bef_fit")
  )
}

validate_bef_fit_re <- function(x) {
  class <- "bef_invalid_fit"
  .bef_require_inherits(x, "bef_fit_re", "x", class)
  if (!identical(class(x), c("bef_fit_re", "bef_fit"))) {
    .bef_abort_validate(
      "`bef_fit_re` objects must have class vector c(\"bef_fit_re\", \"bef_fit\").",
      class
    )
  }

  validate_bef_fit(x)
  .bef_require_exact_fields(
    x$metadata, .bef_fit_re_metadata_fields(), "`bef_fit_re$metadata`", class
  )
  .bef_require_exact_fields(
    x$posterior, .bef_generated_quantity_fields(), "`bef_fit_re$posterior`",
    class
  )

  if (!identical(x$metadata$model_family, "RE")) {
    .bef_abort_validate(
      "`bef_fit_re$metadata$model_family` must be \"RE\".",
      class
    )
  }

  K <- as.integer(x$metadata$data_list$K)
  n_draws <- prod(dim(x$draws)[seq_len(2L)])
  .bef_validate_postprocess_metadata_attrs(x$metadata, class)
  .bef_validate_theta_summary(x$metadata$theta_summary, K, class)
  .bef_validate_theta_rep_draws(x$metadata$theta_rep_draws, K, n_draws, class)
  .bef_validate_generated_quantities(x$posterior, K, n_draws, class)
  if (!identical(x$metadata$theta_rep_draws, x$posterior$theta_rep)) {
    .bef_abort_validate(
      "`metadata$theta_rep_draws` must be identical to `posterior$theta_rep`.",
      class
    )
  }
  .bef_validate_theta_summary_consistency(x$metadata$theta_summary, x$posterior, class)

  x
}

.bef_fit_new <- function(draws,
                         metadata,
                         posterior = list(),
                         cmdstan_fit = NULL,
                         class = "bef_fit") {
  out <- list(
    draws = draws,
    metadata = metadata,
    posterior = posterior
  )

  if (!is.null(cmdstan_fit)) {
    out$cmdstan_fit <- cmdstan_fit
  }

  structure(out, class = class)
}

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
                               sampler_diagnostics_failed = character()) {
  structure(
    list(
      rhat = rhat,
      ess_bulk = ess_bulk,
      ess_tail = ess_tail,
      divergences = divergences,
      max_treedepth = max_treedepth,
      effective_params_summary = effective_params_summary,
      model_family = model_family,
      stan_file_sha256 = stan_file_sha256,
      runtime_seconds = runtime_seconds,
      diagnostic_skipped = diagnostic_skipped,
      sampler_diagnostics_failed = sampler_diagnostics_failed
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

  if (!identical(x$model_family, "RE")) {
    .bef_abort_validate(
      "`bef_diagnostic$model_family` must be \"RE\" for bayesEfron v0.1.",
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

  x
}

.bef_abort_validate <- function(message, class, ...) {
  if (identical(class, "bef_invalid_args")) {
    .bef_abort_invalid_args(message, ..., validate = TRUE)
  } else if (identical(class, "bef_invalid_fit")) {
    .bef_abort_invalid_fit(message, ...)
  } else {
    .bef_abort(
      message,
      class,
      ...,
      extra_class = "bayesEfron_validate_error"
    )
  }
}

.bef_bef_data_fields <- function() {
  c("theta_hat", "sigma", "names", "source")
}

.bef_fit_fields <- function() {
  c("draws", "metadata", "posterior")
}

.bef_universal_metadata_fields <- function() {
  c(
    "model_family",
    "grid_method",
    "seed",
    "cmdstan_version",
    "stan_file_sha256",
    "data_list",
    "runtime_seconds",
    "mean_g_summary",
    "var_g_summary",
    "effective_params_summary",
    "log_marginal_likelihood_summary"
  )
}

.bef_fit_re_metadata_fields <- function() {
  c(
    "model_family",
    "grid_method",
    "seed",
    "cmdstan_version",
    "stan_file_sha256",
    "data_list",
    "runtime_seconds",
    "mean_g_summary",
    "var_g_summary",
    "theta_summary",
    "theta_rep_draws",
    "effective_params_summary",
    "log_marginal_likelihood_summary"
  )
}

.bef_generated_quantity_fields <- function() {
  c(
    "mean_g",
    "var_g",
    "sd_g",
    "theta_map",
    "theta_mean",
    "theta_sd",
    "theta_rep",
    "effective_params",
    "log_marginal_likelihood"
  )
}

.bef_diagnostic_fields <- function() {
  c(
    "rhat",
    "ess_bulk",
    "ess_tail",
    "divergences",
    "max_treedepth",
    "effective_params_summary",
    "model_family",
    "stan_file_sha256",
    "runtime_seconds",
    "diagnostic_skipped",
    "sampler_diagnostics_failed"
  )
}

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
      "`metadata$model_family` must be \"RE\" for bayesEfron v0.1.",
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
    data_list, c("K", "theta_hat", "sigma", "L", "grid", "M", "B"),
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
    c("rhat", "ess_bulk", "ess_tail", "divergences", "max_treedepth"),
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
}

.bef_diagnostic_skipped_fields <- function() {
  c("rhat", "ess_bulk", "ess_tail", "sampler_diagnostics")
}

.bef_sampler_diagnostic_failure_fields <- function() {
  c("rhat", "ess_bulk", "ess_tail", "divergences")
}

.bef_validate_diagnostic_skip_consistency <- function(diagnostics,
                                                      diagnostic_skipped,
                                                      class) {
  scalar_skipped <- intersect(diagnostic_skipped, c("rhat", "ess_bulk", "ess_tail"))
  for (field in scalar_skipped) {
    if (!any(is.na(diagnostics[[field]]))) {
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
        function(field) any(is.na(diagnostics[[field]])),
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

.bef_validate_theta_summary_consistency <- function(theta_summary, posterior, class) {
  .bef_expect_close(
    theta_summary$mean,
    colMeans(posterior$theta_mean),
    "`metadata$theta_summary$mean` must match posterior `theta_mean` column means.",
    class
  )
  .bef_expect_close(
    theta_summary$sd,
    colMeans(posterior$theta_sd),
    "`metadata$theta_summary$sd` must match posterior `theta_sd` column means.",
    class
  )
  .bef_expect_close(
    theta_summary$map,
    colMeans(posterior$theta_map),
    "`metadata$theta_summary$map` must match posterior `theta_map` column means.",
    class
  )
}

.bef_expect_close <- function(x, y, message, class) {
  tolerance <- sqrt(.Machine$double.eps)
  if (length(x) != length(y) || any(abs(x - y) > tolerance)) {
    .bef_abort_validate(message, class)
  }
}

.bef_validate_theta_rep_draws <- function(x, K, n_draws, class) {
  if (!is.matrix(x) || !is.numeric(x) || ncol(x) != K ||
      nrow(x) != n_draws || any(!is.finite(x))) {
    .bef_abort_validate(
      "`metadata$theta_rep_draws` must be a finite numeric matrix with one row per draw and one column per site.",
      class
    )
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
  }

  for (field in c("theta_map", "theta_mean", "theta_sd", "theta_rep")) {
    value <- posterior[[field]]
    if (!is.matrix(value) || ncol(value) != K || nrow(value) != n_draws) {
      .bef_abort_validate(
        sprintf("`bef_fit_re$posterior$%s` must have one row per draw and one column per site.", field),
        class
      )
    }
  }
  for (field in setdiff(
    .bef_generated_quantity_fields(),
    c("theta_map", "theta_mean", "theta_sd", "theta_rep")
  )) {
    if (length(posterior[[field]]) != n_draws) {
      .bef_abort_validate(
        sprintf("`bef_fit_re$posterior$%s` must have one value per draw.", field),
        class
      )
    }
  }
  if (any(posterior$theta_sd < 0)) {
    .bef_abort_validate("`bef_fit_re$posterior$theta_sd` must be non-negative.", class)
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
