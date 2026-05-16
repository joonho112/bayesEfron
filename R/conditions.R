.bef_abort <- function(message,
                       class,
                       ...,
                       parent = NULL,
                       extra_class = NULL,
                       call = rlang::caller_env()) {
  rlang::abort(
    message,
    class = .bef_condition_classes(class, extra_class = extra_class),
    ...,
    parent = parent,
    call = call
  )
}

.bef_warn <- function(message,
                      class,
                      ...,
                      parent = NULL,
                      call = rlang::caller_env()) {
  rlang::warn(
    message,
    class = .bef_condition_classes(class),
    ...,
    parent = parent,
    call = call
  )
}

.bef_condition_classes <- function(class, extra_class = NULL) {
  classes <- switch(
    class,
    bef_invalid_args = c("bef_invalid_args", "bef_pipeline_error", "bef_error"),
    bef_compile_failed = c("bef_compile_failed", "bef_pipeline_error", "bef_error"),
    bef_sampling_failed = c("bef_sampling_failed", "bef_pipeline_error", "bef_error"),
    bef_sampling_partial = c("bef_sampling_partial", "bef_pipeline_error", "bef_error"),
    bef_sampling_unusable = c("bef_sampling_unusable", "bef_pipeline_error", "bef_error"),
    bef_extraction_failed = c("bef_extraction_failed", "bef_pipeline_error", "bef_error"),
    bef_invalid_fit = c("bef_invalid_fit", "bef_pipeline_error", "bef_error"),
    bef_cache_perm_violation = c("bef_cache_perm_violation", "bef_cache_error", "bef_error"),
    bef_lock_timeout = c("bef_lock_timeout", "bef_cache_error", "bef_error"),
    bef_cache_format_mismatch = c("bef_cache_format_mismatch", "bef_cache_error", "bef_error"),
    bef_cache_version_drift = c("bef_cache_version_drift", "bef_cache_error", "bef_error"),
    bef_family_contract_violation = c("bef_family_contract_violation", "bef_firewall_error", "bef_error"),
    bef_err_grid_oracle_required = c("bef_err_grid_oracle_required", "bef_grid_error", "bef_error"),
    bef_grid_rank_deficient = c("bef_grid_rank_deficient", "bef_grid_error", "bef_error"),
    bef_diagnostic_skipped = c("bef_diagnostic_skipped", "bef_warning"),
    bef_sampler_diagnostics_failed = c("bef_sampler_diagnostics_failed", "bef_warning"),
    bef_cache_corruption = c("bef_cache_corruption", "bef_cache_warning", "bef_warning"),
    bef_method_unsupported = c("bef_method_unsupported", "bef_unsupported_method", "bef_warning"),
    bef_family_contract_unexpected_field = c("bef_family_contract_unexpected_field", "bef_firewall_warning", "bef_warning"),
    NULL
  )

  if (is.null(classes)) {
    rlang::abort(sprintf("Unknown bayesEfron condition class: %s", class))
  }
  if (!is.null(extra_class)) {
    classes <- c(classes[1L], extra_class, classes[-1L])
  }
  classes
}

.bef_condition_catalog <- function() {
  c(
    "bef_invalid_args",
    "bef_compile_failed",
    "bef_sampling_failed",
    "bef_sampling_partial",
    "bef_sampling_unusable",
    "bef_extraction_failed",
    "bef_invalid_fit",
    "bef_cache_perm_violation",
    "bef_lock_timeout",
    "bef_cache_format_mismatch",
    "bef_cache_version_drift",
    "bef_family_contract_violation",
    "bef_err_grid_oracle_required",
    "bef_grid_rank_deficient",
    "bef_diagnostic_skipped",
    "bef_sampler_diagnostics_failed",
    "bef_cache_corruption",
    "bef_method_unsupported",
    "bef_family_contract_unexpected_field"
  )
}

.bef_abort_invalid_args <- function(message, ..., parent = NULL, validate = FALSE) {
  extra_class <- if (isTRUE(validate)) "bayesEfron_validate_error" else NULL
  .bef_abort(
    message,
    "bef_invalid_args",
    ...,
    parent = parent,
    extra_class = extra_class
  )
}

.bef_abort_compile_failed <- function(message, ..., parent = NULL) {
  .bef_abort(message, "bef_compile_failed", ..., parent = parent)
}

.bef_abort_sampling_failed <- function(message, ..., parent = NULL) {
  .bef_abort(message, "bef_sampling_failed", ..., parent = parent)
}

.bef_abort_sampling_partial <- function(message, ..., parent = NULL) {
  .bef_abort(message, "bef_sampling_partial", ..., parent = parent)
}

.bef_abort_sampling_unusable <- function(message, ..., parent = NULL) {
  .bef_abort(message, "bef_sampling_unusable", ..., parent = parent)
}

.bef_abort_extraction_failed <- function(message, ..., parent = NULL) {
  .bef_abort(message, "bef_extraction_failed", ..., parent = parent)
}

.bef_abort_invalid_fit <- function(message, ..., parent = NULL) {
  .bef_abort(
    message,
    "bef_invalid_fit",
    ...,
    parent = parent,
    extra_class = "bayesEfron_validate_error"
  )
}

.bef_abort_cache_perm_violation <- function(message, ..., parent = NULL) {
  .bef_abort(message, "bef_cache_perm_violation", ..., parent = parent)
}

.bef_abort_lock_timeout <- function(message, ..., parent = NULL) {
  .bef_abort(message, "bef_lock_timeout", ..., parent = parent)
}

.bef_abort_cache_format_mismatch <- function(message, ..., parent = NULL) {
  .bef_abort(message, "bef_cache_format_mismatch", ..., parent = parent)
}

.bef_abort_cache_version_drift <- function(message, ..., parent = NULL) {
  .bef_abort(message, "bef_cache_version_drift", ..., parent = parent)
}

.bef_abort_family_contract_violation <- function(message, ..., parent = NULL) {
  .bef_abort(message, "bef_family_contract_violation", ..., parent = parent)
}

.bef_abort_grid_oracle_required <- function(message,
                                           ...,
                                           parent = NULL,
                                           required_arg = "theta_true") {
  .bef_abort(
    message,
    "bef_err_grid_oracle_required",
    ...,
    parent = parent,
    required_arg = required_arg
  )
}

.bef_abort_grid_rank_deficient <- function(message, ..., parent = NULL) {
  .bef_abort(message, "bef_grid_rank_deficient", ..., parent = parent)
}

.bef_warn_diagnostic_skipped <- function(message, ..., parent = NULL) {
  .bef_warn(message, "bef_diagnostic_skipped", ..., parent = parent)
}

.bef_warn_sampler_diagnostics_failed <- function(message, ..., parent = NULL) {
  .bef_warn(message, "bef_sampler_diagnostics_failed", ..., parent = parent)
}

.bef_warn_cache_corruption <- function(message, ..., parent = NULL) {
  .bef_warn(message, "bef_cache_corruption", ..., parent = parent)
}

.bef_warn_method_unsupported <- function(message, ..., parent = NULL) {
  .bef_warn(message, "bef_method_unsupported", ..., parent = parent)
}

.bef_warn_family_contract_unexpected_field <- function(message, ..., parent = NULL) {
  .bef_warn(message, "bef_family_contract_unexpected_field", ..., parent = parent)
}
