condition_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

condition_expected_classes <- function() {
  list(
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
    bef_family_contract_violation = c(
      "bef_family_contract_violation", "bef_firewall_error", "bef_error"
    ),
    bef_err_grid_oracle_required = c(
      "bef_err_grid_oracle_required", "bef_grid_error", "bef_error"
    ),
    bef_grid_rank_deficient = c("bef_grid_rank_deficient", "bef_grid_error", "bef_error"),
    bef_diagnostic_skipped = c("bef_diagnostic_skipped", "bef_warning"),
    bef_sampler_diagnostics_failed = c("bef_sampler_diagnostics_failed", "bef_warning"),
    bef_cache_corruption = c("bef_cache_corruption", "bef_cache_warning", "bef_warning"),
    bef_method_unsupported = c(
      "bef_method_unsupported", "bef_unsupported_method", "bef_warning"
    ),
    bef_family_contract_unexpected_field = c(
      "bef_family_contract_unexpected_field", "bef_firewall_warning", "bef_warning"
    )
  )
}

condition_capture_error <- function(expr) {
  tryCatch(force(expr), error = identity)
}

condition_capture_warning <- function(expr) {
  captured <- NULL
  withCallingHandlers(
    force(expr),
    warning = function(w) {
      captured <<- w
      invokeRestart("muffleWarning")
    }
  )
  captured
}

expect_condition_prefix <- function(cnd, expected) {
  expect_identical(class(cnd)[seq_along(expected)], expected)
}

test_that("typed condition catalog locks every condition class chain", {
  expected <- condition_expected_classes()

  expect_identical(condition_ns(".bef_condition_catalog")(), names(expected))
  for (class in names(expected)) {
    expect_identical(
      condition_ns(".bef_condition_classes")(class),
      expected[[class]],
      info = class
    )
  }

  expect_identical(
    condition_ns(".bef_condition_classes")(
      "bef_invalid_args",
      extra_class = "bayesEfron_validate_error"
    ),
    c(
      "bef_invalid_args",
      "bayesEfron_validate_error",
      "bef_pipeline_error",
      "bef_error"
    )
  )

  err <- condition_capture_error(
    condition_ns(".bef_condition_classes")("bef_not_real")
  )
  expect_s3_class(err, "rlang_error")
  expect_match(conditionMessage(err), "Unknown bayesEfron condition class")
})

test_that("abort helpers emit their documented typed conditions", {
  abort_helpers <- list(
    bef_invalid_args = ".bef_abort_invalid_args",
    bef_compile_failed = ".bef_abort_compile_failed",
    bef_sampling_failed = ".bef_abort_sampling_failed",
    bef_sampling_partial = ".bef_abort_sampling_partial",
    bef_sampling_unusable = ".bef_abort_sampling_unusable",
    bef_extraction_failed = ".bef_abort_extraction_failed",
    bef_invalid_fit = ".bef_abort_invalid_fit",
    bef_cache_perm_violation = ".bef_abort_cache_perm_violation",
    bef_lock_timeout = ".bef_abort_lock_timeout",
    bef_cache_format_mismatch = ".bef_abort_cache_format_mismatch",
    bef_cache_version_drift = ".bef_abort_cache_version_drift",
    bef_family_contract_violation = ".bef_abort_family_contract_violation",
    bef_err_grid_oracle_required = ".bef_abort_grid_oracle_required",
    bef_grid_rank_deficient = ".bef_abort_grid_rank_deficient"
  )

  expected <- condition_expected_classes()
  for (class in names(abort_helpers)) {
    helper <- abort_helpers[[class]]
    err <- condition_capture_error(
      do.call(
        condition_ns(helper),
        list(
          message = "condition probe",
          condition_probe = helper,
          required_arg = "theta_true"
        )
      )
    )
    expected_prefix <- expected[[class]]
    if (identical(class, "bef_invalid_fit")) {
      expected_prefix <- c(
        "bef_invalid_fit",
        "bayesEfron_validate_error",
        "bef_pipeline_error",
        "bef_error"
      )
    }

    expect_condition_prefix(err, expected_prefix)
    expect_equal(err$condition_probe, helper)
    if (identical(class, "bef_err_grid_oracle_required")) {
      expect_equal(err$required_arg, "theta_true")
    }
  }

  err <- condition_capture_error(
    condition_ns(".bef_abort_invalid_args")(
      "validation probe",
      validate = TRUE,
      condition_probe = "validate"
    )
  )
  expect_condition_prefix(
    err,
    c(
      "bef_invalid_args",
      "bayesEfron_validate_error",
      "bef_pipeline_error",
      "bef_error"
    )
  )
  expect_equal(err$condition_probe, "validate")

  parent <- simpleError("parent probe")
  err <- condition_capture_error(
    condition_ns(".bef_abort_compile_failed")("compile probe", parent = parent)
  )
  expect_condition_prefix(err, expected$bef_compile_failed)
  expect_identical(err$parent, parent)
})

test_that("warning helpers emit their documented typed conditions", {
  warning_helpers <- list(
    bef_diagnostic_skipped = ".bef_warn_diagnostic_skipped",
    bef_sampler_diagnostics_failed = ".bef_warn_sampler_diagnostics_failed",
    bef_cache_corruption = ".bef_warn_cache_corruption",
    bef_method_unsupported = ".bef_warn_method_unsupported",
    bef_family_contract_unexpected_field = ".bef_warn_family_contract_unexpected_field"
  )

  expected <- condition_expected_classes()
  for (class in names(warning_helpers)) {
    helper <- warning_helpers[[class]]
    warning <- condition_capture_warning(
      do.call(
        condition_ns(helper),
        list(message = "condition probe", condition_probe = helper)
      )
    )

    expect_condition_prefix(warning, expected[[class]])
    expect_equal(warning$condition_probe, helper)
  }
})
