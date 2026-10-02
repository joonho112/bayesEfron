# The classes of the errors and warnings the package raises. Every error
# inherits from "bef_error" and every warning from "bef_warning", so a caller
# can catch all of them with one handler.

condition_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

condition_expected_classes <- function() {
  list(
    bef_invalid_args = c("bef_invalid_args", "bef_pipeline_error", "bef_error"),
    bef_compile_failed = c("bef_compile_failed", "bef_pipeline_error", "bef_error"),
    bef_sampling_failed = c("bef_sampling_failed", "bef_pipeline_error", "bef_error"),
    bef_sampling_partial = c("bef_sampling_partial", "bef_pipeline_error", "bef_error"),
    bef_extraction_failed = c("bef_extraction_failed", "bef_pipeline_error", "bef_error"),
    bef_invalid_fit = c("bef_invalid_fit", "bef_pipeline_error", "bef_error"),
    bef_diagnostic_skipped = c("bef_diagnostic_skipped", "bef_warning"),
    bef_sampler_diagnostics_failed = c("bef_sampler_diagnostics_failed", "bef_warning")
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

test_that("each condition class has its documented parents", {
  expected <- condition_expected_classes()
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
      extra_class = "bef_validate_error"
    ),
    c("bef_invalid_args", "bef_validate_error", "bef_pipeline_error", "bef_error")
  )

  err <- condition_capture_error(
    condition_ns(".bef_condition_classes")("bef_not_real")
  )
  expect_s3_class(err, "rlang_error")
  expect_match(conditionMessage(err), "Unknown bayesEfron condition class")
})

test_that("the error helpers raise their classes and keep extra fields", {
  abort_helpers <- list(
    bef_invalid_args = ".bef_abort_invalid_args",
    bef_compile_failed = ".bef_abort_compile_failed",
    bef_sampling_failed = ".bef_abort_sampling_failed",
    bef_sampling_partial = ".bef_abort_sampling_partial",
    bef_extraction_failed = ".bef_abort_extraction_failed",
    bef_invalid_fit = ".bef_abort_invalid_fit"
  )

  expected <- condition_expected_classes()
  for (class in names(abort_helpers)) {
    helper <- abort_helpers[[class]]
    err <- condition_capture_error(
      do.call(
        condition_ns(helper),
        list(message = "condition probe", condition_probe = helper)
      )
    )
    expected_prefix <- expected[[class]]
    if (identical(class, "bef_invalid_fit")) {
      expected_prefix <- c(
        "bef_invalid_fit", "bef_validate_error", "bef_pipeline_error", "bef_error"
      )
    }

    expect_condition_prefix(err, expected_prefix)
    expect_equal(err$condition_probe, helper)
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
    c("bef_invalid_args", "bef_validate_error", "bef_pipeline_error", "bef_error")
  )
  expect_equal(err$condition_probe, "validate")

  parent <- simpleError("parent probe")
  err <- condition_capture_error(
    condition_ns(".bef_abort_compile_failed")("compile probe", parent = parent)
  )
  expect_condition_prefix(err, expected$bef_compile_failed)
  expect_identical(err$parent, parent)
})

test_that("the warning helpers raise their classes", {
  warning_helpers <- list(
    bef_diagnostic_skipped = ".bef_warn_diagnostic_skipped",
    bef_sampler_diagnostics_failed = ".bef_warn_sampler_diagnostics_failed"
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

test_that("the grid errors keep their classes", {
  err <- condition_capture_error(
    make_efron_grid(c(-1, 0, 2), c(0.2, 0.3, 0.4), grid_method = "paper_simulation")
  )
  expect_condition_prefix(
    err,
    c("bef_err_grid_oracle_required", "bef_grid_error", "bef_error")
  )

  # A basis that has lost a column cannot be told apart from an intercept.
  grid <- make_efron_grid(c(-1, 0, 2, 3), c(0.2, 0.3, 0.4, 0.5), M = 4L)
  grid$B[, 4L] <- grid$B[, 1L]
  err <- condition_capture_error(
    condition_ns(".bef_validate_grid_return")(grid)
  )
  expect_condition_prefix(
    err,
    c("bef_grid_rank_deficient", "bef_grid_error", "bef_error")
  )
  expect_identical(err$expected_rank, 5L)
  expect_lt(err$rank, 5L)
})
