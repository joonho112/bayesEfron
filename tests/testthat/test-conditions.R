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
  grid <- make_efron_grid(c(-1, 0, 2, 3), c(0.2, 0.3, 0.4, 0.5), M = 4L, grid_method = "paper_realdata")
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

test_that("missing values, infinite values and the length of the input are reported separately", {
  problem <- getFromNamespace(".bef_vector_problem", "bayesEfron")
  expect_null(problem(c(0.1, 0.2, 0.3, 0.4, 0.5), "theta_hat", min_len = 5L))
  expect_identical(problem(c(1, NA, 3, NA, 5), "theta_hat", min_len = 5L),
                   "`theta_hat` has 2 missing values (the first at position 2).")
  expect_identical(problem(c(1, Inf, 3, 4, 5), "theta_hat", min_len = 5L),
                   "`theta_hat` has 1 infinite value (the first at position 2).")
  expect_identical(problem(1:4, "theta_hat", min_len = 5L),
                   "`theta_hat` must have at least 5 values, one for each site; it has 4.")
  expect_identical(problem(rep(0.2, 3L), "sigma", len = 5L, positive = TRUE),
                   "`sigma` must have one value for each of the 5 estimates; it has 3.")
  expect_identical(problem(c(0.2, 0, -1, 0.2, 0.2), "sigma", len = 5L, positive = TRUE),
                   "`sigma` must be positive; it has 2 values that are not (the first at position 2).")
  expect_identical(problem("a", "theta_hat"), "`theta_hat` must be a numeric vector, not character.")
  expect_identical(problem(NULL, "sigma"), "`sigma` is missing; give one value for each site.")

  # The message of another argument is one sentence without a parent error.
  err <- tryCatch(
    getFromNamespace(".bef_validate_fit_args", "bayesEfron")(1:5 / 10, rep(0.2, 5L), L = 20),
    error = identity
  )
  expect_s3_class(err, "bef_invalid_args")
  expect_null(err$parent)
  expect_identical(conditionMessage(err),
                   "`L` must be an integer scalar in [51, 401], or `NULL`: element 1 is not >= 51.")
})

test_that("every condition class that the package raises is registered", {
  # Walk the body of every function of the package and collect the classes
  # given as text to .bef_warn() and .bef_abort(). A class that is raised
  # and not registered makes the condition itself fail with an internal
  # error, which is how a failed fit stopped compare_efron_grids() in
  # version 0.3.0.
  ns <- asNamespace("bayesEfron")
  registered <- getFromNamespace(".bef_condition_classes", "bayesEfron")
  found <- character()
  walk <- function(expr) {
    if (is.call(expr)) {
      head <- expr[[1L]]
      if (is.name(head) && as.character(head) %in% c(".bef_warn", ".bef_abort")) {
        args <- as.list(expr)[-1L]
        class_arg <- if ("class" %in% names(args)) args[["class"]] else if (length(args) >= 2L) args[[2L]]
        if (is.character(class_arg)) found <<- c(found, class_arg)
      }
      # An empty argument, as in x[, 1], is skipped.
      for (i in seq_along(expr)) {
        if (!identical(expr[[i]], quote(expr = ))) walk(expr[[i]])
      }
    }
  }
  for (name in ls(ns, all.names = TRUE)) {
    object <- get(name, envir = ns)
    if (is.function(object) && !is.primitive(object)) walk(body(object))
  }
  found <- unique(found)
  expect_true("bef_grid_comparison_fit_failed" %in% found)
  expect_gte(length(found), 5L)
  for (class in found) {
    expect_no_error(registered(class))
  }
  expect_identical(registered("bef_grid_comparison_fit_failed"),
                   c("bef_grid_comparison_fit_failed", "bef_warning"))
  expect_error(registered("bef_not_a_class"), "Unknown bayesEfron condition class")
})
