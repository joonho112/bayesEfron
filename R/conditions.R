# The frame to name in an error or a warning: the outermost call into this
# package, so that the message names the function that the user called and
# not an internal helper.
.bef_user_frame <- function() {
  namespace <- topenv(environment(.bef_user_frame))
  for (i in seq_len(sys.nframe())) {
    fn <- sys.function(i)
    if (is.function(fn) && identical(topenv(environment(fn)), namespace)) {
      return(sys.frame(i))
    }
  }
  NULL
}

.bef_abort <- function(message,
                       class,
                       ...,
                       parent = NULL,
                       extra_class = NULL,
                       call = .bef_user_frame()) {
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
                      call = .bef_user_frame()) {
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
    bef_extraction_failed = c("bef_extraction_failed", "bef_pipeline_error", "bef_error"),
    bef_invalid_fit = c("bef_invalid_fit", "bef_pipeline_error", "bef_error"),
    bef_diagnostic_skipped = c("bef_diagnostic_skipped", "bef_warning"),
    bef_sampler_diagnostics_failed = c("bef_sampler_diagnostics_failed", "bef_warning"),
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

.bef_abort_invalid_args <- function(message, ..., parent = NULL, validate = FALSE) {
  extra_class <- if (isTRUE(validate)) "bef_validate_error" else NULL
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

.bef_abort_extraction_failed <- function(message, ..., parent = NULL) {
  .bef_abort(message, "bef_extraction_failed", ..., parent = parent)
}

.bef_abort_invalid_fit <- function(message, ..., parent = NULL) {
  .bef_abort(
    message,
    "bef_invalid_fit",
    ...,
    parent = parent,
    extra_class = "bef_validate_error"
  )
}

.bef_warn_diagnostic_skipped <- function(message, ..., parent = NULL) {
  .bef_warn(message, "bef_diagnostic_skipped", ..., parent = parent)
}

.bef_warn_sampler_diagnostics_failed <- function(message, ..., parent = NULL) {
  .bef_warn(message, "bef_sampler_diagnostics_failed", ..., parent = parent)
}
