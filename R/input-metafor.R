#' Convert input data to a bayesEfron data object
#'
#' @description
#' Adapt the user's effect-size data into the standalone `bef_data`
#' class that [bayes_efron_fit()] consumes. Three input shapes are
#' supported, all returning a validated `bef_data` object that
#' carries `theta_hat`, `sigma`, optional site labels, and a `source`
#' attribute recording the input shape that produced it.
#'
#' @details
#' Supported input shapes:
#'
#' | Input class | Required fields | Source label |
#' |:------------|:----------------|:-------------|
#' | `bef_data` | (already converted; revalidated and returned) | unchanged |
#' | named `list` | `theta_hat`, `sigma`; optional `names` | `"list"` |
#' | `escalc` from [metafor::escalc()] | `yi`, `vi` (variance); optional row labels | `"metafor::escalc"` |
#'
#' For the `escalc` path, `theta_hat` is taken from `yi` and `sigma`
#' from `sqrt(vi)` after a strict-positivity check on `vi`. Row names
#' on the `escalc` object are propagated as site labels when present
#' and non-default.
#'
#' For the `list` path, an explicit `names` element overrides any
#' names attribute on `theta_hat`. Site labels are dropped when they
#' are missing, empty, or contain `NA`.
#'
#' Unsupported input classes raise a typed `bef_error` via
#' [as_bef_data.default()] with the offending class name in the
#' message.
#'
#' @param x Object to convert. v0.1 supports `bef_data`, `list` with
#'   `theta_hat` and `sigma`, and [metafor::escalc()] objects.
#' @param ... Reserved for future expansion; must be empty in v0.1.
#'
#' @return A validated `bef_data` object with a `source` attribute set
#'   to one of `"list"`, `"metafor::escalc"`, or unchanged for already
#'   converted inputs.
#'
#' @seealso
#'   * [bayes_efron_fit()], which calls `as_bef_data()` internally on
#'     its `theta_hat` / `sigma` arguments.
#'   * [metafor::escalc()] for computing effect sizes and sampling
#'     variances from study-level summary data.
#'
#' @examples
#' # Plain list input.
#' dat_list <- as_bef_data(list(
#'   theta_hat = c(-0.21, 0.04, 0.19, 0.38, 0.61),
#'   sigma     = c( 0.18, 0.15, 0.22, 0.19, 0.24)
#' ))
#' dat_list
#'
#' # Optional site labels (must satisfy the minimum-length-5 constraint).
#' dat_named <- as_bef_data(list(
#'   theta_hat = c(
#'     site_1 = -0.21, site_2 = 0.04, site_3 = 0.19,
#'     site_4 =  0.38, site_5 = 0.61
#'   ),
#'   sigma = c(0.18, 0.15, 0.22, 0.19, 0.24)
#' ))
#'
#' # metafor::escalc() bridge.
#' if (requireNamespace("metafor", quietly = TRUE)) {
#'   esc <- metafor::escalc(
#'     measure = "MD",
#'     m1i  = c(0.10, 0.40, 0.55, 0.30, 0.20),
#'     sd1i = c(0.30, 0.30, 0.35, 0.28, 0.32),
#'     n1i  = c( 60,   55,   62,   58,   65),
#'     m2i  = c(0.05, 0.10, 0.15, 0.08, 0.12),
#'     sd2i = c(0.32, 0.34, 0.36, 0.30, 0.33),
#'     n2i  = c( 60,   55,   62,   58,   65)
#'   )
#'   dat_esc <- as_bef_data(esc)
#'   dat_esc
#' }
#'
#' @export
as_bef_data <- function(x, ...) {
  if (inherits(x, "bef_data")) {
    .bef_check_as_bef_data_dots(list(...))
    return(validate_bef_data(x))
  }

  UseMethod("as_bef_data")
}

#' @rdname as_bef_data
#' @export
as_bef_data.default <- function(x, ...) {
  .bef_check_as_bef_data_dots(list(...))

  class_label <- paste(class(x), collapse = "/")
  if (!nzchar(class_label)) {
    class_label <- typeof(x)
  }

  .bef_abort_as_bef_data(
    sprintf(
      "`as_bef_data()` does not support objects of class <%s> in bayesEfron v0.1.",
      class_label
    ),
    arg = "x",
    predicate = "list, escalc, or bef_data"
  )
}

#' @rdname as_bef_data
#' @export
as_bef_data.list <- function(x, ...) {
  .bef_check_as_bef_data_dots(list(...))
  .bef_require_as_bef_data_fields(x, c("theta_hat", "sigma"), arg = "x")

  theta_hat <- x$theta_hat
  sigma <- x$sigma
  site_names <- .bef_list_site_names(x, theta_hat)

  validate_bef_data(
    new_bef_data(
      theta_hat = theta_hat,
      sigma = sigma,
      names = site_names,
      source = "list"
    )
  )
}

#' @rdname as_bef_data
#' @export
as_bef_data.escalc <- function(x, ...) {
  .bef_check_as_bef_data_dots(list(...))
  .bef_require_as_bef_data_fields(x, c("yi", "vi"), arg = "x")

  theta_hat <- x[["yi"]]
  vi <- x[["vi"]]
  .bef_check_escalc_variance(vi)

  validate_bef_data(
    new_bef_data(
      theta_hat = theta_hat,
      sigma = sqrt(vi),
      names = .bef_escalc_site_names(x),
      source = "metafor::escalc"
    )
  )
}

.bef_check_as_bef_data_dots <- function(dots) {
  if (length(dots) == 0L) {
    return(invisible(TRUE))
  }

  dot_names <- names(dots)
  dot_names <- dot_names[nzchar(dot_names)]
  unsupported <- if (length(dot_names) > 0L) {
    paste(sprintf("`%s`", dot_names), collapse = ", ")
  } else {
    "unnamed arguments"
  }

  .bef_abort_as_bef_data(
    sprintf("`...` is closed in bayesEfron v0.1; unsupported arguments: %s.", unsupported),
    arg = "...",
    predicate = "empty dots"
  )
}

.bef_require_as_bef_data_fields <- function(x, fields, arg) {
  missing <- setdiff(fields, names(x))
  if (length(missing) == 0L) {
    return(invisible(TRUE))
  }

  .bef_abort_as_bef_data(
    sprintf(
      "`%s` must contain fields: %s.",
      arg,
      paste(sprintf("`%s`", fields), collapse = ", ")
    ),
    arg = arg,
    predicate = paste(sprintf("field %s", fields), collapse = "; "),
    missing_fields = missing
  )
}

.bef_list_site_names <- function(x, theta_hat) {
  if ("names" %in% names(x)) {
    return(x$names)
  }

  candidate <- names(theta_hat)
  if (.bef_is_complete_label_vector(candidate, length(theta_hat))) {
    return(as.character(candidate))
  }

  NULL
}

.bef_escalc_site_names <- function(x) {
  candidate <- row.names(x)
  if (!.bef_is_complete_label_vector(candidate, NROW(x))) {
    return(NULL)
  }
  if (identical(candidate, as.character(seq_len(NROW(x))))) {
    return(NULL)
  }

  as.character(candidate)
}

.bef_check_escalc_variance <- function(vi) {
  if (!is.numeric(vi) ||
      length(vi) == 0L ||
      any(!is.finite(vi)) ||
      any(vi <= 0)) {
    .bef_abort_as_bef_data(
      "`x$vi` must be a strictly positive finite numeric vector.",
      arg = "x$vi",
      predicate = "strictly positive finite numeric vector"
    )
  }

  invisible(TRUE)
}

.bef_is_complete_label_vector <- function(x, n) {
  is.character(x) &&
    length(x) == n &&
    !anyNA(x) &&
    all(nzchar(x))
}

.bef_abort_as_bef_data <- function(message, ..., parent = NULL) {
  .bef_abort_invalid_args(
    message,
    ...,
    module = "as-bef-data",
    stage = 4L,
    parent = parent
  )
}
