#' Extract estimates and standard errors for a fit
#'
#' @description
#' [bayes_efron_fit()] takes a vector of estimates and a vector of their
#' standard errors. `as_bef_data()` takes the two vectors out of the object
#' that usually holds them, a list, a data frame or the result of
#' `metafor::escalc()`, and checks them. The result is a `bef_data` object,
#' whose elements `theta_hat` and `sigma` are the arguments of
#' [bayes_efron_fit()].
#'
#' @details
#' The methods differ in where they look for the two vectors.
#'
#' For a list, the elements `theta_hat` and `sigma` are used. Site labels
#' are taken from an element `names` or, if there is none, from the names
#' of `theta_hat`.
#'
#' For a data frame, the estimates are the column named `theta_hat`, `yi`,
#' `estimate`, `effect` or `y`, and the standard errors are the column
#' named `sigma`, `se`, `std_error`, `stderr` or `sei`; in each case the
#' first of these names that the data frame has. The arguments `theta_hat`
#' and `sigma` name other columns. If the column holds variances, as the
#' column `vi` of a metafor data set does, name it in `sigma` and set
#' `variance = TRUE`. Site labels are taken from the column named in `site`
#' or, without it, from the row names if they have been set.
#'
#' For an object of class `escalc`, the estimates are the column `yi` and
#' the standard errors are the square roots of the column `vi`. Row names
#' that have been set are used as site labels.
#'
#' A `bef_data` object is checked again and returned.
#'
#' At least five sites are needed, and the function stops if a value is
#' missing or a standard error is not positive. The site labels are shown
#' when the object is printed; a fitted model numbers the sites in the
#' order in which they were given.
#'
#' @param x A list, a data frame, an `escalc` object from the metafor
#'   package or a `bef_data` object.
#' @param ... For `as_bef_data()`, these dots are for future extensions and
#'   must be empty. The `print()` method passes them on to `format()`;
#'   `format()` and `summary()` do not use them.
#' @param theta_hat,sigma For a data frame, character strings naming the
#'   columns of the estimates and of the standard errors. `NULL` (the
#'   default) means that the columns are found by name, as described in
#'   Details.
#' @param variance Logical. For a data frame, `TRUE` if the column named in
#'   `sigma` holds variances; their square roots are then taken. Defaults
#'   to `FALSE`.
#' @param site For a data frame, a character string naming a column of
#'   site labels. Defaults to `NULL`.
#' @param object A `bef_data` object.
#' @inheritParams bef_fit
#'
#' @return A list of class `bef_data` with the elements `theta_hat`, the
#'   estimates; `sigma`, the standard errors; `names`, the site labels or
#'   `NULL`; and `source`, which is `"list"`, `"data.frame"` or
#'   `"metafor::escalc"` according to the class of `x`.
#'
#' @seealso [bayes_efron_fit()], [raudenbush1985]
#'
#' @examples
#' # A data frame: the columns `yi` and `sei` are found by name.
#' dat <- as_bef_data(raudenbush1985)
#' dat
#'
#' # The same from the variances, with the authors as site labels.
#' as_bef_data(raudenbush1985, sigma = "vi", variance = TRUE, site = "author")
#'
#' # A list.
#' as_bef_data(list(theta_hat = raudenbush1985$yi, sigma = raudenbush1985$sei))
#'
#' # The two vectors that bayes_efron_fit() takes.
#' str(dat[c("theta_hat", "sigma")])
#'
#' @examplesIf requireNamespace("metafor", quietly = TRUE)
#' # An escalc object from the metafor package.
#' es <- metafor::escalc(measure = "GEN", yi = yi, vi = vi, data = raudenbush1985)
#' as_bef_data(es)
#'
#' @order 1
#' @export
as_bef_data <- function(x, ...) {
  if (inherits(x, "bef_data")) {
    .bef_check_as_bef_data_dots(list(...))
    return(validate_bef_data(x))
  }

  UseMethod("as_bef_data")
}

#' @order 2
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
      "`as_bef_data()` has no method for an object of class <%s>.",
      class_label
    ),
    arg = "x"
  )
}

#' @order 3
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

#' @order 4
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
    sprintf("Unused arguments: %s.", unsupported),
    arg = "..."
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
      arg = "x$vi"
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
    parent = parent
  )
}

#' @order 5
#' @rdname as_bef_data
#' @export
as_bef_data.data.frame <- function(x,
                                   ...,
                                   theta_hat = NULL,
                                   sigma = NULL,
                                   variance = FALSE,
                                   site = NULL) {
  .bef_check_as_bef_data_dots(list(...))

  theta_hat <- .bef_match_df_column(
    x, theta_hat, c("theta_hat", "yi", "estimate", "effect", "y"), "theta_hat"
  )
  sigma <- .bef_match_df_column(
    x, sigma, c("sigma", "se", "std_error", "stderr", "sei"), "sigma"
  )
  if (!is.logical(variance) || length(variance) != 1L || is.na(variance)) {
    .bef_abort_as_bef_data(
      "`variance` must be a single non-missing logical.",
      arg = "variance"
    )
  }

  values <- x[[theta_hat]]
  errors <- x[[sigma]]
  if (!is.numeric(values) || !is.numeric(errors)) {
    .bef_abort_as_bef_data(
      sprintf("Columns `%s` and `%s` must both be numeric.", theta_hat, sigma),
      arg = "x"
    )
  }
  # Missing values are an error. Dropping them would change the number of
  # sites without the caller noticing.
  if (anyNA(values) || anyNA(errors)) {
    .bef_abort_as_bef_data(
      sprintf(
        "Columns `%s` and `%s` must not contain missing values (found %d and %d).",
        theta_hat, sigma, sum(is.na(values)), sum(is.na(errors))
      ),
      arg = "x"
    )
  }
  if (isTRUE(variance)) {
    if (any(errors < 0)) {
      .bef_abort_as_bef_data(
        sprintf("Variance column `%s` must be non-negative.", sigma),
        arg = "x"
      )
    }
    errors <- sqrt(errors)
  }

  labels <- NULL
  if (!is.null(site)) {
    site <- .bef_match_df_column(x, site, character(), "site")
    labels <- as.character(x[[site]])
    if (!.bef_is_complete_label_vector(labels, nrow(x))) {
      .bef_abort_as_bef_data(
        sprintf("Site column `%s` must have one non-missing label per row.", site),
        arg = "x"
      )
    }
  } else if (!is.null(rownames(x)) && !identical(rownames(x), as.character(seq_len(nrow(x))))) {
    labels <- rownames(x)
  }

  validate_bef_data(new_bef_data(
    theta_hat = as.numeric(values),
    sigma = as.numeric(errors),
    names = labels,
    source = "data.frame"
  ))
}

# The column to use: the one named by the caller, or the first conventional
# name that the data frame has.
.bef_match_df_column <- function(x, given, conventional, arg) {
  if (!is.null(given)) {
    if (!.bef_is_string(given) || !given %in% names(x)) {
      .bef_abort_as_bef_data(
        sprintf("`%s` must name a column of `x`; got \"%s\".", arg,
                paste(given, collapse = ", ")),
        arg = arg
      )
    }
    return(given)
  }
  found <- intersect(conventional, names(x))
  if (length(found) == 0L) {
    .bef_abort_as_bef_data(
      sprintf(
        "Cannot infer the `%s` column. Name it explicitly, or use one of: %s.",
        arg, paste(conventional, collapse = ", ")
      ),
      arg = arg
    )
  }
  found[[1L]]
}
