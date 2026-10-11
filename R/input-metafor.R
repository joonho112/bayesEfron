#' Extract estimates and standard errors from a data frame or list
#'
#' @description
#' [bayes_efron_fit()] needs a vector of effect estimates and a vector of
#' their standard errors. `as_bef_data()` extracts the two vectors from a
#' data frame, a list or an `escalc` object from the metafor package and
#' checks them. The result is a `bef_data` object, which can be passed to
#' [bayes_efron_fit()] in place of the two vectors.
#'
#' @details
#' The methods differ in where they look for the estimates and the
#' standard errors.
#'
#' For a list, the elements `theta_hat` and `sigma` are used. Site labels
#' are taken from an element `names` or, if there is none, from the names
#' of `theta_hat`.
#'
#' For a data frame, the estimates are the column named `theta_hat`, `yi`,
#' `estimate`, `effect` or `y`, and the standard errors are the column
#' named `sigma`, `se`, `std_error`, `stderr` or `sei`. If the data frame
#' has more than one of these columns, the first in the list is used. If
#' the data frame has none
#' of the columns for standard errors but a column `vi`, `var` or
#' `variance`, as a metafor data set does, the square roots of that column
#' are used, with a message. Use the arguments `theta_hat` and `sigma` to
#' specify columns with other names. For a column of variances with
#' another name, specify it as `sigma` and set `variance = TRUE`. If
#' `sigma` names a column whose name begins with `vi` or `var` and
#' `variance` is not `TRUE`, a warning is issued, because the name
#' suggests that the column contains variances. Site labels are taken
#' from the column named in `site`. Without `site`, the row names are
#' used if they have been set.
#'
#' For an `escalc` object, the estimates are the column `yi` and the
#' standard errors are the square roots of the column `vi`. Site labels
#' are the study labels given to `metafor::escalc()` in its argument
#' `slab` or, without them, row names that have been set.
#'
#' A `bef_data` object is checked again and returned.
#'
#' At least 5 sites are required. An error is thrown if a value is missing,
#' if a standard error is not positive, or if a site label occurs more
#' than once. The site labels are stored in a fitted model: they are shown
#' in its tables and plots, and sites can be selected by their labels.
#'
#' @param x For `as_bef_data()`, a data frame, a list, an `escalc` object
#'   from the metafor package or a `bef_data` object; see Details. For
#'   `print()`, `format()` and `as.data.frame()`, a `bef_data` object.
#' @param ... For `as_bef_data()`, these dots are for future extensions and
#'   must be empty. The `print()` method passes them on to `format()`; the
#'   other methods ignore them.
#' @param theta_hat,sigma For a data frame, the names of the columns that
#'   contain the estimates and the standard errors, as character strings.
#'   If `NULL` (the default), the columns are found by name; see Details.
#' @param variance For a data frame, set to `TRUE` if the column specified
#'   in `sigma` contains variances. Their square roots are then used as
#'   the standard errors.
#' @param site For a data frame, an optional character string naming a
#'   column of site labels. If `NULL` (the default), row names that have
#'   been set are used as labels.
#' @param object A `bef_data` object.
#' @param use_cli For `format()`: if `TRUE` and the cli package is
#'   installed, the first line is in bold; if `FALSE`, the result is plain
#'   text. If `NULL` (the default), the value of
#'   `getOption("bayesEfron.use_cli", TRUE)` is used.
#'
#' @return An object of class `bef_data`, which is a list with the
#'   following components.
#'   \describe{
#'     \item{`theta_hat`}{The estimates.}
#'     \item{`sigma`}{The standard errors.}
#'     \item{`names`}{The site labels, or `NULL` if there are none.}
#'     \item{`source`}{`"list"`, `"data.frame"` or `"metafor::escalc"`,
#'       depending on the class of `x`.}
#'   }
#'
#' @seealso [bayes_efron_fit()], [krw_discrimination]
#'
#' @examples
#' # A data frame: the columns `yi` and `sei` are found by name.
#' dat <- as_bef_data(krw_discrimination)
#' dat
#'
#' # The same from the variances, with firm IDs as site labels.
#' as_bef_data(krw_discrimination, sigma = "vi", variance = TRUE, site = "firm_id")
#'
#' # A list.
#' as_bef_data(list(theta_hat = krw_discrimination$yi, sigma = krw_discrimination$sei))
#'
#' # The two vectors that bayes_efron_fit() takes.
#' str(dat[c("theta_hat", "sigma")])
#'
#' @examplesIf requireNamespace("metafor", quietly = TRUE)
#' # An escalc object from the metafor package.
#' es <- metafor::escalc(measure = "GEN", yi = yi, vi = vi, data = krw_discrimination)
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
  # metafor keeps study labels in the attribute `slab` of the estimates.
  slab <- attr(x[["yi"]], "slab", exact = TRUE)
  if (!is.null(slab) && length(slab) == NROW(x) && !anyNA(slab)) {
    slab <- as.character(slab)
    if (all(nzchar(slab)) && !identical(slab, as.character(seq_len(NROW(x))))) {
      return(slab)
    }
  }
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
  if (!is.logical(variance) || length(variance) != 1L || is.na(variance)) {
    .bef_abort_as_bef_data(
      "`variance` must be a single non-missing logical.",
      arg = "variance"
    )
  }
  standard_errors <- c("sigma", "se", "std_error", "stderr", "sei")
  variances <- c("vi", "var", "variance")
  if (is.null(sigma) && !any(standard_errors %in% names(x)) && any(variances %in% names(x))) {
    # A data frame with estimates and sampling variances only, as
    # metafor::escalc() returns it.
    sigma <- intersect(variances, names(x))[[1L]]
    if (!isTRUE(variance)) {
      message(sprintf(
        "No column of standard errors was found; using the square roots of the variances in `%s`.",
        sigma
      ))
      variance <- TRUE
    }
  } else {
    sigma <- .bef_match_df_column(x, sigma, standard_errors, "sigma")
    if (!isTRUE(variance) && grepl("^(vi|var)", sigma, ignore.case = TRUE)) {
      warning(
        sprintf(
          paste0(
            "Column `%s` is used as the standard errors, but its name suggests ",
            "variances. If it holds variances, set `variance = TRUE`."
          ),
          sigma
        ),
        call. = FALSE
      )
    }
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
