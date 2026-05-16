.bef_targets_schema <- function() {
  c(
    "target_id",
    "tier",
    "description",
    "tolerance_type",
    "tolerance_value",
    "expected_value",
    "expected_value_lower",
    "expected_value_upper",
    "fixture_path",
    "chapter_ref",
    "release_blocking",
    "status"
  )
}

.bef_runtime_env_schema <- function() {
  c(
    "target_id",
    "R_version_string",
    "splines_pkg_version",
    "RDS_serialization_version",
    "platform",
    "BLAS_LAPACK"
  )
}

.bef_testthat_file <- function(..., mustWork = TRUE) {
  relative_path <- file.path(...)
  path <- testthat::test_path(relative_path)
  if (!file.exists(path)) {
    required <- if (isTRUE(mustWork)) {
      file.path("tests", "testthat", relative_path)
    } else {
      c("DESCRIPTION", "NAMESPACE")
    }
    root <- bef_test_source_root(required = required)
    path <- file.path(root, "tests", "testthat", relative_path)
  }
  if (isTRUE(mustWork) && !file.exists(path)) {
    stop(sprintf("Testthat file does not exist: %s", path), call. = FALSE)
  }
  normalizePath(path, mustWork = mustWork)
}

.bef_source_file <- function(...,
                             required = c("DESCRIPTION", "NAMESPACE"),
                             mustWork = TRUE) {
  root <- bef_test_source_root(required = required)
  path <- file.path(root, ...)
  if (isTRUE(mustWork) && !file.exists(path)) {
    stop(sprintf("Source file does not exist: %s", path), call. = FALSE)
  }
  normalizePath(path, mustWork = mustWork)
}

.bef_read_csv <- function(path) {
  utils::read.csv(
    path,
    na.strings = c("", "NA"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.bef_parse_logical <- function(x, field) {
  if (is.logical(x)) {
    return(x)
  }
  if (is.character(x) && all(is.na(x) | x %in% c("TRUE", "FALSE"))) {
    return(ifelse(is.na(x), NA, x == "TRUE"))
  }
  stop(sprintf("`%s` must parse as logical.", field), call. = FALSE)
}

.bef_validate_targets <- function(x) {
  if (!identical(names(x), .bef_targets_schema())) {
    stop("verification_targets.csv schema does not match App-D section D.2.", call. = FALSE)
  }
  if (anyDuplicated(x$target_id)) {
    stop("verification_targets.csv contains duplicate target_id values.", call. = FALSE)
  }
  if (!all(x$tier %in% 0:3)) {
    stop("verification_targets.csv has invalid tier values.", call. = FALSE)
  }
  if (!all(x$tolerance_type %in% c("absolute", "relative", "interval", "BYTE"))) {
    stop("verification_targets.csv has invalid tolerance_type values.", call. = FALSE)
  }
  if (!all(x$status %in% c("active", "deferred", "obsolete"))) {
    stop("verification_targets.csv has invalid status values.", call. = FALSE)
  }

  x$release_blocking <- .bef_parse_logical(x$release_blocking, "release_blocking")

  scalar_tolerance <- x$tolerance_type %in% c("absolute", "relative")
  if (any(scalar_tolerance & is.na(x$tolerance_value))) {
    stop("absolute/relative targets require tolerance_value.", call. = FALSE)
  }
  if (any(x$tolerance_type %in% c("interval", "BYTE") & !is.na(x$tolerance_value))) {
    stop("interval/BYTE targets must leave tolerance_value as NA.", call. = FALSE)
  }
  if (any(x$tolerance_type == "interval" &
    (is.na(x$expected_value_lower) | is.na(x$expected_value_upper)))) {
    stop("interval targets require expected_value_lower and expected_value_upper.", call. = FALSE)
  }

  x
}

.bef_load_targets <- local({
  cache <- NULL

  function(path = .bef_testthat_file("verification_targets.csv"), refresh = FALSE) {
    if (!is.null(cache) && !isTRUE(refresh) && missing(path)) {
      return(cache)
    }
    out <- .bef_validate_targets(.bef_read_csv(path))
    if (missing(path)) {
      cache <<- out
    }
    out
  }
})

.bef_target <- function(target_id,
                        targets = .bef_load_targets(),
                        statuses = NULL,
                        require_one = TRUE) {
  if (!is.character(target_id) || length(target_id) != 1L || is.na(target_id)) {
    stop("`target_id` must be a single non-missing character value.", call. = FALSE)
  }
  hit <- targets[targets$target_id == target_id, , drop = FALSE]
  if (isTRUE(require_one) && nrow(hit) != 1L) {
    stop(sprintf("Expected exactly one verification target for `%s`.", target_id), call. = FALSE)
  }
  if (!is.null(statuses) && nrow(hit) > 0L && !hit$status %in% statuses) {
    stop(sprintf("Target `%s` has status `%s`.", target_id, hit$status), call. = FALSE)
  }
  if (!isTRUE(require_one)) {
    return(hit)
  }
  as.list(hit[1L, , drop = FALSE])
}

.bef_target_fixture_path <- function(target_id,
                                     targets = .bef_load_targets(),
                                     mustWork = TRUE) {
  target <- .bef_target(target_id, targets = targets)
  fixture_path <- target$fixture_path
  if (!is.character(fixture_path) || length(fixture_path) != 1L ||
    is.na(fixture_path) || !nzchar(fixture_path)) {
    stop(sprintf("Target `%s` has no fixture_path.", target_id), call. = FALSE)
  }
  .bef_testthat_file(fixture_path, mustWork = mustWork)
}

.bef_validate_runtime_env <- function(x, targets = .bef_load_targets()) {
  if (!identical(names(x), .bef_runtime_env_schema())) {
    stop("verification_runtime_env.csv schema does not match Ch08 section 8.2.B.", call. = FALSE)
  }
  if (anyDuplicated(x$target_id)) {
    stop("verification_runtime_env.csv contains duplicate target_id values.", call. = FALSE)
  }

  byte_ids <- targets$target_id[targets$tolerance_type == "BYTE"]
  if (!setequal(x$target_id, byte_ids)) {
    stop("verification_runtime_env.csv must contain exactly one row per BYTE target.", call. = FALSE)
  }
  x
}

.bef_load_runtime_env <- local({
  cache <- NULL

  function(path = .bef_testthat_file("verification_runtime_env.csv"),
           targets = .bef_load_targets(),
           refresh = FALSE) {
    if (!is.null(cache) && !isTRUE(refresh) && missing(path) && missing(targets)) {
      return(cache)
    }
    out <- .bef_validate_runtime_env(.bef_read_csv(path), targets = targets)
    if (missing(path) && missing(targets)) {
      cache <<- out
    }
    out
  }
})

.bef_byte_runtime_env <- function(target_id,
                                  targets = .bef_load_targets(),
                                  runtime_env = .bef_load_runtime_env(targets = targets),
                                  require_one = TRUE) {
  target <- .bef_target(target_id, targets = targets)
  if (!identical(target$tolerance_type, "BYTE")) {
    stop(sprintf("Target `%s` is not a BYTE target.", target_id), call. = FALSE)
  }
  hit <- runtime_env[runtime_env$target_id == target_id, , drop = FALSE]
  if (isTRUE(require_one) && nrow(hit) != 1L) {
    stop(sprintf("Expected exactly one runtime-env row for `%s`.", target_id), call. = FALSE)
  }
  if (!isTRUE(require_one)) {
    return(hit)
  }
  as.list(hit[1L, , drop = FALSE])
}

.bef_rds_serialization_version <- function() {
  if (getRversion() >= "3.5.0") {
    return("3")
  }
  "2"
}

.bef_current_runtime_env <- function(target_id = NA_character_) {
  versions <- extSoftVersion()
  data.frame(
    target_id = target_id,
    R_version_string = R.version.string,
    splines_pkg_version = as.character(utils::packageVersion("splines")),
    RDS_serialization_version = .bef_rds_serialization_version(),
    platform = R.version$platform,
    BLAS_LAPACK = paste(
      c("BLAS", "LAPACK"),
      as.character(versions[c("BLAS", "LAPACK")]),
      sep = "=",
      collapse = " | "
    ),
    row.names = NULL,
    check.names = FALSE
  )
}
