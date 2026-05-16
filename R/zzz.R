.bayesEfron_cache <- NULL
.bayesEfron_msgs_emitted <- NULL
.bayesEfron_session <- NULL

.onLoad <- function(libname, pkgname) {
  .bef_initialize_session_state(pkgname = pkgname)
}

.bef_initialize_session_state <- function(pkgname = "bayesEfron") {
  .bayesEfron_cache <<- new.env(hash = TRUE, parent = emptyenv())
  .bayesEfron_msgs_emitted <<- new.env(hash = TRUE, parent = emptyenv())
  .bayesEfron_session <<- new.env(hash = TRUE, parent = emptyenv())
  .bef_stamp_session_state(pkgname = pkgname)
  invisible(TRUE)
}

.bef_ensure_session_state <- function(pkgname = "bayesEfron") {
  if (!is.environment(.bayesEfron_cache)) {
    .bayesEfron_cache <<- new.env(hash = TRUE, parent = emptyenv())
  }
  if (!is.environment(.bayesEfron_msgs_emitted)) {
    .bayesEfron_msgs_emitted <<- new.env(hash = TRUE, parent = emptyenv())
  }
  if (!is.environment(.bayesEfron_session)) {
    .bayesEfron_session <<- new.env(hash = TRUE, parent = emptyenv())
  }
  if (is.null(.bayesEfron_session$session_id)) {
    .bef_stamp_session_state(pkgname = pkgname)
  }
  invisible(TRUE)
}

.bef_stamp_session_state <- function(pkgname = "bayesEfron") {
  .bayesEfron_session$session_id <- paste(
    Sys.getpid(),
    format(Sys.time(), "%Y%m%d%H%M%OS6"),
    sep = "-"
  )
  .bayesEfron_session$cache_format_version <- .bef_cache_format_version()
  .bayesEfron_session$cache_root <- .bef_cache_root()
  .bayesEfron_session$cache_dir <- .bef_cache_dir(create = FALSE)
  .bayesEfron_session$package_version <- .bef_package_version(pkgname)
  invisible(TRUE)
}

.bef_package_version <- function(pkgname = "bayesEfron") {
  tryCatch(
    as.character(utils::packageVersion(pkgname)),
    error = function(err) NA_character_
  )
}
