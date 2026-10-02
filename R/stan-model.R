# Locating, compiling and reusing the Stan model.
#
# The compiled model is kept in a per-user cache directory and reused across
# sessions. Its file name carries a hash of the Stan file, the CmdStan version
# and the machine architecture, so a change in any of them gives a new file.

# Models already attached in this session, keyed by executable path.
.bef_models <- new.env(parent = emptyenv())

.bef_check_cmdstanr_installed <- function() {
  rlang::check_installed(
    "cmdstanr",
    reason = "to compile or run the bayesEfron Stan model"
  )
}

.bef_cmdstan_model <- function(...) {
  cmdstanr::cmdstan_model(...)
}

.bef_cmdstan_version <- function() {
  as.character(cmdstanr::cmdstan_version())
}

.bef_is_windows <- function() {
  identical(.Platform$OS.type, "windows")
}

.bef_stan_file <- function(model_family = "RE") {
  system.file("stan", "efron_re.stan", package = "bayesEfron", mustWork = TRUE)
}

.bef_stan_file_sha256 <- function(stan_file) {
  digest::digest(stan_file, algo = "sha256", file = TRUE)
}

.bef_cache_dir <- function(create = TRUE) {
  dir <- Sys.getenv("BAYESEFRON_CACHE_ROOT", unset = "")
  if (!nzchar(dir)) {
    dir <- tools::R_user_dir("bayesEfron", which = "cache")
  }
  if (isTRUE(create)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  dir
}

.bef_exe_file <- function(stan_file) {
  name <- sprintf(
    "efron_re-%s-cmdstan%s-%s",
    substr(.bef_stan_file_sha256(stan_file), 1L, 12L),
    .bef_cmdstan_version(),
    R.version$arch
  )
  if (.bef_is_windows()) {
    name <- paste0(name, ".exe")
  }
  file.path(.bef_cache_dir(), name)
}

# cmdstanr attaches to an executable without checking that it runs, so a
# damaged file would only fail at sampling time. On Unix the executable is
# asked for its build information; on Windows it needs CmdStan's libraries on
# the PATH to start, so only an empty file is detected there.
.bef_exe_ok <- function(exe) {
  if (!file.exists(exe) || file.size(exe) == 0) {
    return(FALSE)
  }
  if (.bef_is_windows()) {
    return(TRUE)
  }
  status <- tryCatch(
    suppressWarnings(system2(exe, "info", stdout = FALSE, stderr = FALSE)),
    error = function(err) 1L
  )
  identical(status, 0L)
}

.bef_model <- function(model_family = "RE",
                       force_recompile = FALSE,
                       quiet = TRUE) {
  .bef_check_cmdstanr_installed()
  stan_file <- .bef_stan_file(model_family)
  exe <- .bef_exe_file(stan_file)

  if (isTRUE(force_recompile)) {
    unlink(exe)
  } else if (!is.null(.bef_models[[exe]]) && file.exists(exe)) {
    return(.bef_models[[exe]])
  }

  model <- NULL
  if (file.exists(exe)) {
    # compile = FALSE: installing the package gives the Stan file a new
    # modification time, and cmdstanr would otherwise recompile after every
    # installation.
    model <- tryCatch(
      .bef_cmdstan_model(stan_file, exe_file = exe, compile = FALSE),
      error = function(err) NULL
    )
    if (is.null(model) || !.bef_exe_ok(exe)) {
      unlink(exe)
      model <- NULL
    }
  }
  if (is.null(model)) {
    model <- .bef_compile_model(stan_file, exe, quiet = quiet)
  }

  .bef_models[[exe]] <- model
  model
}

.bef_compile_model <- function(stan_file, exe, quiet = TRUE) {
  # Compile under a name that only this process uses, then rename, so that no
  # other session can see a partly written file under the final name. The
  # process id goes before the extension to keep ".exe" last on Windows.
  tmp <- sub("(\\.exe)?$", sprintf("-%d\\1", Sys.getpid()), exe)
  model <- tryCatch(
    .bef_cmdstan_model(
      stan_file,
      exe_file = tmp,
      stanc_options = list("O1"),
      quiet = quiet
    ),
    error = function(err) {
      .bef_abort_compile_failed(
        "CmdStan failed to compile the bayesEfron Stan model.",
        stan_file = stan_file,
        parent = err
      )
    }
  )
  if (!file.exists(tmp)) {
    .bef_abort_compile_failed(
      sprintf("The compiled model could not be written to `%s`.", dirname(exe)),
      cache_dir = dirname(exe)
    )
  }
  # If the rename fails, another session is running the existing file (this
  # happens on Windows); keep using the private copy in this session.
  if (isTRUE(suppressWarnings(file.rename(tmp, exe)))) {
    model$exe_file(exe)
  } else {
    exe <- tmp
  }
  if (!.bef_exe_ok(exe)) {
    unlink(exe)
    .bef_abort_compile_failed(
      "The compiled bayesEfron model does not run.",
      exe_file = exe
    )
  }
  model
}

#' Compile the Stan model ahead of time, or remove the compiled model
#'
#' @description
#' The first call to [bayes_efron_fit()] compiles the Stan model, which can
#' take up to a minute. `bayes_efron_compile()` does that step on its own,
#' so that a later fit starts sampling at once. `bayes_efron_clear_cache()`
#' deletes the compiled model.
#'
#' @details
#' The compiled model is stored in `tools::R_user_dir("bayesEfron", "cache")`
#' and reused in later sessions. The file name includes a hash of the Stan
#' file, the CmdStan version and the machine architecture, so the model is
#' compiled again after an update of the package's Stan file or of CmdStan.
#' Set the environment variable `BAYESEFRON_CACHE_ROOT` to store it somewhere
#' else, for example on a cluster where the home directory is not writable.
#'
#' `bayes_efron_compile()` needs the cmdstanr package and a working CmdStan
#' installation.
#'
#' @param model_family Character string. `"RE"` is the only model available.
#' @param quiet Logical. If `TRUE` (the default), the compiler output is not
#'   shown.
#' @param force_recompile Logical. If `TRUE`, the model is compiled even when
#'   a compiled copy exists. Defaults to `FALSE`.
#' @param scope Character string. `"all"` (the default) deletes the cached
#'   files. The values `"session"`, `"compiled_models"` and
#'   `"lock_only"` are accepted for scripts written for earlier versions and
#'   give a warning; `"lock_only"` no longer does anything.
#'
#' @return `bayes_efron_compile()` returns the `cmdstanr::CmdStanModel`
#'   object, invisibly. `bayes_efron_clear_cache()` returns the number of
#'   files it deleted, invisibly.
#'
#' @seealso [bayes_efron_fit()]
#'
#' @examples
#' \dontrun{
#' bayes_efron_compile()
#'
#' # Start again from the Stan file.
#' bayes_efron_clear_cache()
#' bayes_efron_compile()
#' }
#'
#' @export
bayes_efron_compile <- function(model_family = "RE",
                                quiet = TRUE,
                                force_recompile = FALSE) {
  .bef_check_compile_arg(
    checkmate::assert_choice(model_family, choices = "RE"),
    arg = "model_family"
  )
  .bef_check_compile_arg(checkmate::assert_flag(quiet), arg = "quiet")
  .bef_check_compile_arg(
    checkmate::assert_flag(force_recompile),
    arg = "force_recompile"
  )
  invisible(.bef_model(
    model_family = model_family,
    force_recompile = force_recompile,
    quiet = quiet
  ))
}

#' @rdname bayes_efron_compile
#' @export
bayes_efron_clear_cache <- function(scope = "all") {
  old_scopes <- c("session", "compiled_models", "lock_only")
  .bef_check_compile_arg(
    checkmate::assert_choice(scope, choices = c("all", old_scopes)),
    arg = "scope"
  )
  if (scope %in% old_scopes) {
    warning(
      sprintf(
        "`scope = \"%s\"` is deprecated; call `bayes_efron_clear_cache()` without arguments.",
        scope
      ),
      call. = FALSE
    )
  }

  removed <- 0L
  dir <- .bef_cache_dir(create = FALSE)
  if (scope %in% c("all", "compiled_models") && dir.exists(dir)) {
    # The "v1" subdirectory was used by earlier versions of the package.
    files <- c(
      list.files(dir, pattern = "^efron_re-", full.names = TRUE),
      list.files(file.path(dir, "v1"), full.names = TRUE, all.files = TRUE,
                 no.. = TRUE)
    )
    unlink(file.path(dir, "v1"), recursive = TRUE)
    unlink(files)
    removed <- sum(!file.exists(files))
  }
  if (!identical(scope, "lock_only")) {
    rm(list = ls(.bef_models, all.names = TRUE), envir = .bef_models)
  }

  message(sprintf(
    "bayesEfron: removed %d cached file%s from %s.",
    removed, if (removed == 1L) "" else "s", dir
  ))
  invisible(removed)
}

.bef_check_compile_arg <- function(expr, arg) {
  tryCatch(
    {
      force(expr)
      invisible(TRUE)
    },
    error = function(err) {
      .bef_abort_invalid_args(
        sprintf("`%s` is not valid: %s", arg, conditionMessage(err)),
        arg = arg,
        parent = err
      )
    }
  )
}
