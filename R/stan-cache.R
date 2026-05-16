.bef_cache_format_version <- function() {
  "v1"
}

.bef_cache_root <- function() {
  root <- Sys.getenv("BAYESEFRON_CACHE_ROOT", unset = "")
  if (!nzchar(root)) {
    root <- tools::R_user_dir("bayesEfron", which = "cache")
  }

  normalizePath(path.expand(root), winslash = "/", mustWork = FALSE)
}

.bef_cache_dir <- function(create = FALSE) {
  path <- file.path(.bef_cache_root(), .bef_cache_format_version())
  if (isTRUE(create)) {
    path <- .bef_ensure_cache_dir(path)
  }
  path
}

.bef_ensure_cache_dir <- function(path = .bef_cache_dir(create = FALSE)) {
  if (!dir.exists(path)) {
    ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE, mode = "0700")
    if (!isTRUE(ok) && !dir.exists(path)) {
      .bef_abort_cache_perm_violation(
        sprintf("Failed to create bayesEfron cache directory: %s", path),
        cache_dir = path
      )
    }
  }

  if (!identical(.Platform$OS.type, "windows")) {
    suppressWarnings(Sys.chmod(path, mode = "0700"))
  }

  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  .bef_validate_cache_dir_owner(path)
  path
}

.bef_cache_lock_path <- function() {
  file.path(.bef_cache_dir(create = FALSE), ".lock")
}

.bef_cache_lock_timeout <- function() {
  600L
}

.bef_acquire_lock <- function(timeout = .bef_cache_lock_timeout(),
                              poll_interval = 0.5,
                              lock_path = .bef_cache_lock_path(),
                              now = Sys.time,
                              sleep = Sys.sleep) {
  timeout <- .bef_lock_number(timeout, "timeout", lower = 0)
  poll_interval <- .bef_lock_number(
    poll_interval, "poll_interval", lower = 0
  )
  .bef_ensure_cache_dir(dirname(lock_path))

  started_at <- now()
  owner <- NULL
  repeat {
    payload <- .bef_lock_payload(acquired_at = now())
    if (isTRUE(.bef_try_create_lock(lock_path, payload))) {
      return(invisible(lock_path))
    }

    owner <- .bef_read_lock_payload(lock_path)
    if (.bef_lock_is_stale(owner)) {
      unlink(lock_path, force = TRUE)
      next
    }

    elapsed <- as.numeric(difftime(now(), started_at, units = "secs"))
    if (is.finite(elapsed) && elapsed >= timeout) {
      .bef_abort_lock_timeout(
        sprintf(
          "Timed out after %.3f seconds waiting for bayesEfron cache lock.",
          timeout
        ),
        lock_path = lock_path,
        timeout = timeout,
        poll_interval = poll_interval,
        started_at = .bef_format_lock_time(started_at),
        elapsed = elapsed,
        owner_pid = .bef_lock_field(owner, "pid"),
        owner_host = .bef_lock_field(owner, "nodename"),
        owner_started_at = .bef_lock_field(owner, "process_start_time"),
        owner_acquired_at = .bef_lock_field(owner, "acquired_at"),
        parent = .bef_lock_owner_condition(lock_path, owner)
      )
    }

    sleep(poll_interval)
  }
}

.bef_release_lock <- function(lock_path = .bef_cache_lock_path()) {
  if (file.exists(lock_path) && .bef_lock_is_owned_by_current_process(lock_path)) {
    unlink(lock_path, force = TRUE)
  }
  invisible(TRUE)
}

.bef_try_create_lock <- function(lock_path, payload) {
  con <- tryCatch(
    suppressWarnings(file(lock_path, open = "wx")),
    error = function(err) NULL
  )
  if (is.null(con)) {
    return(FALSE)
  }
  close_con <- TRUE
  on.exit({
    if (isTRUE(close_con)) {
      try(close(con), silent = TRUE)
    }
  }, add = TRUE)

  tryCatch(
    {
      writeLines(.bef_lock_payload_lines(payload), con, useBytes = TRUE)
      close(con)
      close_con <- FALSE
      TRUE
    },
    error = function(err) {
      try(close(con), silent = TRUE)
      close_con <<- FALSE
      if (.bef_lock_is_owned_by_current_process(lock_path)) {
        unlink(lock_path, force = TRUE)
      }
      .bef_abort_cache_perm_violation(
        sprintf("Failed to write bayesEfron cache lock: %s", lock_path),
        lock_path = lock_path,
        parent = err
      )
    }
  )
}

.bef_write_lock_payload <- function(lock_path, payload = .bef_lock_payload()) {
  writeLines(.bef_lock_payload_lines(payload), lock_path, useBytes = TRUE)
  invisible(payload)
}

.bef_lock_payload <- function(pid = Sys.getpid(),
                              nodename = Sys.info()[["nodename"]],
                              start_time = .bef_process_start_time(pid),
                              acquired_at = Sys.time()) {
  list(
    pid = as.character(pid),
    nodename = .bef_lock_string(nodename, "unknown-host"),
    process_start_time = .bef_lock_string(start_time, "unknown-start-time"),
    acquired_at = .bef_format_lock_time(acquired_at)
  )
}

.bef_lock_payload_lines <- function(payload) {
  c(
    payload$pid,
    payload$nodename,
    payload$process_start_time,
    payload$acquired_at
  )
}

.bef_read_lock_payload <- function(lock_path = .bef_cache_lock_path()) {
  if (!file.exists(lock_path)) {
    return(.bef_malformed_lock_payload("missing"))
  }

  lines <- tryCatch(
    readLines(lock_path, warn = FALSE),
    error = function(err) character()
  )
  if (length(lines) != 4L || any(!nzchar(lines))) {
    return(.bef_malformed_lock_payload("malformed"))
  }

  list(
    pid = lines[[1L]],
    nodename = lines[[2L]],
    process_start_time = lines[[3L]],
    acquired_at = lines[[4L]],
    malformed = FALSE
  )
}

.bef_malformed_lock_payload <- function(reason) {
  list(
    pid = NA_character_,
    nodename = NA_character_,
    process_start_time = NA_character_,
    acquired_at = NA_character_,
    malformed = TRUE,
    malformed_reason = reason
  )
}

.bef_lock_is_stale <- function(payload,
                               current_nodename = Sys.info()[["nodename"]]) {
  if (!is.list(payload) || isTRUE(payload$malformed)) {
    return(FALSE)
  }
  if (!identical(payload$nodename, .bef_lock_string(current_nodename, "unknown-host"))) {
    return(FALSE)
  }

  pid <- suppressWarnings(as.integer(payload$pid))
  if (is.na(pid) || pid < 1L) {
    return(FALSE)
  }
  if (!.bef_process_alive(pid)) {
    return(TRUE)
  }

  current_start_time <- .bef_process_start_time(pid)
  if (!.bef_lock_start_time_known(current_start_time) ||
      !.bef_lock_start_time_known(payload$process_start_time)) {
    return(FALSE)
  }

  !identical(current_start_time, payload$process_start_time)
}

.bef_lock_is_owned_by_current_process <- function(lock_path = .bef_cache_lock_path()) {
  payload <- .bef_read_lock_payload(lock_path)
  if (!is.list(payload) || isTRUE(payload$malformed)) {
    return(FALSE)
  }
  identical(payload$pid, as.character(Sys.getpid())) &&
    identical(payload$nodename, .bef_lock_string(Sys.info()[["nodename"]], "unknown-host")) &&
    identical(payload$process_start_time, .bef_process_start_time(Sys.getpid()))
}

.bef_process_alive <- function(pid) {
  pid <- suppressWarnings(as.integer(pid))
  if (is.na(pid) || pid < 1L) {
    return(FALSE)
  }

  if (identical(.Platform$OS.type, "unix")) {
    proc_path <- file.path("/proc", as.character(pid))
    if (dir.exists("/proc")) {
      return(file.exists(proc_path))
    }
    status <- suppressWarnings(system2(
      "ps", c("-p", as.character(pid)), stdout = FALSE, stderr = FALSE
    ))
    return(identical(status, 0L))
  }

  .bef_windows_process_alive(pid)
}

.bef_windows_process_alive <- function(pid) {
  if (requireNamespace("ps", quietly = TRUE)) {
    alive <- tryCatch(
      ps::ps_is_running(ps::ps_handle(pid)),
      error = function(err) FALSE
    )
    if (isTRUE(alive)) {
      return(TRUE)
    }
  }

  .bef_windows_tasklist_process_alive(pid)
}

.bef_windows_tasklist_process_alive <- function(pid, system2_fun = system2) {
  pid <- suppressWarnings(as.integer(pid))
  if (is.na(pid) || pid < 1L) {
    return(FALSE)
  }

  pid_chr <- as.character(pid)
  output <- tryCatch(
    system2_fun(
      "tasklist",
      c("/FI", paste("PID eq", pid_chr), "/FO", "CSV", "/NH"),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(err) character()
  )
  status <- attr(output, "status", exact = TRUE)
  if (!is.null(status) && !identical(status, 0L)) {
    return(FALSE)
  }

  lines <- trimws(output)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0L || any(grepl("^INFO:", lines, ignore.case = TRUE))) {
    return(FALSE)
  }

  any(vapply(lines, .bef_tasklist_line_has_pid, logical(1), pid = pid_chr))
}

.bef_tasklist_line_has_pid <- function(line, pid) {
  row <- tryCatch(
    utils::read.csv(
      text = line,
      header = FALSE,
      stringsAsFactors = FALSE,
      colClasses = "character"
    ),
    error = function(err) NULL
  )
  if (is.null(row) || nrow(row) < 1L || ncol(row) < 2L) {
    return(FALSE)
  }

  identical(trimws(row[[2L]][[1L]]), pid)
}

.bef_process_start_time <- function(pid = Sys.getpid()) {
  pid <- suppressWarnings(as.integer(pid))
  if (is.na(pid) || pid < 1L) {
    return("unknown-start-time")
  }

  if (identical(.Platform$OS.type, "unix")) {
    linux_start <- .bef_linux_process_start_time(pid)
    if (.bef_lock_start_time_known(linux_start)) {
      return(linux_start)
    }
    ps_start <- .bef_ps_process_start_time(pid)
    if (.bef_lock_start_time_known(ps_start)) {
      return(ps_start)
    }
  }

  "unknown-start-time"
}

.bef_linux_process_start_time <- function(pid) {
  stat_path <- file.path("/proc", as.character(pid), "stat")
  if (!file.exists(stat_path)) {
    return("unknown-start-time")
  }

  stat <- tryCatch(readLines(stat_path, warn = FALSE, n = 1L), error = function(err) "")
  if (!nzchar(stat)) {
    return("unknown-start-time")
  }
  fields <- strsplit(stat, " ", fixed = TRUE)[[1L]]
  if (length(fields) < 22L || !nzchar(fields[[22L]])) {
    return("unknown-start-time")
  }

  paste("linux_ticks", fields[[22L]], sep = ":")
}

.bef_ps_process_start_time <- function(pid) {
  value <- tryCatch(
    system2(
      "ps",
      c("-p", as.character(pid), "-o", "lstart="),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(err) character()
  )
  status <- attr(value, "status", exact = TRUE)
  if (!is.null(status) && !identical(status, 0L)) {
    return("unknown-start-time")
  }
  value <- paste(trimws(value), collapse = " ")
  if (!nzchar(value)) {
    return("unknown-start-time")
  }

  paste("ps_lstart", value, sep = ":")
}

.bef_lock_start_time_known <- function(x) {
  .bef_is_string(x) && !identical(x, "unknown-start-time")
}

.bef_lock_owner_condition <- function(lock_path, owner) {
  rlang::cnd(
    class = "bef_lock_owner",
    message = sprintf("Owner payload for cache lock: %s", lock_path),
    lock_path = lock_path,
    owner = owner
  )
}

.bef_lock_field <- function(owner, field) {
  if (!is.list(owner) || is.null(owner[[field]])) {
    return(NA_character_)
  }
  owner[[field]]
}

.bef_format_lock_time <- function(x) {
  format(x, "%Y-%m-%dT%H:%M:%S%z")
}

.bef_lock_number <- function(x, arg, lower = 0) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower) {
    .bef_abort_lock_invalid(
      sprintf("`%s` must be a finite numeric scalar >= %s.", arg, lower),
      arg = arg,
      predicate = sprintf("finite numeric scalar >= %s", lower)
    )
  }
  as.numeric(x)
}

.bef_lock_string <- function(x, fallback) {
  if (.bef_is_string(x)) {
    as.character(x)
  } else {
    fallback
  }
}

.bef_abort_lock_invalid <- function(message, ..., parent = NULL) {
  .bef_abort_invalid_args(
    message,
    ...,
    module = "cache-lock",
    stage = 6L,
    parent = parent
  )
}

.bef_default_cpp_options <- function() {
  NULL
}

.bef_default_stanc_options <- function() {
  list("O1")
}

#' Clear bayesEfron compilation cache artifacts
#'
#' @description
#' Remove selected in-session and on-disk cache artifacts used by the
#' CmdStan compilation cache. Intended for cache maintenance,
#' stale-lock recovery after an interrupted compile, and forced
#' rebuilds during development.
#'
#' Choose the smallest scope that resolves the issue at hand. The
#' `"lock_only"` scope is the safest and is appropriate for clearing
#' a stale lock left behind by a killed process. `"all"` is the
#' largest scope and will force every subsequent fit in this session
#' (and on this machine, until the cache is repopulated) to recompile
#' from scratch.
#'
#' @details
#' # Scopes
#'
#' | Scope | Removes | Use when |
#' |:------|:--------|:---------|
#' | `"lock_only"` | Stale lock file | Resuming after an interrupted compile that left a lock behind. |
#' | `"session"` | The in-session cache only | Forcing the current session to re-attach to disk-cached binaries. |
#' | `"compiled_models"` | Cache-entry binaries, sidecar JSON, and companion Stan-source copies under the current cache format directory; preserves the cache directory and lock file | Forcing a recompile while keeping the cache root intact. |
#' | `"all"` | The entire cache root, plus the in-session cache | Resetting the cache to a clean slate. |
#'
#' This function does **not** acquire the cache lock. Avoid calling
#' the disk-clearing scopes while another R process is compiling a
#' bayesEfron model in the same cache root.
#'
#' The cache root location is controlled by the environment variable
#' `BAYESEFRON_CACHE_ROOT` (with a sensible per-user default if
#' unset).
#'
#' @param scope Character scalar. One of `"lock_only"` (default),
#'   `"session"`, `"compiled_models"`, or `"all"`.
#'
#' @return Invisibly, a named integer vector with elements
#'   `lock_files`, `session_keys`, `disk_models`, and
#'   `disk_sidecars`. Each element counts the number of artifacts
#'   removed; the contract is fixed even when the scope leaves a
#'   given counter at zero.
#'
#' @seealso
#'   * [bayes_efron_compile()] for repopulating the cache after a
#'     clear.
#'   * [bayes_efron_fit()] for the user-facing fit pipeline that
#'     consumes the cache.
#'
#' @examples
#' # Non-destructive: clear a stale lock if one is present.
#' bayes_efron_clear_cache("lock_only")
#'
#' \dontrun{
#' # Destructive: drop the current session's cache entries.
#' bayes_efron_clear_cache("session")
#'
#' # Most destructive: reset the entire cache root.
#' bayes_efron_clear_cache("all")
#' }
#'
#' @export
bayes_efron_clear_cache <- function(
    scope = c("lock_only", "session", "compiled_models", "all")) {
  scope <- .bef_clear_cache_scope(scope)
  counts <- .bef_clear_cache_counts()

  if (identical(scope, "lock_only")) {
    counts[["lock_files"]] <- .bef_clear_cache_lock()
  } else if (identical(scope, "session")) {
    counts[["session_keys"]] <- as.integer(.bef_cache_clear_session())
  } else if (identical(scope, "compiled_models")) {
    disk_counts <- .bef_clear_compiled_models()
    counts[names(disk_counts)] <- disk_counts
  } else if (identical(scope, "all")) {
    .bef_validate_all_cache_root()
    counts[["session_keys"]] <- as.integer(.bef_cache_clear_session())
    disk_counts <- .bef_clear_all_cache_root()
    counts[names(disk_counts)] <- disk_counts
  }

  message(.bef_clear_cache_message(scope, counts))
  invisible(counts)
}

.bef_clear_cache_scope <- function(scope) {
  choices <- c("lock_only", "session", "compiled_models", "all")
  tryCatch(
    match.arg(scope, choices),
    error = function(err) {
      .bef_abort_invalid_args(
        "`scope` must be one of \"lock_only\", \"session\", \"compiled_models\", or \"all\".",
        arg = "scope",
        predicate = paste(choices, collapse = ", "),
        module = "clear-cache",
        stage = 6L,
        parent = err
      )
    }
  )
}

.bef_clear_cache_counts <- function() {
  c(
    lock_files = 0L,
    session_keys = 0L,
    disk_models = 0L,
    disk_sidecars = 0L
  )
}

.bef_clear_cache_message <- function(scope, counts) {
  sprintf(
    paste0(
      "bayesEfron cache clear (%s): lock_files=%d, session_keys=%d, ",
      "disk_models=%d, disk_sidecars=%d."
    ),
    scope,
    counts[["lock_files"]],
    counts[["session_keys"]],
    counts[["disk_models"]],
    counts[["disk_sidecars"]]
  )
}

.bef_clear_cache_lock <- function(lock_path = .bef_cache_lock_path()) {
  cache_dir <- dirname(lock_path)
  if (dir.exists(cache_dir)) {
    .bef_validate_cache_dir_owner(cache_dir)
  }

  .bef_unlink_cache_paths(
    lock_path,
    artifact = "lock_file",
    count = TRUE,
    recursive = TRUE
  )
}

.bef_clear_compiled_models <- function(cache_dir = .bef_cache_dir(create = FALSE)) {
  counts <- c(disk_models = 0L, disk_sidecars = 0L)
  if (!dir.exists(cache_dir)) {
    return(counts)
  }
  .bef_validate_cache_dir_owner(cache_dir)

  artifacts <- .bef_cache_disk_artifacts(cache_dir)
  counts[["disk_models"]] <- .bef_unlink_cache_paths(
    artifacts$models,
    artifact = "compiled_model",
    count = TRUE
  )
  counts[["disk_sidecars"]] <- .bef_unlink_cache_paths(
    artifacts$sidecars,
    artifact = "cache_sidecar",
    count = TRUE
  )
  .bef_unlink_cache_paths(
    artifacts$sources,
    artifact = "stan_source_copy",
    count = FALSE
  )

  counts
}

.bef_validate_all_cache_root <- function(root = .bef_cache_root()) {
  if (!file.exists(root)) {
    return(invisible(root))
  }
  .bef_validate_cache_dir_owner(root)

  cache_dir <- file.path(root, .bef_cache_format_version())
  if (file.exists(cache_dir)) {
    .bef_validate_cache_dir_owner(cache_dir)
  }
  invisible(root)
}

.bef_clear_all_cache_root <- function(root = .bef_cache_root()) {
  counts <- c(lock_files = 0L, disk_models = 0L, disk_sidecars = 0L)
  if (!file.exists(root)) {
    return(counts)
  }

  cache_dir <- file.path(root, .bef_cache_format_version())
  if (dir.exists(cache_dir)) {
    artifacts <- .bef_cache_disk_artifacts(cache_dir)
    counts[["lock_files"]] <- as.integer(file.exists(artifacts$lock))
    counts[["disk_models"]] <- length(artifacts$models)
    counts[["disk_sidecars"]] <- length(artifacts$sidecars)
  }

  .bef_unlink_cache_paths(
    root,
    artifact = "cache_root",
    count = FALSE,
    recursive = TRUE
  )
  counts
}

.bef_cache_disk_artifacts <- function(cache_dir = .bef_cache_dir(create = FALSE)) {
  if (!dir.exists(cache_dir)) {
    return(list(
      lock = file.path(cache_dir, ".lock"),
      models = character(),
      sidecars = character(),
      sources = character()
    ))
  }

  entries <- list.files(
    cache_dir,
    all.files = TRUE,
    no.. = TRUE,
    full.names = TRUE
  )
  basenames <- basename(entries)
  regular <- .bef_regular_existing_files(entries)

  list(
    lock = file.path(cache_dir, ".lock"),
    models = entries[basenames %in% basenames[.bef_is_cache_model_name(basenames)] &
      entries %in% regular],
    sidecars = entries[grepl("^[0-9a-f]{64}[.]meta[.]json$", basenames) &
      entries %in% regular],
    sources = entries[grepl("^[0-9a-f]{64}[.]stan$", basenames) &
      entries %in% regular]
  )
}

.bef_is_cache_model_name <- function(x) {
  grepl("^[0-9a-f]{64}$", x)
}

.bef_regular_existing_files <- function(paths) {
  if (length(paths) == 0L) {
    return(character())
  }
  info <- file.info(paths)
  paths[!is.na(info$isdir) & !info$isdir]
}

.bef_unlink_cache_paths <- function(paths,
                                    artifact,
                                    count = TRUE,
                                    recursive = FALSE) {
  paths <- unique(as.character(paths))
  paths <- paths[nzchar(paths)]
  if (length(paths) == 0L) {
    return(0L)
  }

  existed <- file.exists(paths)
  targets <- paths[existed]
  if (length(targets) == 0L) {
    return(0L)
  }

  unlink(targets, force = TRUE, recursive = recursive)
  failed <- file.exists(targets)
  if (any(failed)) {
    .bef_abort_cache_perm_violation(
      sprintf("Failed to remove bayesEfron cache artifact `%s`.", artifact),
      artifact = artifact,
      paths = targets[failed]
    )
  }

  if (isTRUE(count)) {
    as.integer(length(targets))
  } else {
    0L
  }
}

.bef_cache_key <- function(stan_file = NULL,
                           model_family = "RE",
                           cpp_options = .bef_default_cpp_options(),
                           stanc_options = .bef_default_stanc_options(),
                           cmdstan_version = NULL,
                           cmdstanr_version = NULL,
                           arch = R.version$arch,
                           compiler = Sys.getenv("CXX", unset = "default"),
                           cache_format_version = .bef_cache_format_version(),
                           makevars_snapshot = NULL,
                           os_major = NULL) {
  stan_file_sha256 <- .bef_stan_file_sha256(stan_file)
  cmdstan_version <- .bef_resolve_cmdstan_version(cmdstan_version)
  cmdstanr_version <- .bef_resolve_cmdstanr_version(cmdstanr_version)
  model_family <- .bef_cache_key_model_family(model_family)
  cache_format_version <- .bef_cache_key_scalar(
    cache_format_version, "cache_format_version", "non-empty character scalar"
  )
  arch <- .bef_cache_key_scalar(arch, "arch", "non-empty character scalar")
  compiler <- .bef_cache_key_scalar(
    compiler, "compiler", "non-empty character scalar"
  )

  cpp_options_normalized <- .bef_normalize_opts(cpp_options)
  stanc_options_normalized <- .bef_normalize_opts(stanc_options)
  makevars_snapshot <- if (is.null(makevars_snapshot)) {
    .bef_makevars_snapshot()
  } else {
    makevars_snapshot
  }
  makevars_normalized <- .bef_normalize_makevars(makevars_snapshot)
  os_major <- .bef_cache_key_scalar(
    if (is.null(os_major)) {
      paste(R.version$platform, .bef_os_major(), sep = "|")
    } else {
      os_major
    },
    "os_major",
    "non-empty character scalar"
  )

  payload <- c(
    stan_file_sha256 = stan_file_sha256,
    cmdstan_version = cmdstan_version,
    cmdstanr_version = cmdstanr_version,
    arch = arch,
    compiler = compiler,
    cache_format_version = cache_format_version,
    model_family = model_family,
    cpp_options_sha256 = .bef_hash_text(cpp_options_normalized),
    stanc_options_sha256 = .bef_hash_text(stanc_options_normalized),
    makevars_sha256 = .bef_hash_text(makevars_normalized),
    os_major = os_major
  )

  payload_json <- .bef_json_payload(payload)
  key <- .bef_hash_text(payload_json)

  list(
    key = key,
    payload = payload,
    cmdstan_v_pre = cmdstan_version,
    payload_json = payload_json,
    provenance = list(
      stan_file = normalizePath(stan_file, winslash = "/", mustWork = TRUE),
      cpp_options = cpp_options,
      stanc_options = stanc_options,
      cpp_options_normalized = cpp_options_normalized,
      stanc_options_normalized = stanc_options_normalized,
      makevars_snapshot = makevars_snapshot,
      makevars_normalized = makevars_normalized,
      os_major = os_major
    )
  )
}

.bef_stan_file_sha256 <- function(stan_file) {
  if (!.bef_is_string(stan_file) ||
      !file.exists(stan_file) ||
      !isFALSE(file.info(stan_file)$isdir)) {
    .bef_abort_cache_key_invalid(
      "`stan_file` must be an existing Stan source file.",
      arg = "stan_file",
      predicate = "existing file path"
    )
  }

  tryCatch(
    digest::digest(stan_file, algo = "sha256", file = TRUE),
    error = function(err) {
      .bef_abort_cache_key_invalid(
        sprintf("Failed to hash Stan source file: %s", stan_file),
        arg = "stan_file",
        predicate = "readable file path",
        parent = err
      )
    }
  )
}

.bef_resolve_cmdstan_version <- function(cmdstan_version) {
  if (!is.null(cmdstan_version)) {
    return(.bef_cache_key_scalar(
      as.character(cmdstan_version),
      "cmdstan_version",
      "non-empty character scalar"
    ))
  }

  rlang::check_installed(
    "cmdstanr",
    reason = "to compile or run the bayesEfron Stan model"
  )
  tryCatch(
    as.character(cmdstanr::cmdstan_version()),
    error = function(err) {
      .bef_abort_compile_failed(
        "Failed to determine the CmdStan version for the cache key.",
        parent = err
      )
    }
  )
}

.bef_resolve_cmdstanr_version <- function(cmdstanr_version) {
  if (!is.null(cmdstanr_version)) {
    return(.bef_cache_key_scalar(
      as.character(cmdstanr_version),
      "cmdstanr_version",
      "non-empty character scalar"
    ))
  }

  rlang::check_installed(
    "cmdstanr",
    reason = "to compile or run the bayesEfron Stan model"
  )
  .bef_cache_key_scalar(
    as.character(utils::packageVersion("cmdstanr")),
    "cmdstanr_version",
    "installed package version"
  )
}

.bef_cache_key_model_family <- function(model_family) {
  model_family <- .bef_cache_key_scalar(
    model_family,
    "model_family",
    '"RE"'
  )
  if (!identical(model_family, "RE")) {
    .bef_abort_cache_key_invalid(
      "`model_family` must be \"RE\" for bayesEfron v0.1.",
      arg = "model_family",
      predicate = '"RE"'
    )
  }
  model_family
}

.bef_normalize_opts <- function(opts) {
  .bef_require_jsonlite()
  if (is.null(opts) || (is.list(opts) && length(opts) == 0L)) {
    return("{}")
  }
  if (!is.list(opts)) {
    .bef_abort_cache_key_invalid(
      "`cpp_options` and `stanc_options` must be NULL or lists.",
      arg = "options",
      predicate = "NULL or list"
    )
  }

  normalized <- .bef_sort_option_list(opts, arg = "options")
  .bef_json(normalized)
}

.bef_sort_option_list <- function(x, arg) {
  nms <- names(x)
  has_names <- !is.null(nms)
  if (has_names) {
    named <- nzchar(nms)
    if (any(named) && !all(named)) {
      .bef_abort_cache_key_invalid(
        "`cpp_options` and `stanc_options` cannot mix named and unnamed entries.",
        arg = arg,
        predicate = "all named or all unnamed"
      )
    }
    if (all(named)) {
      x <- x[order(nms)]
      return(lapply(x, .bef_sort_json_value))
    }
  }

  x <- lapply(x, .bef_sort_json_value)
  x[order(vapply(x, .bef_json, character(1)))]
}

.bef_sort_json_value <- function(x) {
  if (!is.list(x) || inherits(x, "data.frame")) {
    return(x)
  }

  nms <- names(x)
  if (!is.null(nms) && all(nzchar(nms))) {
    x <- x[order(nms)]
    return(lapply(x, .bef_sort_json_value))
  }
  if (!is.null(nms) && any(nzchar(nms))) {
    .bef_abort_cache_key_invalid(
      "Nested cache-key option lists cannot mix named and unnamed entries.",
      arg = "options",
      predicate = "nested lists all named or all unnamed"
    )
  }

  x <- lapply(x, .bef_sort_json_value)
  x[order(vapply(x, .bef_json, character(1)))]
}

.bef_makevars_snapshot <- function(makevars_path = file.path(path.expand("~"), ".R", "Makevars")) {
  env_names <- c(
    "CXXFLAGS",
    "CFLAGS",
    "LDFLAGS",
    "PKG_CXXFLAGS",
    "PKG_CPPFLAGS"
  )
  config_names <- c("CXX17", "CXX17FLAGS", "CXX17PICFLAGS", "CXX17STD")

  env_vars <- Sys.getenv(env_names, unset = "")
  names(env_vars) <- env_names

  list(
    env_vars = as.list(env_vars),
    makevars_bytes_sha256 = .bef_makevars_file_sha256(makevars_path),
    r_cmd_config = as.list(.bef_r_cmd_config(config_names))
  )
}

.bef_makevars_file_sha256 <- function(path) {
  if (!file.exists(path) || isTRUE(file.info(path)$isdir)) {
    return("absent")
  }
  tryCatch(
    digest::digest(path, algo = "sha256", file = TRUE),
    error = function(err) "unavailable"
  )
}

.bef_r_cmd_config <- function(keys) {
  out <- vapply(
    keys,
    function(key) {
      value <- tryCatch(
        suppressWarnings(system2(
          file.path(R.home("bin"), "R"),
          c("CMD", "config", key),
          stdout = TRUE,
          stderr = FALSE
        )),
        error = function(err) "unavailable"
      )
      status <- attr(value, "status", exact = TRUE)
      if (!is.null(status) && !identical(status, 0L)) {
        return("unavailable")
      }
      if (length(value) == 0L || all(!nzchar(value))) {
        return("unavailable")
      }
      paste(value, collapse = "\n")
    },
    character(1)
  )
  names(out) <- keys
  out
}

.bef_normalize_makevars <- function(snapshot = NULL) {
  .bef_require_jsonlite()
  if (is.null(snapshot)) {
    snapshot <- .bef_makevars_snapshot()
  }
  if (!is.list(snapshot)) {
    .bef_abort_cache_key_invalid(
      "`makevars_snapshot` must be NULL or a list.",
      arg = "makevars_snapshot",
      predicate = "NULL or list"
    )
  }

  required <- c("env_vars", "makevars_bytes_sha256", "r_cmd_config")
  missing <- setdiff(required, names(snapshot))
  if (length(missing) > 0L) {
    .bef_abort_cache_key_invalid(
      "`makevars_snapshot` is missing required fields.",
      arg = "makevars_snapshot",
      predicate = paste(required, collapse = ", "),
      missing_fields = missing
    )
  }

  .bef_json(.bef_sort_json_value(snapshot))
}

.bef_os_major <- function() {
  sysname <- Sys.info()[["sysname"]]
  if (is.na(sysname) || !nzchar(sysname)) {
    sysname <- .Platform$OS.type
  }

  major <- "unknown"
  if (identical(sysname, "Darwin")) {
    product_version <- tryCatch(
      system2("sw_vers", "-productVersion", stdout = TRUE, stderr = FALSE),
      error = function(err) character()
    )
    if (length(product_version) > 0L && nzchar(product_version[[1L]])) {
      major <- strsplit(product_version[[1L]], ".", fixed = TRUE)[[1L]][[1L]]
    }
  } else if (identical(.Platform$OS.type, "windows")) {
    release <- Sys.info()[["release"]]
    if (!is.na(release) && nzchar(release)) {
      major <- sub("^([0-9]+).*$", "\\1", release)
    }
  } else {
    os <- R.version$os
    if (!is.null(os) && nzchar(os) && grepl("[0-9]", os)) {
      major <- sub("^.*?([0-9]+).*$", "\\1", os)
    }
  }

  paste(sysname, major, sep = "-")
}

.bef_json_payload <- function(payload) {
  .bef_json(as.list(payload))
}

.bef_json <- function(x) {
  .bef_require_jsonlite()
  as.character(jsonlite::toJSON(
    x,
    auto_unbox = TRUE,
    null = "null",
    digits = NA
  ))
}

.bef_hash_text <- function(x) {
  digest::digest(x, algo = "sha256", serialize = FALSE)
}

.bef_require_jsonlite <- function() {
  rlang::check_installed(
    "jsonlite",
    reason = "to compute deterministic bayesEfron cache keys and validate cache sidecars"
  )
}

.bef_cache_key_scalar <- function(x, arg, predicate) {
  if (!.bef_is_string(x)) {
    .bef_abort_cache_key_invalid(
      sprintf("`%s` must be a %s.", arg, predicate),
      arg = arg,
      predicate = predicate
    )
  }
  as.character(x)
}

.bef_abort_cache_key_invalid <- function(message, ..., parent = NULL) {
  .bef_abort_invalid_args(
    message,
    ...,
    module = "cache-key",
    stage = 6L,
    parent = parent
  )
}

.bef_model <- function(model_name = "RE",
                       cpp_options = .bef_default_cpp_options(),
                       stanc_options = .bef_default_stanc_options(),
                       force_recompile = FALSE,
                       stan_file = NULL,
                       cache_dir = .bef_cache_dir(create = TRUE),
                       lock_timeout = .bef_cache_lock_timeout(),
                       key_info = NULL,
                       cmdstan_model_fun = .bef_cmdstan_model,
                       cmdstan_version_fun = .bef_cmdstan_version,
                       acquire_lock_fun = .bef_acquire_lock,
                       release_lock_fun = .bef_release_lock,
                       check_installed = TRUE,
                       cache_key_args = list()) {
  force_recompile <- .bef_model_logical(force_recompile, "force_recompile")
  model_family <- .bef_model_family(model_name)
  if (isTRUE(check_installed)) {
    .bef_check_cmdstanr_installed()
  }

  if (is.null(stan_file)) {
    stan_file <- .bef_stan_file(model_family)
  }
  key_info <- if (is.null(key_info)) {
    do.call(
      .bef_cache_key,
      c(
        list(
          stan_file = stan_file,
          model_family = model_family,
          cpp_options = cpp_options,
          stanc_options = stanc_options
        ),
        cache_key_args
      )
    )
  } else {
    .bef_validate_key_info(key_info)
  }
  key <- key_info$key

  if (!isTRUE(force_recompile) && .bef_cache_exists(key)) {
    return(.bef_cache_get(key))
  }
  if (isTRUE(force_recompile)) {
    .bef_cache_remove(key)
  }

  cache_dir <- .bef_ensure_cache_dir(cache_dir)
  paths <- .bef_cache_entry_paths(key, cache_dir = cache_dir)
  lock_path <- .bef_cache_lock_path_for_dir(cache_dir)
  acquire_lock_fun(timeout = lock_timeout, lock_path = lock_path)
  on.exit(release_lock_fun(lock_path), add = TRUE, after = FALSE)

  if (!isTRUE(force_recompile) && .bef_cache_exists(key)) {
    return(.bef_cache_get(key))
  }

  if (isTRUE(force_recompile)) {
    .bef_delete_cache_entry(paths)
  } else {
    disk_model <- .bef_try_disk_model(
      key_info = key_info,
      paths = paths,
      cmdstan_model_fun = cmdstan_model_fun
    )
    if (!is.null(disk_model)) {
      .bef_cache_set(key, disk_model)
      return(disk_model)
    }
  }

  model <- .bef_compile_model(
    key_info = key_info,
    paths = paths,
    cpp_options = cpp_options,
    stanc_options = stanc_options,
    cmdstan_model_fun = cmdstan_model_fun,
    cmdstan_version_fun = cmdstan_version_fun
  )
  .bef_cache_set(key, model)
  model
}

.bef_try_disk_model <- function(key_info, paths, cmdstan_model_fun) {
  if (!file.exists(paths$meta)) {
    return(NULL)
  }

  read <- .bef_read_sidecar(paths$meta)
  if (!isTRUE(read$valid)) {
    .bef_delete_cache_entry(paths)
    return(NULL)
  }

  validation <- tryCatch(
    .bef_validate_sidecar(
      read$meta,
      key_info = key_info,
      paths = paths,
      warn_corruption = TRUE
    ),
    bef_cache_format_mismatch = function(err) {
      .bef_delete_cache_entry(paths)
      NULL
    }
  )
  if (is.null(validation) || !isTRUE(validation$valid)) {
    .bef_delete_cache_entry(paths)
    return(NULL)
  }

  tryCatch(
    cmdstan_model_fun(
      stan_file = paths$stan,
      compile = FALSE,
      exe_file = paths$exe
    ),
    error = function(err) {
      .bef_delete_cache_entry(paths)
      NULL
    }
  )
}

.bef_compile_model <- function(key_info,
                               paths,
                               cpp_options,
                               stanc_options,
                               cmdstan_model_fun,
                               cmdstan_version_fun) {
  staging_dir <- .bef_staging_dir(paths$cache_dir)
  unlink_staging <- TRUE
  on.exit({
    if (isTRUE(unlink_staging)) {
      unlink(staging_dir, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)

  stage_paths <- .bef_stage_entry_paths(key_info$key, staging_dir)
  .bef_copy_file(
    key_info$provenance$stan_file,
    stage_paths$stan,
    artifact = "stan_source"
  )

  started <- Sys.time()
  model <- tryCatch(
    cmdstan_model_fun(
      stan_file = stage_paths$stan,
      compile = TRUE,
      dir = staging_dir,
      cpp_options = cpp_options,
      stanc_options = stanc_options
    ),
    error = function(err) {
      .bef_abort_compile_failed(
        "CmdStan failed to compile the bayesEfron Stan model.",
        key = key_info$key,
        stan_file = key_info$provenance$stan_file,
        parent = err
      )
    }
  )
  compile_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  .bef_normalize_stage_executable(model, stage_paths)
  if (!file.exists(stage_paths$exe) || !isFALSE(file.info(stage_paths$exe)$isdir)) {
    .bef_abort_compile_failed(
      "CmdStan compilation did not create the expected bayesEfron executable.",
      key = key_info$key,
      expected_exe = stage_paths$exe
    )
  }

  cmdstan_v_post <- tryCatch(
    .bef_sidecar_string(cmdstan_version_fun(), "cmdstan_v_post"),
    error = function(err) {
      .bef_abort_compile_failed(
        "Failed to re-capture CmdStan version after compilation.",
        key = key_info$key,
        parent = err
      )
    }
  )
  if (!identical(key_info$cmdstan_v_pre, cmdstan_v_post)) {
    .bef_abort_cache_version_drift(
      "CmdStan version changed between cache-key capture and compilation.",
      key = key_info$key,
      cmdstan_v_pre = key_info$cmdstan_v_pre,
      cmdstan_v_post = cmdstan_v_post
    )
  }

  meta <- .bef_build_sidecar(
    key_info,
    stage_paths$exe,
    compiled_at = Sys.time(),
    cmdstan_v_post = cmdstan_v_post,
    compile_seconds = compile_seconds
  )
  .bef_write_sidecar(meta, stage_paths$meta)
  .bef_promote_staged_entry(stage_paths, paths)
  unlink(staging_dir, recursive = TRUE, force = TRUE)
  unlink_staging <- FALSE

  tryCatch(
    cmdstan_model_fun(
      stan_file = paths$stan,
      compile = FALSE,
      exe_file = paths$exe
    ),
    error = function(err) {
      .bef_delete_cache_entry(paths)
      .bef_abort_compile_failed(
        "Compiled bayesEfron model could not be reattached from the final cache path.",
        key = key_info$key,
        stan_file = paths$stan,
        exe_file = paths$exe,
        parent = err
      )
    }
  )
}

.bef_copy_file <- function(from, to, artifact) {
  ok <- file.copy(from, to, overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
  if (!isTRUE(ok)) {
    .bef_abort_cache_perm_violation(
      sprintf("Failed to stage bayesEfron cache artifact `%s`.", artifact),
      artifact = artifact,
      from = from,
      to = to
    )
  }
  invisible(to)
}

.bef_normalize_stage_executable <- function(model, stage_paths) {
  if (file.exists(stage_paths$exe) && isFALSE(file.info(stage_paths$exe)$isdir)) {
    return(invisible(stage_paths$exe))
  }

  exe_file <- .bef_model_exe_file(model)
  if (is.null(exe_file) ||
      !file.exists(exe_file) ||
      !isFALSE(file.info(exe_file)$isdir)) {
    return(invisible(FALSE))
  }

  if (identical(
    normalizePath(exe_file, winslash = "/", mustWork = TRUE),
    normalizePath(stage_paths$exe, winslash = "/", mustWork = FALSE)
  )) {
    return(invisible(stage_paths$exe))
  }

  .bef_copy_file(exe_file, stage_paths$exe, artifact = "compiled_executable")
}

.bef_model_exe_file <- function(model) {
  if (!is.null(model) && is.function(model$exe_file)) {
    value <- tryCatch(model$exe_file(), error = function(err) NULL)
    if (.bef_is_string(value)) {
      return(as.character(value))
    }
  }
  NULL
}

.bef_cache_lock_path_for_dir <- function(cache_dir) {
  file.path(.bef_ensure_cache_dir(cache_dir), ".lock")
}

.bef_cmdstan_model <- function(...) {
  cmdstanr::cmdstan_model(...)
}

.bef_cmdstan_version <- function() {
  as.character(cmdstanr::cmdstan_version())
}

.bef_check_cmdstanr_installed <- function() {
  rlang::check_installed(
    "cmdstanr",
    reason = "to compile or run the bayesEfron Stan model"
  )
}

.bef_stan_file <- function(model_family = "RE") {
  model_family <- .bef_model_family(model_family)
  file_name <- switch(model_family, RE = "efron_re.stan")
  path <- system.file("stan", file_name, package = "bayesEfron", mustWork = FALSE)
  if (!nzchar(path)) {
    path <- file.path("inst", "stan", file_name)
  }
  if (!file.exists(path) || !isFALSE(file.info(path)$isdir)) {
    .bef_abort_compile_failed(
      sprintf("Cannot locate bayesEfron Stan source file `%s`.", file_name),
      model_family = model_family,
      stan_file = path
    )
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.bef_model_family <- function(model_name) {
  if (!.bef_is_string(model_name)) {
    .bef_abort_model_cache_invalid(
      "`model_name` must be a non-empty character scalar.",
      arg = "model_name",
      predicate = '"RE" or "efron_re"'
    )
  }
  model_name <- as.character(model_name)
  if (identical(model_name, "RE") || identical(model_name, "efron_re")) {
    return("RE")
  }

  .bef_abort_model_cache_invalid(
    "`model_name` must be \"RE\" for bayesEfron v0.1.",
    arg = "model_name",
    predicate = '"RE" or "efron_re"'
  )
}

.bef_model_logical <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .bef_abort_model_cache_invalid(
      sprintf("`%s` must be TRUE or FALSE.", arg),
      arg = arg,
      predicate = "logical scalar"
    )
  }
  isTRUE(x)
}

.bef_abort_model_cache_invalid <- function(message, ..., parent = NULL) {
  .bef_abort_invalid_args(
    message,
    ...,
    module = "model-cache",
    stage = 6L,
    parent = parent
  )
}

.bef_cache_entry_paths <- function(key,
                                   cache_dir = .bef_cache_dir(create = TRUE)) {
  key <- .bef_cache_key_hex(key)
  cache_dir <- .bef_ensure_cache_dir(cache_dir)

  list(
    key = key,
    cache_dir = cache_dir,
    stan = file.path(cache_dir, paste0(key, ".stan")),
    exe = file.path(cache_dir, key),
    meta = file.path(cache_dir, paste0(key, ".meta.json"))
  )
}

.bef_sidecar_path <- function(key, cache_dir = .bef_cache_dir(create = FALSE)) {
  file.path(cache_dir, paste0(.bef_cache_key_hex(key), ".meta.json"))
}

.bef_staging_dir <- function(cache_dir = .bef_cache_dir(create = TRUE)) {
  cache_dir <- .bef_ensure_cache_dir(cache_dir)
  for (i in seq_len(100L)) {
    path <- tempfile(".staging-", tmpdir = cache_dir)
    if (!dir.exists(path)) {
      return(.bef_ensure_cache_dir(path))
    }
  }

  .bef_abort_cache_perm_violation(
    "Failed to allocate a unique bayesEfron cache staging directory.",
    cache_dir = cache_dir
  )
}

.bef_stage_entry_paths <- function(key, staging_dir) {
  key <- .bef_cache_key_hex(key)
  staging_dir <- .bef_ensure_cache_dir(staging_dir)

  list(
    key = key,
    staging_dir = staging_dir,
    stan = file.path(staging_dir, paste0(key, ".stan")),
    exe = file.path(staging_dir, key),
    meta = file.path(staging_dir, paste0(key, ".meta.json"))
  )
}

.bef_build_sidecar <- function(key_info,
                               binary_path,
                               compiled_at = Sys.time(),
                               cmdstan_v_post = NULL,
                               compile_seconds = NULL,
                               host = Sys.info()[["nodename"]],
                               r_version = R.version.string) {
  key_info <- .bef_validate_key_info(key_info)
  binary_sha256 <- .bef_binary_sha256(binary_path)
  payload <- key_info$payload
  provenance <- key_info$provenance
  if (is.null(cmdstan_v_post)) {
    cmdstan_v_post <- key_info$cmdstan_v_pre
  }

  meta <- list(
    cache_format_version = payload[["cache_format_version"]],
    stan_file_sha256 = payload[["stan_file_sha256"]],
    cmdstan_version = .bef_sidecar_string(cmdstan_v_post, "cmdstan_v_post"),
    cmdstanr_version = payload[["cmdstanr_version"]],
    arch = payload[["arch"]],
    compiler = payload[["compiler"]],
    compiled_at = .bef_format_sidecar_time(compiled_at),
    cpp_options = .bef_sidecar_cpp_options(provenance$cpp_options),
    stanc_options = .bef_sidecar_stanc_options(provenance$stanc_options),
    makevars_snapshot = provenance$makevars_snapshot,
    os_major = provenance$os_major,
    binary_sha256 = binary_sha256
  )

  if (!is.null(host)) {
    meta$host <- .bef_sidecar_string(host, "host")
  }
  if (!is.null(r_version)) {
    meta$r_version <- .bef_sidecar_string(r_version, "r_version")
  }
  if (!is.null(compile_seconds)) {
    meta$compile_seconds <- .bef_sidecar_number(
      compile_seconds, "compile_seconds", lower = 0
    )
  }

  meta
}

.bef_write_sidecar <- function(meta, meta_path) {
  .bef_require_jsonlite()
  if (!is.list(meta)) {
    .bef_abort_sidecar_invalid(
      "`meta` must be a sidecar metadata list.",
      arg = "meta",
      predicate = "list"
    )
  }
  if (!.bef_is_string(meta_path)) {
    .bef_abort_sidecar_invalid(
      "`meta_path` must be a non-empty character scalar.",
      arg = "meta_path",
      predicate = "non-empty character scalar"
    )
  }

  dir_path <- .bef_ensure_cache_dir(dirname(meta_path))
  tmp <- tempfile(".meta-", tmpdir = dir_path, fileext = ".json")
  on.exit(unlink(tmp, force = TRUE), add = TRUE)

  tryCatch(
    jsonlite::write_json(
      meta,
      tmp,
      auto_unbox = TRUE,
      null = "null",
      pretty = TRUE
    ),
    error = function(err) {
      .bef_abort_cache_perm_violation(
        sprintf("Failed to write bayesEfron cache sidecar: %s", meta_path),
        meta_path = meta_path,
        parent = err
      )
    }
  )

  ok <- suppressWarnings(file.rename(tmp, meta_path))
  if (!isTRUE(ok)) {
    unlink(meta_path, force = TRUE)
    ok <- suppressWarnings(file.rename(tmp, meta_path))
  }
  if (!isTRUE(ok)) {
    .bef_abort_cache_perm_violation(
      sprintf("Failed to install bayesEfron cache sidecar: %s", meta_path),
      meta_path = meta_path
    )
  }

  invisible(normalizePath(meta_path, winslash = "/", mustWork = TRUE))
}

.bef_read_sidecar <- function(meta_path) {
  .bef_require_jsonlite()
  if (!.bef_is_string(meta_path)) {
    .bef_abort_sidecar_invalid(
      "`meta_path` must be a non-empty character scalar.",
      arg = "meta_path",
      predicate = "non-empty character scalar"
    )
  }
  if (!file.exists(meta_path) || !isFALSE(file.info(meta_path)$isdir)) {
    return(.bef_sidecar_result(
      FALSE,
      reason = "missing",
      meta_path = meta_path,
      message = "Sidecar file is missing."
    ))
  }

  tryCatch(
    {
      meta <- jsonlite::read_json(meta_path, simplifyVector = FALSE)
      .bef_sidecar_result(TRUE, meta = meta, meta_path = meta_path)
    },
    error = function(err) {
      .bef_sidecar_result(
        FALSE,
        reason = "parse",
        meta_path = meta_path,
        message = "Sidecar JSON could not be parsed.",
        parent = err
      )
    }
  )
}

.bef_validate_sidecar <- function(meta,
                                  key_info = NULL,
                                  paths = NULL,
                                  expected_format = .bef_cache_format_version(),
                                  warn_corruption = TRUE) {
  schema <- .bef_validate_sidecar_schema(meta, expected_format = expected_format)
  if (!isTRUE(schema$valid)) {
    return(schema)
  }

  if (!is.null(key_info)) {
    key_info <- .bef_validate_key_info(key_info)
    identity <- .bef_validate_sidecar_identity(meta, key_info)
    if (!isTRUE(identity$valid)) {
      return(identity)
    }
  }

  if (!is.null(paths)) {
    paths <- .bef_validate_cache_paths(paths)
    source_result <- .bef_validate_sidecar_source(meta, paths)
    if (!isTRUE(source_result$valid)) {
      return(source_result)
    }

    binary_result <- .bef_validate_sidecar_binary(
      meta,
      paths,
      warn_corruption = warn_corruption
    )
    if (!isTRUE(binary_result$valid)) {
      return(binary_result)
    }
  }

  .bef_sidecar_result(TRUE, meta = meta, paths = paths)
}

.bef_validate_sidecar_schema <- function(meta,
                                         expected_format = .bef_cache_format_version()) {
  if (!is.list(meta)) {
    return(.bef_sidecar_result(
      FALSE,
      reason = "schema",
      message = "Sidecar metadata is not a JSON object.",
      field = NA_character_
    ))
  }

  missing <- setdiff(.bef_sidecar_required_fields(), names(meta))
  if (length(missing) > 0L) {
    return(.bef_sidecar_result(
      FALSE,
      reason = "schema",
      message = "Sidecar metadata is missing required fields.",
      missing_fields = missing
    ))
  }

  scalar_fields <- c(
    "cache_format_version",
    "stan_file_sha256",
    "cmdstan_version",
    "cmdstanr_version",
    "arch",
    "compiler",
    "compiled_at",
    "os_major",
    "binary_sha256"
  )
  for (field in scalar_fields) {
    if (!.bef_is_string(meta[[field]])) {
      return(.bef_sidecar_result(
        FALSE,
        reason = "schema",
        message = sprintf("Sidecar field `%s` must be a string.", field),
        field = field
      ))
    }
  }

  for (field in c("stan_file_sha256", "binary_sha256")) {
    if (!.bef_is_sha256_hex(meta[[field]])) {
      return(.bef_sidecar_result(
        FALSE,
        reason = "schema",
        message = sprintf("Sidecar field `%s` must be a SHA-256 hex digest.", field),
        field = field
      ))
    }
  }

  if (!identical(meta$cache_format_version, expected_format)) {
    .bef_abort_cache_format_mismatch(
      sprintf(
        "bayesEfron cache sidecar format `%s` does not match expected `%s`.",
        meta$cache_format_version,
        expected_format
      ),
      cache_format_version = meta$cache_format_version,
      expected_cache_format_version = expected_format
    )
  }

  if (!is.list(meta$cpp_options)) {
    return(.bef_sidecar_result(
      FALSE,
      reason = "schema",
      message = "Sidecar field `cpp_options` must be a JSON object or list.",
      field = "cpp_options"
    ))
  }
  if (!is.list(meta$stanc_options)) {
    return(.bef_sidecar_result(
      FALSE,
      reason = "schema",
      message = "Sidecar field `stanc_options` must be a JSON array/list.",
      field = "stanc_options"
    ))
  }

  makevars <- .bef_validate_sidecar_makevars(meta$makevars_snapshot)
  if (!isTRUE(makevars$valid)) {
    return(makevars)
  }

  if (!is.null(meta$host) && !.bef_is_string(meta$host)) {
    return(.bef_sidecar_result(
      FALSE,
      reason = "schema",
      message = "Sidecar field `host` must be a string when present.",
      field = "host"
    ))
  }
  if (!is.null(meta$r_version) && !.bef_is_string(meta$r_version)) {
    return(.bef_sidecar_result(
      FALSE,
      reason = "schema",
      message = "Sidecar field `r_version` must be a string when present.",
      field = "r_version"
    ))
  }
  if (!is.null(meta$compile_seconds) &&
      (!is.numeric(meta$compile_seconds) ||
        length(meta$compile_seconds) != 1L ||
        !is.finite(meta$compile_seconds) ||
        meta$compile_seconds < 0)) {
    return(.bef_sidecar_result(
      FALSE,
      reason = "schema",
      message = "Sidecar field `compile_seconds` must be a non-negative number.",
      field = "compile_seconds"
    ))
  }

  .bef_sidecar_result(TRUE, meta = meta)
}

.bef_validate_sidecar_identity <- function(meta, key_info) {
  payload <- key_info$payload
  scalar_pairs <- c(
    stan_file_sha256 = "stan_file_sha256",
    cmdstan_version = "cmdstan_version",
    cmdstanr_version = "cmdstanr_version",
    arch = "arch",
    compiler = "compiler",
    cache_format_version = "cache_format_version",
    os_major = "os_major"
  )
  for (field in names(scalar_pairs)) {
    payload_field <- scalar_pairs[[field]]
    if (!identical(meta[[field]], payload[[payload_field]])) {
      return(.bef_sidecar_result(
        FALSE,
        reason = "payload_mismatch",
        message = sprintf("Sidecar field `%s` does not match the cache key payload.", field),
        field = field,
        recorded = meta[[field]],
        expected = payload[[payload_field]]
      ))
    }
  }

  option_pairs <- list(
    cpp_options = c(meta_field = "cpp_options", payload_field = "cpp_options_sha256"),
    stanc_options = c(meta_field = "stanc_options", payload_field = "stanc_options_sha256")
  )
  for (pair in option_pairs) {
    computed <- .bef_hash_text(.bef_normalize_opts(meta[[pair[["meta_field"]]]]))
    expected <- payload[[pair[["payload_field"]]]]
    if (!identical(computed, expected)) {
      return(.bef_sidecar_result(
        FALSE,
        reason = "payload_mismatch",
        message = sprintf(
          "Sidecar field `%s` does not match the cache key payload.",
          pair[["meta_field"]]
        ),
        field = pair[["meta_field"]],
        recorded_sha256 = computed,
        expected_sha256 = expected
      ))
    }
  }

  makevars_sha256 <- .bef_hash_text(.bef_normalize_makevars(meta$makevars_snapshot))
  if (!identical(makevars_sha256, payload[["makevars_sha256"]])) {
    return(.bef_sidecar_result(
      FALSE,
      reason = "payload_mismatch",
      message = "Sidecar field `makevars_snapshot` does not match the cache key payload.",
      field = "makevars_snapshot",
      recorded_sha256 = makevars_sha256,
      expected_sha256 = payload[["makevars_sha256"]]
    ))
  }

  .bef_sidecar_result(TRUE, meta = meta)
}

.bef_validate_sidecar_source <- function(meta, paths) {
  if (!file.exists(paths$stan) || !isFALSE(file.info(paths$stan)$isdir)) {
    return(.bef_sidecar_result(
      FALSE,
      reason = "stan_source_missing",
      message = "Cached Stan source copy is missing.",
      stan_path = paths$stan
    ))
  }

  computed <- .bef_binary_sha256(paths$stan)
  if (!identical(computed, meta$stan_file_sha256)) {
    return(.bef_sidecar_result(
      FALSE,
      reason = "stan_file_sha256_mismatch",
      message = "Cached Stan source copy does not match sidecar digest.",
      stan_path = paths$stan,
      recorded_sha256 = meta$stan_file_sha256,
      computed_sha256 = computed
    ))
  }

  .bef_sidecar_result(TRUE, meta = meta)
}

.bef_validate_sidecar_binary <- function(meta, paths, warn_corruption = TRUE) {
  if (!file.exists(paths$exe) || !isFALSE(file.info(paths$exe)$isdir)) {
    return(.bef_sidecar_result(
      FALSE,
      reason = "binary_missing",
      message = "Cached CmdStan executable is missing.",
      binary_path = paths$exe
    ))
  }

  computed <- .bef_binary_sha256(paths$exe)
  if (!identical(computed, meta$binary_sha256)) {
    if (isTRUE(warn_corruption)) {
      .bef_warn_cache_corruption(
        "Cached bayesEfron CmdStan executable failed SHA-256 verification.",
        recorded_sha256 = meta$binary_sha256,
        computed_sha256 = computed,
        key = paths$key,
        cache_dir = paths$cache_dir,
        binary_path = paths$exe
      )
    }
    return(.bef_sidecar_result(
      FALSE,
      reason = "binary_sha256_mismatch",
      message = "Cached CmdStan executable does not match sidecar digest.",
      binary_path = paths$exe,
      recorded_sha256 = meta$binary_sha256,
      computed_sha256 = computed
    ))
  }

  .bef_sidecar_result(TRUE, meta = meta)
}

.bef_delete_cache_entry <- function(paths) {
  paths <- if (.bef_is_string(paths) && .bef_is_sha256_hex(paths)) {
    .bef_cache_entry_paths(paths)
  } else {
    .bef_validate_cache_paths(paths)
  }

  targets <- c(stan = paths$stan, exe = paths$exe, meta = paths$meta)
  existed <- file.exists(targets)
  unlink(targets, force = TRUE, recursive = TRUE)
  removed <- existed & !file.exists(targets)
  invisible(stats::setNames(as.integer(removed), names(targets)))
}

.bef_promote_staged_entry <- function(stage_paths, entry_paths) {
  stage_paths <- .bef_validate_stage_paths(stage_paths)
  entry_paths <- .bef_validate_cache_paths(entry_paths)
  .bef_ensure_cache_dir(entry_paths$cache_dir)

  promoted <- character()
  for (name in c("meta", "stan", "exe")) {
    if (!file.exists(stage_paths[[name]])) {
      unlink(promoted, force = TRUE, recursive = TRUE)
      .bef_delete_cache_entry(entry_paths)
      .bef_abort_cache_perm_violation(
        sprintf("Missing staged cache artifact `%s`: %s", name, stage_paths[[name]]),
        artifact = name,
        path = stage_paths[[name]]
      )
    }
    ok <- suppressWarnings(file.rename(stage_paths[[name]], entry_paths[[name]]))
    if (!isTRUE(ok)) {
      unlink(entry_paths[[name]], force = TRUE)
      ok <- suppressWarnings(file.rename(stage_paths[[name]], entry_paths[[name]]))
    }
    if (!isTRUE(ok)) {
      unlink(promoted, force = TRUE, recursive = TRUE)
      .bef_delete_cache_entry(entry_paths)
      .bef_abort_cache_perm_violation(
        sprintf("Failed to promote staged cache artifact `%s`.", name),
        artifact = name,
        from = stage_paths[[name]],
        to = entry_paths[[name]]
      )
    }
    promoted <- c(promoted, entry_paths[[name]])
  }

  invisible(entry_paths)
}

.bef_binary_sha256 <- function(path) {
  if (!.bef_is_string(path) ||
      !file.exists(path) ||
      !isFALSE(file.info(path)$isdir)) {
    .bef_abort_sidecar_invalid(
      "`path` must be an existing file path.",
      arg = "path",
      predicate = "existing file path"
    )
  }

  tryCatch(
    digest::digest(path, algo = "sha256", file = TRUE),
    error = function(err) {
      .bef_abort_sidecar_invalid(
        sprintf("Failed to hash file contents: %s", path),
        arg = "path",
        predicate = "readable file path",
        parent = err
      )
    }
  )
}

.bef_validate_cache_dir_owner <- function(path = .bef_cache_dir(create = TRUE)) {
  if (identical(.Platform$OS.type, "windows")) {
    return(invisible(path))
  }
  if (!exists("Sys.getuid", envir = baseenv(), mode = "function")) {
    return(invisible(path))
  }

  info <- file.info(path)
  uid <- info[["uid"]]
  current_uid <- get("Sys.getuid", envir = baseenv(), mode = "function")()
  if (length(uid) != 1L ||
      is.na(uid) ||
      is.na(current_uid) ||
      identical(as.integer(uid), as.integer(current_uid))) {
    return(invisible(path))
  }

  .bef_abort_cache_perm_violation(
    sprintf("bayesEfron cache directory is owned by uid %s, not current uid %s.",
            uid, current_uid),
    cache_dir = path,
    owner_uid = uid,
    current_uid = current_uid
  )
}

.bef_validate_key_info <- function(key_info) {
  if (!is.list(key_info)) {
    .bef_abort_sidecar_invalid(
      "`key_info` must be the list returned by `.bef_cache_key()`.",
      arg = "key_info",
      predicate = "cache-key result list"
    )
  }

  required_top <- c("key", "payload", "cmdstan_v_pre", "provenance")
  missing_top <- setdiff(required_top, names(key_info))
  if (length(missing_top) > 0L) {
    .bef_abort_sidecar_invalid(
      "`key_info` is missing required fields.",
      arg = "key_info",
      predicate = paste(required_top, collapse = ", "),
      missing_fields = missing_top
    )
  }
  key_info$key <- .bef_cache_key_hex(key_info$key)

  required_payload <- .bef_cache_key_payload_fields()
  missing_payload <- setdiff(required_payload, names(key_info$payload))
  if (length(missing_payload) > 0L) {
    .bef_abort_sidecar_invalid(
      "`key_info$payload` is missing required fields.",
      arg = "key_info$payload",
      predicate = paste(required_payload, collapse = ", "),
      missing_fields = missing_payload
    )
  }

  required_provenance <- c("cpp_options", "stanc_options", "makevars_snapshot", "os_major")
  missing_provenance <- setdiff(required_provenance, names(key_info$provenance))
  if (length(missing_provenance) > 0L) {
    .bef_abort_sidecar_invalid(
      "`key_info$provenance` is missing required fields.",
      arg = "key_info$provenance",
      predicate = paste(required_provenance, collapse = ", "),
      missing_fields = missing_provenance
    )
  }

  key_info
}

.bef_validate_cache_paths <- function(paths) {
  if (!is.list(paths)) {
    .bef_abort_sidecar_invalid(
      "`paths` must be a cache path list.",
      arg = "paths",
      predicate = "cache path list"
    )
  }

  required <- c("key", "cache_dir", "stan", "exe", "meta")
  missing <- setdiff(required, names(paths))
  if (length(missing) > 0L) {
    .bef_abort_sidecar_invalid(
      "`paths` is missing required fields.",
      arg = "paths",
      predicate = paste(required, collapse = ", "),
      missing_fields = missing
    )
  }
  paths$key <- .bef_cache_key_hex(paths$key)
  for (field in required[-1L]) {
    paths[[field]] <- .bef_sidecar_string(paths[[field]], field)
  }
  paths$cache_dir <- normalizePath(paths$cache_dir, winslash = "/", mustWork = FALSE)
  paths
}

.bef_validate_stage_paths <- function(paths) {
  if (!is.list(paths)) {
    .bef_abort_sidecar_invalid(
      "`stage_paths` must be a staging path list.",
      arg = "stage_paths",
      predicate = "staging path list"
    )
  }

  required <- c("key", "staging_dir", "stan", "exe", "meta")
  missing <- setdiff(required, names(paths))
  if (length(missing) > 0L) {
    .bef_abort_sidecar_invalid(
      "`stage_paths` is missing required fields.",
      arg = "stage_paths",
      predicate = paste(required, collapse = ", "),
      missing_fields = missing
    )
  }
  paths$key <- .bef_cache_key_hex(paths$key)
  for (field in required[-1L]) {
    paths[[field]] <- .bef_sidecar_string(paths[[field]], field)
  }
  paths
}

.bef_validate_sidecar_makevars <- function(makevars) {
  if (!is.list(makevars)) {
    return(.bef_sidecar_result(
      FALSE,
      reason = "schema",
      message = "Sidecar field `makevars_snapshot` must be an object.",
      field = "makevars_snapshot"
    ))
  }

  required <- c("env_vars", "makevars_bytes_sha256", "r_cmd_config")
  missing <- setdiff(required, names(makevars))
  if (length(missing) > 0L) {
    return(.bef_sidecar_result(
      FALSE,
      reason = "schema",
      message = "Sidecar field `makevars_snapshot` is missing required fields.",
      field = "makevars_snapshot",
      missing_fields = missing
    ))
  }
  if (!is.list(makevars$env_vars)) {
    return(.bef_sidecar_result(
      FALSE,
      reason = "schema",
      message = "`makevars_snapshot$env_vars` must be an object.",
      field = "makevars_snapshot.env_vars"
    ))
  }
  if (!.bef_is_string(makevars$makevars_bytes_sha256)) {
    return(.bef_sidecar_result(
      FALSE,
      reason = "schema",
      message = "`makevars_snapshot$makevars_bytes_sha256` must be a string.",
      field = "makevars_snapshot.makevars_bytes_sha256"
    ))
  }
  if (!identical(makevars$makevars_bytes_sha256, "absent") &&
      !identical(makevars$makevars_bytes_sha256, "unavailable") &&
      !.bef_is_sha256_hex(makevars$makevars_bytes_sha256)) {
    return(.bef_sidecar_result(
      FALSE,
      reason = "schema",
      message = "`makevars_snapshot$makevars_bytes_sha256` must be absent, unavailable, or SHA-256.",
      field = "makevars_snapshot.makevars_bytes_sha256"
    ))
  }
  if (!is.list(makevars$r_cmd_config)) {
    return(.bef_sidecar_result(
      FALSE,
      reason = "schema",
      message = "`makevars_snapshot$r_cmd_config` must be an object.",
      field = "makevars_snapshot.r_cmd_config"
    ))
  }

  .bef_sidecar_result(TRUE)
}

.bef_sidecar_required_fields <- function() {
  c(
    "cache_format_version",
    "stan_file_sha256",
    "cmdstan_version",
    "cmdstanr_version",
    "arch",
    "compiler",
    "compiled_at",
    "cpp_options",
    "stanc_options",
    "makevars_snapshot",
    "os_major",
    "binary_sha256"
  )
}

.bef_cache_key_payload_fields <- function() {
  c(
    "stan_file_sha256",
    "cmdstan_version",
    "cmdstanr_version",
    "arch",
    "compiler",
    "cache_format_version",
    "model_family",
    "cpp_options_sha256",
    "stanc_options_sha256",
    "makevars_sha256",
    "os_major"
  )
}

.bef_sidecar_result <- function(valid, reason = NULL, meta = NULL, message = NULL, ...) {
  structure(
    c(
      list(
        valid = isTRUE(valid),
        reason = reason,
        message = message,
        meta = meta
      ),
      list(...)
    ),
    class = "bef_sidecar_validation"
  )
}

.bef_cache_key_hex <- function(key) {
  if (!.bef_is_string(key) || !.bef_is_sha256_hex(key)) {
    .bef_abort_sidecar_invalid(
      "`key` must be a 64-character lowercase SHA-256 hex string.",
      arg = "key",
      predicate = "64-character lowercase SHA-256 hex string"
    )
  }
  as.character(key)
}

.bef_is_sha256_hex <- function(x) {
  .bef_is_string(x) && grepl("^[0-9a-f]{64}$", x)
}

.bef_sidecar_string <- function(x, arg) {
  if (!.bef_is_string(x)) {
    .bef_abort_sidecar_invalid(
      sprintf("`%s` must be a non-empty character scalar.", arg),
      arg = arg,
      predicate = "non-empty character scalar"
    )
  }
  as.character(x)
}

.bef_sidecar_number <- function(x, arg, lower = -Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower) {
    .bef_abort_sidecar_invalid(
      sprintf("`%s` must be a finite numeric scalar >= %s.", arg, lower),
      arg = arg,
      predicate = sprintf("finite numeric scalar >= %s", lower)
    )
  }
  as.numeric(x)
}

.bef_sidecar_cpp_options <- function(x) {
  if (is.null(x)) {
    return(stats::setNames(list(), character()))
  }
  if (!is.list(x)) {
    .bef_abort_sidecar_invalid(
      "Sidecar `cpp_options` must be NULL or a list.",
      arg = "cpp_options",
      predicate = "NULL or list"
    )
  }
  x
}

.bef_sidecar_stanc_options <- function(x) {
  if (is.null(x)) {
    return(list())
  }
  if (!is.list(x)) {
    .bef_abort_sidecar_invalid(
      "Sidecar `stanc_options` must be NULL or a list.",
      arg = "stanc_options",
      predicate = "NULL or list"
    )
  }
  x
}

.bef_format_sidecar_time <- function(x) {
  if (.bef_is_string(x)) {
    return(as.character(x))
  }
  if (inherits(x, "POSIXt")) {
    return(format(x, "%Y-%m-%dT%H:%M:%S%z"))
  }
  .bef_abort_sidecar_invalid(
    "`compiled_at` must be a POSIX time or non-empty character scalar.",
    arg = "compiled_at",
    predicate = "POSIX time or non-empty character scalar"
  )
}

.bef_abort_sidecar_invalid <- function(message, ..., parent = NULL) {
  .bef_abort_invalid_args(
    message,
    ...,
    module = "cache-sidecar",
    stage = 6L,
    parent = parent
  )
}

.bef_cache_env <- function() {
  .bef_ensure_session_state()
  .bayesEfron_cache
}

.bef_session_cache_env <- function() {
  .bef_cache_env()
}

.bef_session_msgs_env <- function() {
  .bef_ensure_session_state()
  .bayesEfron_msgs_emitted
}

.bef_session_state <- function() {
  .bef_ensure_session_state()
  .bayesEfron_session
}

.bef_session_id <- function() {
  .bef_session_state()$session_id
}

.bef_cache_exists <- function(key) {
  exists(key, envir = .bef_cache_env(), inherits = FALSE)
}

.bef_cache_get <- function(key) {
  get(key, envir = .bef_cache_env(), inherits = FALSE)
}

.bef_cache_set <- function(key, value) {
  assign(key, value, envir = .bef_cache_env())
  invisible(value)
}

.bef_cache_remove <- function(key) {
  if (.bef_cache_exists(key)) {
    rm(list = key, envir = .bef_cache_env())
  }
  invisible(TRUE)
}

.bef_cache_keys <- function() {
  ls(envir = .bef_cache_env(), all.names = TRUE)
}

.bef_cache_clear_session <- function() {
  keys <- .bef_cache_keys()
  if (length(keys) > 0L) {
    rm(list = keys, envir = .bef_cache_env())
  }
  invisible(length(keys))
}

.bef_reset_session_cache <- function() {
  .bef_cache_clear_session()
  invisible(TRUE)
}

.bef_reset_session_messages <- function() {
  env <- .bef_session_msgs_env()
  keys <- ls(envir = env, all.names = TRUE)
  if (length(keys) > 0L) {
    rm(list = keys, envir = env)
  }
  invisible(TRUE)
}

.bef_reset_session_state <- function(cache = TRUE, messages = TRUE) {
  if (isTRUE(cache)) {
    .bef_reset_session_cache()
  }
  if (isTRUE(messages)) {
    .bef_reset_session_messages()
  }
  invisible(TRUE)
}
