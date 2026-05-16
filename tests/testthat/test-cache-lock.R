lock_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

with_cache_lock_root <- function(code) {
  old <- Sys.getenv("BAYESEFRON_CACHE_ROOT", unset = NA_character_)
  root <- tempfile("bayesefron-lock-root-")
  Sys.setenv(BAYESEFRON_CACHE_ROOT = root)
  on.exit({
    if (is.na(old)) {
      Sys.unsetenv("BAYESEFRON_CACHE_ROOT")
    } else {
      Sys.setenv(BAYESEFRON_CACHE_ROOT = old)
    }
  }, add = TRUE)

  force(code)
}

lock_payload_for <- function(pid = Sys.getpid(),
                             nodename = Sys.info()[["nodename"]],
                             start_time = lock_ns(".bef_process_start_time")(pid),
                             acquired_at = as.POSIXct(
                               "2026-05-11 12:34:56",
                               tz = "UTC"
                             )) {
  lock_ns(".bef_lock_payload")(
    pid = pid,
    nodename = nodename,
    start_time = start_time,
    acquired_at = acquired_at
  )
}

write_lock_payload <- function(lock_path, payload) {
  dir.create(dirname(lock_path), recursive = TRUE, showWarnings = FALSE)
  lock_ns(".bef_write_lock_payload")(lock_path, payload)
}

expect_lock_timeout <- function(err) {
  expect_s3_class(err, "bef_lock_timeout")
  expect_s3_class(err, "bef_cache_error")
  expect_s3_class(err, "bef_error")
  expect_s3_class(err$parent, "bef_lock_owner")
}

expect_lock_invalid_args <- function(err) {
  expect_s3_class(err, "bef_invalid_args")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
  expect_equal(err$module, "cache-lock")
  expect_equal(err$stage, 6L)
}

test_that("acquire creates a four-line lock and release is idempotent", {
  with_cache_lock_root({
    acquire <- lock_ns(".bef_acquire_lock")
    release <- lock_ns(".bef_release_lock")
    read_payload <- lock_ns(".bef_read_lock_payload")
    lock_path <- lock_ns(".bef_cache_lock_path")()

    expect_false(file.exists(lock_path))
    expect_invisible(acquire(timeout = 0.1, poll_interval = 0.001))
    expect_true(file.exists(lock_path))
    expect_equal(basename(lock_path), ".lock")
    expect_equal(basename(dirname(lock_path)), "v1")

    lines <- readLines(lock_path, warn = FALSE)
    expect_length(lines, 4L)
    expect_true(all(nzchar(lines)))

    payload <- read_payload(lock_path)
    expect_false(isTRUE(payload$malformed))
    expect_equal(payload$pid, as.character(Sys.getpid()))
    expect_equal(payload$nodename, Sys.info()[["nodename"]])
    expect_match(payload$acquired_at, "^\\d{4}-\\d{2}-\\d{2}T")

    expect_invisible(release(lock_path))
    expect_false(file.exists(lock_path))
    expect_invisible(release(lock_path))
  })
})

test_that("an existing live same-process lock times out without truncation", {
  with_cache_lock_root({
    acquire <- lock_ns(".bef_acquire_lock")
    release <- lock_ns(".bef_release_lock")
    lock_path <- lock_ns(".bef_cache_lock_path")()

    acquire(timeout = 0.1, poll_interval = 0.001)
    before <- readLines(lock_path, warn = FALSE)
    err <- tryCatch(
      acquire(timeout = 0, poll_interval = 0),
      error = identity
    )
    after <- readLines(lock_path, warn = FALSE)
    expected_lock_path <- normalizePath(lock_path, winslash = "/", mustWork = TRUE)

    expect_lock_timeout(err)
    expect_identical(after, before)
    expect_equal(err$lock_path, expected_lock_path)
    expect_equal(err$timeout, 0)
    expect_equal(err$poll_interval, 0)
    expect_equal(err$owner_pid, as.character(Sys.getpid()))
    expect_equal(err$owner_host, Sys.info()[["nodename"]])
    expect_equal(err$parent$lock_path, expected_lock_path)
    expect_false(isTRUE(err$parent$owner$malformed))

    release(lock_path)
  })
})

test_that("dead same-host owner locks are recovered", {
  with_cache_lock_root({
    acquire <- lock_ns(".bef_acquire_lock")
    release <- lock_ns(".bef_release_lock")
    read_payload <- lock_ns(".bef_read_lock_payload")
    process_alive <- lock_ns(".bef_process_alive")
    lock_path <- lock_ns(".bef_cache_lock_path")()
    dead_pid <- 99999999L

    skip_if(process_alive(dead_pid), "sentinel PID unexpectedly exists")
    write_lock_payload(
      lock_path,
      lock_payload_for(
        pid = dead_pid,
        start_time = "ps_lstart:dead-owner"
      )
    )

    acquire(timeout = 0, poll_interval = 0)
    payload <- read_payload(lock_path)
    expect_equal(payload$pid, as.character(Sys.getpid()))
    expect_equal(payload$nodename, Sys.info()[["nodename"]])

    release(lock_path)
  })
})

test_that("Windows tasklist PID liveness parser is non-destructive and PID-specific", {
  tasklist_alive <- lock_ns(".bef_windows_tasklist_process_alive")

  live_task <- function(command, args, stdout, stderr) {
    expect_equal(command, "tasklist")
    expect_equal(args, c("/FI", "PID eq 1234", "/FO", "CSV", "/NH"))
    expect_true(isTRUE(stdout))
    expect_true(isTRUE(stderr))
    "\"Rterm.exe\",\"1234\",\"Console\",\"1\",\"50,000 K\""
  }
  expect_true(tasklist_alive(1234L, system2_fun = live_task))

  other_task <- function(command, args, stdout, stderr) {
    "\"Rterm.exe\",\"1234\",\"Console\",\"1\",\"50,000 K\""
  }
  expect_false(tasklist_alive(5678L, system2_fun = other_task))

  no_task <- function(command, args, stdout, stderr) {
    "INFO: No tasks are running which match the specified criteria."
  }
  expect_false(tasklist_alive(1234L, system2_fun = no_task))

  tasklist_error <- function(command, args, stdout, stderr) {
    out <- character()
    attr(out, "status") <- 1L
    out
  }
  expect_false(tasklist_alive(1234L, system2_fun = tasklist_error))
})

test_that("same-host locks with mismatched live process start time are recovered", {
  with_cache_lock_root({
    acquire <- lock_ns(".bef_acquire_lock")
    release <- lock_ns(".bef_release_lock")
    read_payload <- lock_ns(".bef_read_lock_payload")
    process_start_time <- lock_ns(".bef_process_start_time")
    start_time_known <- lock_ns(".bef_lock_start_time_known")
    lock_path <- lock_ns(".bef_cache_lock_path")()
    current_start <- process_start_time(Sys.getpid())

    skip_if_not(
      start_time_known(current_start),
      "process start time is unavailable on this platform"
    )
    write_lock_payload(
      lock_path,
      lock_payload_for(
        pid = Sys.getpid(),
        start_time = "ps_lstart:not-this-process"
      )
    )

    acquire(timeout = 0, poll_interval = 0)
    payload <- read_payload(lock_path)
    expect_equal(payload$pid, as.character(Sys.getpid()))
    expect_equal(payload$process_start_time, current_start)

    release(lock_path)
  })
})

test_that("foreign-host locks are never recovered locally", {
  with_cache_lock_root({
    acquire <- lock_ns(".bef_acquire_lock")
    lock_path <- lock_ns(".bef_cache_lock_path")()

    write_lock_payload(
      lock_path,
      lock_payload_for(
        pid = 99999999L,
        nodename = "foreign-host-for-bayesefron-tests",
        start_time = "ps_lstart:foreign"
      )
    )
    before <- readLines(lock_path, warn = FALSE)
    err <- tryCatch(
      acquire(timeout = 0, poll_interval = 0),
      error = identity
    )
    after <- readLines(lock_path, warn = FALSE)

    expect_lock_timeout(err)
    expect_identical(after, before)
    expect_equal(err$owner_host, "foreign-host-for-bayesefron-tests")
    unlink(lock_path, force = TRUE)
  })
})

test_that("malformed lock payloads are treated as held", {
  with_cache_lock_root({
    acquire <- lock_ns(".bef_acquire_lock")
    release <- lock_ns(".bef_release_lock")
    read_payload <- lock_ns(".bef_read_lock_payload")
    lock_path <- lock_ns(".bef_cache_lock_path")()

    dir.create(dirname(lock_path), recursive = TRUE, showWarnings = FALSE)
    writeLines(c(as.character(Sys.getpid()), Sys.info()[["nodename"]]), lock_path)
    before <- readLines(lock_path, warn = FALSE)
    malformed <- read_payload(lock_path)
    expect_true(isTRUE(malformed$malformed))

    err <- tryCatch(
      acquire(timeout = 0, poll_interval = 0),
      error = identity
    )
    after <- readLines(lock_path, warn = FALSE)

    expect_lock_timeout(err)
    expect_identical(after, before)
    expect_true(is.na(err$owner_pid))
    expect_invisible(release(lock_path))
    expect_true(file.exists(lock_path))
    unlink(lock_path, force = TRUE)
  })
})

test_that("release only removes locks owned by the current process", {
  with_cache_lock_root({
    release <- lock_ns(".bef_release_lock")
    lock_path <- lock_ns(".bef_cache_lock_path")()

    write_lock_payload(
      lock_path,
      lock_payload_for(
        pid = Sys.getpid(),
        nodename = "foreign-host-for-bayesefron-tests"
      )
    )
    expect_invisible(release(lock_path))
    expect_true(file.exists(lock_path))

    writeLines("partial", lock_path)
    expect_invisible(release(lock_path))
    expect_true(file.exists(lock_path))

    unlink(lock_path, force = TRUE)
  })
})

test_that("lock argument validation uses typed invalid-argument errors", {
  err <- tryCatch(
    lock_ns(".bef_acquire_lock")(timeout = -1),
    error = identity
  )
  expect_lock_invalid_args(err)
  expect_equal(err$arg, "timeout")

  err <- tryCatch(
    lock_ns(".bef_acquire_lock")(poll_interval = NA_real_),
    error = identity
  )
  expect_lock_invalid_args(err)
  expect_equal(err$arg, "poll_interval")
})
