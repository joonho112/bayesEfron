clear_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

with_clear_cache_root <- function(code) {
  old <- Sys.getenv("BAYESEFRON_CACHE_ROOT", unset = NA_character_)
  root <- tempfile("bayesefron-clear-root-")
  Sys.setenv(BAYESEFRON_CACHE_ROOT = root)
  clear_ns(".bef_cache_clear_session")()
  on.exit({
    clear_ns(".bef_cache_clear_session")()
    unlink(root, recursive = TRUE, force = TRUE)
    if (is.na(old)) {
      Sys.unsetenv("BAYESEFRON_CACHE_ROOT")
    } else {
      Sys.setenv(BAYESEFRON_CACHE_ROOT = old)
    }
  }, add = TRUE)

  force(code)
}

clear_counts <- function(lock_files = 0L,
                         session_keys = 0L,
                         disk_models = 0L,
                         disk_sidecars = 0L) {
  c(
    lock_files = as.integer(lock_files),
    session_keys = as.integer(session_keys),
    disk_models = as.integer(disk_models),
    disk_sidecars = as.integer(disk_sidecars)
  )
}

clear_key <- function(hex_digit) {
  paste(rep(hex_digit, 64L), collapse = "")
}

clear_write_entry <- function(key, cache_dir = clear_ns(".bef_cache_dir")(create = TRUE)) {
  paths <- list(
    exe = file.path(cache_dir, key),
    stan = file.path(cache_dir, paste0(key, ".stan")),
    meta = file.path(cache_dir, paste0(key, ".meta.json"))
  )
  writeBin(charToRaw("fake compiled executable\n"), paths$exe)
  writeBin(charToRaw("parameters { real y; }\n"), paths$stan)
  writeLines("{}", paths$meta, useBytes = TRUE)
  paths
}

test_that("clear-cache scope is validated with typed invalid-args errors", {
  with_clear_cache_root({
    err <- tryCatch(bayes_efron_clear_cache("bogus"), error = identity)

    expect_s3_class(err, "bef_invalid_args")
    expect_s3_class(err, "bef_pipeline_error")
    expect_s3_class(err, "bef_error")
    expect_equal(err$module, "clear-cache")
    expect_equal(err$stage, 6L)
    expect_equal(err$arg, "scope")
  })
})

test_that("default lock-only clearing removes only the cache lock", {
  with_clear_cache_root({
    cache_dir <- clear_ns(".bef_cache_dir")(create = TRUE)
    lock_path <- clear_ns(".bef_cache_lock_path")()
    paths <- clear_write_entry(clear_key("a"), cache_dir)
    clear_ns(".bef_cache_set")("tier1-key", list(model = "cached"))
    writeLines("malformed lock", lock_path, useBytes = TRUE)

    out <- NULL
    expect_message(
      out <- bayes_efron_clear_cache(),
      "lock_only"
    )

    expect_identical(out, clear_counts(lock_files = 1L))
    expect_false(file.exists(lock_path))
    expect_true(file.exists(paths$exe))
    expect_true(file.exists(paths$stan))
    expect_true(file.exists(paths$meta))
    expect_true(clear_ns(".bef_cache_exists")("tier1-key"))

    again <- suppressMessages(bayes_efron_clear_cache("lock_only"))
    expect_identical(again, clear_counts())
  })
})

test_that("lock-only clearing is a no-op when the cache directory is absent", {
  with_clear_cache_root({
    root <- clear_ns(".bef_cache_root")()

    out <- suppressMessages(bayes_efron_clear_cache("lock_only"))

    expect_identical(out, clear_counts())
    expect_false(file.exists(root))
  })
})

test_that("session clearing only drops Tier 1 models", {
  with_clear_cache_root({
    cache_dir <- clear_ns(".bef_cache_dir")(create = TRUE)
    lock_path <- clear_ns(".bef_cache_lock_path")()
    paths <- clear_write_entry(clear_key("b"), cache_dir)
    writeLines("malformed lock", lock_path, useBytes = TRUE)
    clear_ns(".bef_cache_set")("first", list(model = 1))
    clear_ns(".bef_cache_set")("second", list(model = 2))

    out <- suppressMessages(bayes_efron_clear_cache("session"))

    expect_identical(out, clear_counts(session_keys = 2L))
    expect_false(clear_ns(".bef_cache_exists")("first"))
    expect_false(clear_ns(".bef_cache_exists")("second"))
    expect_true(file.exists(lock_path))
    expect_true(file.exists(paths$exe))
    expect_true(file.exists(paths$stan))
    expect_true(file.exists(paths$meta))
  })
})

test_that("compiled-model clearing removes entry artifacts and preserves the lock", {
  with_clear_cache_root({
    cache_dir <- clear_ns(".bef_cache_dir")(create = TRUE)
    lock_path <- clear_ns(".bef_cache_lock_path")()
    paths_a <- clear_write_entry(clear_key("a"), cache_dir)
    paths_b <- clear_write_entry(clear_key("b"), cache_dir)
    orphan_source <- file.path(cache_dir, paste0(clear_key("c"), ".stan"))
    staging_dir <- file.path(cache_dir, ".staging-leftover")
    dir.create(staging_dir)
    writeBin(charToRaw("orphan source\n"), orphan_source)
    writeLines("malformed lock", lock_path, useBytes = TRUE)
    clear_ns(".bef_cache_set")("tier1-key", list(model = "cached"))

    out <- suppressMessages(bayes_efron_clear_cache("compiled_models"))

    expect_identical(out, clear_counts(disk_models = 2L, disk_sidecars = 2L))
    expect_true(dir.exists(cache_dir))
    expect_true(file.exists(lock_path))
    expect_true(dir.exists(staging_dir))
    expect_true(clear_ns(".bef_cache_exists")("tier1-key"))
    expect_false(file.exists(paths_a$exe))
    expect_false(file.exists(paths_a$stan))
    expect_false(file.exists(paths_a$meta))
    expect_false(file.exists(paths_b$exe))
    expect_false(file.exists(paths_b$stan))
    expect_false(file.exists(paths_b$meta))
    expect_false(file.exists(orphan_source))

    again <- suppressMessages(bayes_efron_clear_cache("compiled_models"))
    expect_identical(again, clear_counts())
  })
})

test_that("all clearing removes session state and the full cache root", {
  with_clear_cache_root({
    root <- clear_ns(".bef_cache_root")()
    cache_dir <- clear_ns(".bef_cache_dir")(create = TRUE)
    lock_path <- clear_ns(".bef_cache_lock_path")()
    paths <- clear_write_entry(clear_key("d"), cache_dir)
    writeLines("malformed lock", lock_path, useBytes = TRUE)
    writeLines("residual", file.path(root, "residual.txt"), useBytes = TRUE)
    dir.create(file.path(cache_dir, ".staging-leftover"))
    clear_ns(".bef_cache_set")("first", list(model = 1))
    clear_ns(".bef_cache_set")("second", list(model = 2))

    out <- suppressMessages(bayes_efron_clear_cache("all"))

    expect_identical(
      out,
      clear_counts(
        lock_files = 1L,
        session_keys = 2L,
        disk_models = 1L,
        disk_sidecars = 1L
      )
    )
    expect_false(file.exists(root))
    expect_false(file.exists(paths$exe))
    expect_false(clear_ns(".bef_cache_exists")("first"))
    expect_false(clear_ns(".bef_cache_exists")("second"))

    again <- suppressMessages(bayes_efron_clear_cache("all"))
    expect_identical(again, clear_counts())
    expect_false(file.exists(root))
  })
})
