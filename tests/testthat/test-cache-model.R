model_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

with_cache_model_root <- function(code) {
  old <- Sys.getenv("BAYESEFRON_CACHE_ROOT", unset = NA_character_)
  root <- tempfile("bayesefron-model-root-")
  Sys.setenv(BAYESEFRON_CACHE_ROOT = root)
  model_ns(".bef_cache_clear_session")()
  on.exit({
    model_ns(".bef_cache_clear_session")()
    if (is.na(old)) {
      Sys.unsetenv("BAYESEFRON_CACHE_ROOT")
    } else {
      Sys.setenv(BAYESEFRON_CACHE_ROOT = old)
    }
  }, add = TRUE)

  force(code)
}

model_fixture_stan <- function(contents = "parameters { real y; }\nmodel { y ~ normal(0, 1); }\n") {
  path <- tempfile(fileext = ".stan")
  writeBin(charToRaw(contents), path)
  path
}

model_fixture_makevars <- function() {
  list(
    env_vars = list(
      CXXFLAGS = "",
      CFLAGS = "",
      LDFLAGS = "",
      PKG_CXXFLAGS = "",
      PKG_CPPFLAGS = ""
    ),
    makevars_bytes_sha256 = "absent",
    r_cmd_config = list(
      CXX17 = "clang++",
      CXX17FLAGS = "-O2",
      CXX17PICFLAGS = "-fPIC",
      CXX17STD = "-std=gnu++17"
    )
  )
}

model_key_fixture <- function(stan_file = model_fixture_stan(),
                              cpp_options = list(stan_threads = TRUE),
                              stanc_options = list("O1"),
                              cmdstan_version = "2.34.1") {
  model_ns(".bef_cache_key")(
    stan_file = stan_file,
    cmdstan_version = cmdstan_version,
    cmdstanr_version = "0.7.1",
    arch = "aarch64",
    compiler = "clang++",
    cpp_options = cpp_options,
    stanc_options = stanc_options,
    makevars_snapshot = model_fixture_makevars(),
    os_major = "aarch64-apple-darwin22|Darwin-14"
  )
}

new_fake_cmdstan <- function(fail_compile = FALSE,
                             fail_reattach = FALSE,
                             compile_exe_basename = NULL) {
  env <- new.env(parent = emptyenv())
  env$calls <- list()

  factory <- function(stan_file,
                      compile,
                      dir = NULL,
                      cpp_options = NULL,
                      stanc_options = NULL,
                      exe_file = NULL,
                      ...) {
    call <- list(
      stan_file = stan_file,
      compile = compile,
      dir = dir,
      cpp_options = cpp_options,
      stanc_options = stanc_options,
      exe_file = exe_file
    )
    env$calls[[length(env$calls) + 1L]] <- call

    if (isTRUE(compile)) {
      if (isTRUE(fail_compile)) {
        stop("fake compile failure", call. = FALSE)
      }
      exe_name <- if (is.null(compile_exe_basename)) {
        sub("[.]stan$", "", basename(stan_file))
      } else {
        compile_exe_basename
      }
      compiled_exe <- file.path(dir, exe_name)
      writeBin(charToRaw("fake compiled executable\n"), compiled_exe)
      return(structure(
        list(
          kind = "compiled",
          stan_file = stan_file,
          exe_path = compiled_exe,
          exe_file = function() compiled_exe
        ),
        class = "fake_cmdstan_model"
      ))
    }

    if (isTRUE(fail_reattach)) {
      stop("fake reattach failure", call. = FALSE)
    }
    if (!file.exists(stan_file) || !file.exists(exe_file)) {
      stop("fake reattach missing artifact", call. = FALSE)
    }
    structure(
      list(
        kind = "reattached",
        stan_file = stan_file,
        exe_path = exe_file,
        exe_file = function() exe_file
      ),
      class = "fake_cmdstan_model"
    )
  }

  list(
    factory = factory,
    calls = function() env$calls
  )
}

write_model_cache_entry <- function(key_info, paths, exe_bytes = "cached executable\n") {
  file.copy(key_info$provenance$stan_file, paths$stan)
  writeBin(charToRaw(exe_bytes), paths$exe)
  meta <- model_ns(".bef_build_sidecar")(
    key_info,
    paths$exe,
    compiled_at = as.POSIXct("2026-05-11 12:34:56", tz = "UTC"),
    cmdstan_v_post = key_info$cmdstan_v_pre,
    compile_seconds = 1
  )
  model_ns(".bef_write_sidecar")(meta, paths$meta)
  invisible(meta)
}

call_bef_model <- function(key_info,
                           fake = new_fake_cmdstan(),
                           version = function() key_info$cmdstan_v_pre,
                           force_recompile = FALSE,
                           acquire_lock_fun = model_ns(".bef_acquire_lock"),
                           release_lock_fun = model_ns(".bef_release_lock")) {
  model_ns(".bef_model")(
    "RE",
    cpp_options = key_info$provenance$cpp_options,
    stanc_options = key_info$provenance$stanc_options,
    force_recompile = force_recompile,
    stan_file = key_info$provenance$stan_file,
    key_info = key_info,
    cmdstan_model_fun = fake$factory,
    cmdstan_version_fun = version,
    acquire_lock_fun = acquire_lock_fun,
    release_lock_fun = release_lock_fun,
    check_installed = FALSE
  )
}

expect_model_cache_invalid_args <- function(err) {
  expect_s3_class(err, "bef_invalid_args")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
  expect_equal(err$module, "model-cache")
  expect_equal(err$stage, 6L)
}

test_that("Tier 1 hit returns without locking or disk access", {
  with_cache_model_root({
    key_info <- model_key_fixture()
    fake <- new_fake_cmdstan()
    cached <- structure(list(kind = "tier1"), class = "fake_cmdstan_model")
    model_ns(".bef_cache_set")(key_info$key, cached)

    acquire <- function(...) stop("lock should not be acquired", call. = FALSE)
    out <- call_bef_model(key_info, fake = fake, acquire_lock_fun = acquire)

    expect_identical(out, cached)
    expect_length(fake$calls(), 0L)
  })
})

test_that("Tier 1 is rechecked after lock acquisition", {
  with_cache_model_root({
    key_info <- model_key_fixture()
    fake <- new_fake_cmdstan()
    cached <- structure(list(kind = "tier1-after-lock"), class = "fake_cmdstan_model")
    acquired <- 0L
    released <- 0L
    acquire <- function(timeout, lock_path) {
      acquired <<- acquired + 1L
      model_ns(".bef_cache_set")(key_info$key, cached)
      invisible(lock_path)
    }
    release <- function(lock_path) {
      released <<- released + 1L
      invisible(TRUE)
    }

    out <- call_bef_model(
      key_info,
      fake = fake,
      acquire_lock_fun = acquire,
      release_lock_fun = release
    )

    expect_identical(out, cached)
    expect_equal(acquired, 1L)
    expect_equal(released, 1L)
    expect_length(fake$calls(), 0L)
  })
})

test_that("valid Tier 2 sidecar reattaches and memoizes without compiling", {
  with_cache_model_root({
    key_info <- model_key_fixture()
    paths <- model_ns(".bef_cache_entry_paths")(key_info$key)
    write_model_cache_entry(key_info, paths)
    fake <- new_fake_cmdstan()

    out <- call_bef_model(key_info, fake = fake)
    calls <- fake$calls()

    expect_equal(out$kind, "reattached")
    expect_length(calls, 1L)
    expect_false(isTRUE(calls[[1L]]$compile))
    expect_equal(calls[[1L]]$stan_file, paths$stan)
    expect_equal(calls[[1L]]$exe_file, paths$exe)
    expect_identical(model_ns(".bef_cache_get")(key_info$key), out)
  })
})

test_that("missing Tier 2 entry cold compiles, promotes, reattaches, and memoizes", {
  with_cache_model_root({
    key_info <- model_key_fixture()
    paths <- model_ns(".bef_cache_entry_paths")(key_info$key)
    fake <- new_fake_cmdstan()

    out <- call_bef_model(key_info, fake = fake)
    calls <- fake$calls()

    expect_equal(out$kind, "reattached")
    expect_length(calls, 2L)
    expect_true(isTRUE(calls[[1L]]$compile))
    expect_false(isTRUE(calls[[2L]]$compile))
    expect_equal(calls[[1L]]$cpp_options, list(stan_threads = TRUE))
    expect_equal(calls[[1L]]$stanc_options, list("O1"))
    expect_true(all(file.exists(c(paths$stan, paths$exe, paths$meta))))
    expect_false(any(grepl("^\\.staging-", list.files(paths$cache_dir))))

    read <- model_ns(".bef_read_sidecar")(paths$meta)
    expect_true(read$valid)
    expect_equal(read$meta$cmdstan_version, "2.34.1")
    expect_equal(read$meta$cpp_options, list(stan_threads = TRUE))
    expect_identical(model_ns(".bef_cache_get")(key_info$key), out)
  })
})

test_that("compile output can be normalized from returned model exe_file", {
  with_cache_model_root({
    key_info <- model_key_fixture()
    paths <- model_ns(".bef_cache_entry_paths")(key_info$key)
    fake <- new_fake_cmdstan(compile_exe_basename = "cmdstan-output")

    out <- call_bef_model(key_info, fake = fake)

    expect_equal(out$kind, "reattached")
    expect_true(file.exists(paths$exe))
    expect_equal(readLines(paths$exe, warn = FALSE), "fake compiled executable")
  })
})

test_that("invalid sidecar parse deletes stale artifacts and recompiles", {
  with_cache_model_root({
    key_info <- model_key_fixture()
    paths <- model_ns(".bef_cache_entry_paths")(key_info$key)
    dir.create(dirname(paths$meta), recursive = TRUE, showWarnings = FALSE)
    writeLines("{bad-json", paths$meta)
    writeLines("old stan", paths$stan)
    writeLines("old exe", paths$exe)
    fake <- new_fake_cmdstan()

    out <- call_bef_model(key_info, fake = fake)

    expect_equal(out$kind, "reattached")
    expect_length(fake$calls(), 2L)
    expect_true(model_ns(".bef_read_sidecar")(paths$meta)$valid)
    expect_equal(readLines(paths$stan, warn = FALSE), readLines(key_info$provenance$stan_file, warn = FALSE))
  })
})

test_that("binary checksum mismatch warns, deletes, and recompiles", {
  with_cache_model_root({
    key_info <- model_key_fixture()
    paths <- model_ns(".bef_cache_entry_paths")(key_info$key)
    write_model_cache_entry(key_info, paths)
    writeBin(charToRaw("corrupt executable\n"), paths$exe)
    fake <- new_fake_cmdstan()
    warning <- NULL

    out <- withCallingHandlers(
      call_bef_model(key_info, fake = fake),
      bef_cache_corruption = function(w) {
        warning <<- w
        invokeRestart("muffleWarning")
      }
    )

    expect_equal(out$kind, "reattached")
    expect_s3_class(warning, "bef_cache_corruption")
    expect_length(fake$calls(), 2L)
    expect_true(model_ns(".bef_read_sidecar")(paths$meta)$valid)
  })
})

test_that("reattach failure deletes entry and recompiles", {
  with_cache_model_root({
    key_info <- model_key_fixture()
    paths <- model_ns(".bef_cache_entry_paths")(key_info$key)
    write_model_cache_entry(key_info, paths)
    fail_once <- new.env(parent = emptyenv())
    fail_once$n <- 0L
    fake <- new_fake_cmdstan()
    factory <- function(...) {
      args <- list(...)
      if (!isTRUE(args$compile)) {
        fail_once$n <- fail_once$n + 1L
        if (identical(fail_once$n, 1L)) {
          stop("first reattach fails", call. = FALSE)
        }
      }
      fake$factory(...)
    }

    out <- model_ns(".bef_model")(
      "RE",
      cpp_options = key_info$provenance$cpp_options,
      stanc_options = key_info$provenance$stanc_options,
      stan_file = key_info$provenance$stan_file,
      key_info = key_info,
      cmdstan_model_fun = factory,
      cmdstan_version_fun = function() key_info$cmdstan_v_pre,
      check_installed = FALSE
    )

    expect_equal(out$kind, "reattached")
    expect_equal(fail_once$n, 2L)
    expect_length(fake$calls(), 2L)
  })
})

test_that("compile failures release lock, remove staging, and leave no entry", {
  with_cache_model_root({
    key_info <- model_key_fixture()
    paths <- model_ns(".bef_cache_entry_paths")(key_info$key)
    fake <- new_fake_cmdstan(fail_compile = TRUE)

    err <- tryCatch(call_bef_model(key_info, fake = fake), error = identity)

    expect_s3_class(err, "bef_compile_failed")
    expect_s3_class(err, "bef_pipeline_error")
    expect_false(file.exists(model_ns(".bef_cache_lock_path")()))
    expect_false(any(file.exists(c(paths$stan, paths$exe, paths$meta))))
    expect_false(any(grepl("^\\.staging-", list.files(paths$cache_dir))))
  })
})

test_that("CmdStan version drift aborts without final artifacts", {
  with_cache_model_root({
    key_info <- model_key_fixture(cmdstan_version = "2.34.1")
    paths <- model_ns(".bef_cache_entry_paths")(key_info$key)
    fake <- new_fake_cmdstan()

    err <- tryCatch(
      call_bef_model(key_info, fake = fake, version = function() "2.35.0"),
      error = identity
    )

    expect_s3_class(err, "bef_cache_version_drift")
    expect_s3_class(err, "bef_cache_error")
    expect_equal(err$cmdstan_v_pre, "2.34.1")
    expect_equal(err$cmdstan_v_post, "2.35.0")
    expect_false(file.exists(model_ns(".bef_cache_lock_path")()))
    expect_false(any(file.exists(c(paths$stan, paths$exe, paths$meta))))
    expect_false(any(grepl("^\\.staging-", list.files(paths$cache_dir))))
  })
})

test_that("force_recompile bypasses Tier 1 and Tier 2 entries", {
  with_cache_model_root({
    key_info <- model_key_fixture()
    paths <- model_ns(".bef_cache_entry_paths")(key_info$key)
    write_model_cache_entry(key_info, paths, exe_bytes = "old executable\n")
    model_ns(".bef_cache_set")(key_info$key, list(kind = "tier1"))
    fake <- new_fake_cmdstan()

    out <- call_bef_model(key_info, fake = fake, force_recompile = TRUE)

    expect_equal(out$kind, "reattached")
    expect_length(fake$calls(), 2L)
    expect_equal(readLines(paths$exe, warn = FALSE), "fake compiled executable")
    expect_identical(model_ns(".bef_cache_get")(key_info$key), out)
  })
})

test_that("format mismatch self-heals by deleting and recompiling", {
  with_cache_model_root({
    key_info <- model_key_fixture()
    paths <- model_ns(".bef_cache_entry_paths")(key_info$key)
    meta <- write_model_cache_entry(key_info, paths)
    meta$cache_format_version <- "v2"
    model_ns(".bef_write_sidecar")(meta, paths$meta)
    fake <- new_fake_cmdstan()

    out <- call_bef_model(key_info, fake = fake)

    expect_equal(out$kind, "reattached")
    expect_length(fake$calls(), 2L)
    read <- model_ns(".bef_read_sidecar")(paths$meta)
    expect_true(read$valid)
    expect_equal(read$meta$cache_format_version, "v1")
  })
})

test_that("invalid model-cache arguments are typed", {
  err <- tryCatch(
    model_ns(".bef_model")("HE", check_installed = FALSE),
    error = identity
  )
  expect_model_cache_invalid_args(err)
  expect_equal(err$arg, "model_name")

  err <- tryCatch(
    model_ns(".bef_model")("RE", force_recompile = NA, check_installed = FALSE),
    error = identity
  )
  expect_model_cache_invalid_args(err)
  expect_equal(err$arg, "force_recompile")
})
