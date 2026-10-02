# The compiled model: where it is stored, when it is compiled and when an
# existing file is reused. cmdstanr is replaced by a stand-in, so nothing is
# compiled here.

stan_ns <- function(name) getFromNamespace(name, "bayesEfron")

# A stand-in for cmdstanr::cmdstan_model(). When asked to compile it writes a
# small file at `exe_file`, as CmdStan would; on Unix the file is a script that
# exits with status 0, which is what the check of a compiled model looks for.
fake_cmdstan <- function(calls = new.env(parent = emptyenv()),
                         write_exe = TRUE,
                         compile_error = NULL) {
  calls$log <- list()
  function(stan_file, exe_file, compile = TRUE, ...) {
    calls$log[[length(calls$log) + 1L]] <- list(
      exe_file = exe_file, compile = compile, extra = list(...)
    )
    if (isTRUE(compile)) {
      if (!is.null(compile_error)) {
        stop(compile_error, call. = FALSE)
      }
      if (isTRUE(write_exe)) {
        writeLines(c("#!/bin/sh", "exit 0"), exe_file)
        Sys.chmod(exe_file, "0755")
      }
    }
    path <- exe_file
    structure(
      list(exe_file = function(new = NULL) {
        if (!is.null(new)) path <<- new
        path
      }),
      class = "fake_cmdstan_model"
    )
  }
}

local_fake_cmdstan <- function(calls = new.env(parent = emptyenv()),
                               ...,
                               env = parent.frame()) {
  withr::local_envvar(
    c(BAYESEFRON_CACHE_ROOT = withr::local_tempdir(.local_envir = env)),
    .local_envir = env
  )
  testthat::local_mocked_bindings(
    .bef_cmdstan_model = fake_cmdstan(calls, ...),
    .bef_cmdstan_version = function() "2.99.0",
    .bef_check_cmdstanr_installed = function() invisible(TRUE),
    .package = "bayesEfron",
    .env = env
  )
  models <- stan_ns(".bef_models")
  rm(list = ls(models, all.names = TRUE), envir = models)
  withr::defer(rm(list = ls(models, all.names = TRUE), envir = models), envir = env)
  calls
}

compiled_calls <- function(calls) {
  sum(vapply(calls$log, function(call) isTRUE(call$compile), logical(1)))
}

test_that("the executable name carries the Stan file hash, the CmdStan version and the architecture", {
  local_fake_cmdstan()
  stan_file <- stan_ns(".bef_stan_file")("RE")
  sha <- stan_ns(".bef_stan_file_sha256")(stan_file)
  exe <- stan_ns(".bef_exe_file")(stan_file)

  expect_match(sha, "^[0-9a-f]{64}$")
  expect_identical(sha, digest::digest(stan_file, algo = "sha256", file = TRUE))
  expect_identical(
    basename(exe),
    sprintf(
      "efron_re-%s-cmdstan2.99.0-%s%s",
      substr(sha, 1L, 12L),
      R.version$arch,
      if (identical(.Platform$OS.type, "windows")) ".exe" else ""
    )
  )
  # dirname() returns forward slashes on Windows, where the temporary
  # directory is written with backslashes.
  expect_identical(
    normalizePath(dirname(exe)),
    normalizePath(Sys.getenv("BAYESEFRON_CACHE_ROOT"))
  )
})

test_that("the cache directory is created, and defaults to the user cache directory", {
  root <- file.path(withr::local_tempdir(), "not", "there", "yet")
  withr::local_envvar(c(BAYESEFRON_CACHE_ROOT = root))
  expect_false(dir.exists(root))
  expect_identical(stan_ns(".bef_cache_dir")(create = FALSE), root)
  expect_false(dir.exists(root))
  expect_identical(stan_ns(".bef_cache_dir")(), root)
  expect_true(dir.exists(root))

  withr::local_envvar(c(BAYESEFRON_CACHE_ROOT = NA_character_))
  expect_identical(
    stan_ns(".bef_cache_dir")(create = FALSE),
    tools::R_user_dir("bayesEfron", which = "cache")
  )
})

test_that("the first call compiles, and the second call in a session reuses the model", {
  calls <- local_fake_cmdstan()

  first <- stan_ns(".bef_model")("RE")
  expect_length(calls$log, 1L)
  expect_true(calls$log[[1L]]$compile)
  expect_identical(calls$log[[1L]]$extra$stanc_options, list("O1"))
  expect_true(file.exists(first$exe_file()))
  expect_identical(first$exe_file(), stan_ns(".bef_exe_file")(stan_ns(".bef_stan_file")()))

  second <- stan_ns(".bef_model")("RE")
  expect_identical(second, first)
  expect_length(calls$log, 1L)
})

test_that("compilation uses a process-specific name and leaves only the final file", {
  calls <- local_fake_cmdstan()
  model <- stan_ns(".bef_model")("RE")

  compiled_as <- calls$log[[1L]]$exe_file
  expect_identical(
    compiled_as,
    sub("([.]exe)?$", sprintf("-%d\\1", Sys.getpid()), model$exe_file())
  )
  expect_identical(
    list.files(Sys.getenv("BAYESEFRON_CACHE_ROOT")),
    basename(model$exe_file())
  )
})

test_that("an existing executable is attached with compile = FALSE", {
  calls <- local_fake_cmdstan()
  stan_ns(".bef_model")("RE")
  models <- stan_ns(".bef_models")
  rm(list = ls(models, all.names = TRUE), envir = models)

  # A new session: nothing is attached yet, but the file is there.
  model <- stan_ns(".bef_model")("RE")
  expect_length(calls$log, 2L)
  expect_false(calls$log[[2L]]$compile)
  expect_identical(calls$log[[2L]]$exe_file, model$exe_file())
  expect_identical(compiled_calls(calls), 1L)
})

test_that("force_recompile compiles even when the model is attached", {
  calls <- local_fake_cmdstan()
  stan_ns(".bef_model")("RE")
  stan_ns(".bef_model")("RE", force_recompile = TRUE)
  expect_identical(compiled_calls(calls), 2L)
})

test_that("an empty executable is replaced by a fresh compile", {
  calls <- local_fake_cmdstan()
  exe <- stan_ns(".bef_exe_file")(stan_ns(".bef_stan_file")())
  file.create(exe)

  model <- stan_ns(".bef_model")("RE")
  expect_identical(compiled_calls(calls), 1L)
  expect_gt(file.size(model$exe_file()), 0)
})

test_that("an executable that does not run is replaced by a fresh compile", {
  # On Windows only an empty file is detected.
  skip_on_os("windows")
  calls <- local_fake_cmdstan()
  exe <- stan_ns(".bef_exe_file")(stan_ns(".bef_stan_file")())
  writeBin(as.raw(rep(c(1L, 200L, 37L, 99L), 256L)), exe)
  Sys.chmod(exe, "0755")
  expect_false(stan_ns(".bef_exe_ok")(exe))

  model <- stan_ns(".bef_model")("RE")
  expect_identical(compiled_calls(calls), 1L)
  expect_true(stan_ns(".bef_exe_ok")(model$exe_file()))
})

test_that("on Windows the names end in .exe, with the process id before the extension", {
  calls <- local_fake_cmdstan()
  testthat::local_mocked_bindings(
    .bef_is_windows = function() TRUE,
    .package = "bayesEfron"
  )

  model <- stan_ns(".bef_model")("RE")
  expect_match(basename(model$exe_file()), "[.]exe$")
  expect_match(
    basename(calls$log[[1L]]$exe_file),
    sprintf("-%d[.]exe$", Sys.getpid())
  )
})

test_that("a failed compile is a bef_compile_failed error", {
  local_fake_cmdstan(compile_error = "stanc says no")
  err <- tryCatch(stan_ns(".bef_model")("RE"), error = identity)
  expect_s3_class(err, "bef_compile_failed")
  expect_match(conditionMessage(err$parent), "stanc says no")

  # cmdstanr does not signal an error when it cannot write the executable.
  local_fake_cmdstan(write_exe = FALSE)
  err <- tryCatch(stan_ns(".bef_model")("RE"), error = identity)
  expect_s3_class(err, "bef_compile_failed")
  expect_match(conditionMessage(err), "could not be written")
})

test_that("a missing cmdstanr is reported before anything else happens", {
  testthat::local_mocked_bindings(
    .bef_check_cmdstanr_installed = function() {
      rlang::abort("cmdstanr is required", class = "rlib_error_package_not_found")
    },
    .package = "bayesEfron"
  )
  expect_error(stan_ns(".bef_model")("RE"), class = "rlib_error_package_not_found")
  expect_error(bayes_efron_compile(), class = "rlib_error_package_not_found")
})
