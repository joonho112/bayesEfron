sidecar_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

with_cache_sidecar_root <- function(code) {
  old <- Sys.getenv("BAYESEFRON_CACHE_ROOT", unset = NA_character_)
  root <- tempfile("bayesefron-sidecar-root-")
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

sidecar_fixture_stan <- function(contents = "parameters { real y; }\nmodel { y ~ normal(0, 1); }\n") {
  path <- tempfile(fileext = ".stan")
  writeBin(charToRaw(contents), path)
  path
}

sidecar_fixture_makevars <- function() {
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

sidecar_key_fixture <- function(stan_file = sidecar_fixture_stan(),
                                cpp_options = list(stan_threads = TRUE),
                                stanc_options = list("O1")) {
  sidecar_ns(".bef_cache_key")(
    stan_file = stan_file,
    cmdstan_version = "2.34.1",
    cmdstanr_version = "0.7.1",
    arch = "aarch64",
    compiler = "clang++",
    cpp_options = cpp_options,
    stanc_options = stanc_options,
    makevars_snapshot = sidecar_fixture_makevars(),
    os_major = "aarch64-apple-darwin22|Darwin-14"
  )
}

write_sidecar_fixture <- function(key_info = sidecar_key_fixture()) {
  paths <- sidecar_ns(".bef_cache_entry_paths")(key_info$key)
  file.copy(key_info$provenance$stan_file, paths$stan)
  writeBin(charToRaw("fake cmdstan executable bytes\n"), paths$exe)
  meta <- sidecar_ns(".bef_build_sidecar")(
    key_info,
    paths$exe,
    compiled_at = as.POSIXct("2026-05-11 12:34:56", tz = "UTC"),
    cmdstan_v_post = "2.34.1",
    compile_seconds = 207,
    host = "alabaster.local",
    r_version = "R version fixture"
  )
  sidecar_ns(".bef_write_sidecar")(meta, paths$meta)
  list(key_info = key_info, paths = paths, meta = meta)
}

expect_sidecar_invalid_args <- function(err) {
  expect_s3_class(err, "bef_invalid_args")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
  expect_equal(err$module, "cache-sidecar")
  expect_equal(err$stage, 6L)
}

expect_cache_format_mismatch <- function(err) {
  expect_s3_class(err, "bef_cache_format_mismatch")
  expect_s3_class(err, "bef_cache_error")
  expect_s3_class(err, "bef_error")
}

test_that("cache entry paths and staging layout follow the v1 contract", {
  with_cache_sidecar_root({
    key <- sidecar_key_fixture()$key
    paths <- sidecar_ns(".bef_cache_entry_paths")(key)

    expect_equal(basename(paths$cache_dir), "v1")
    expect_equal(basename(paths$stan), paste0(key, ".stan"))
    expect_equal(basename(paths$exe), key)
    expect_equal(basename(paths$meta), paste0(key, ".meta.json"))
    expect_equal(sidecar_ns(".bef_sidecar_path")(key), paths$meta)

    staging_dir <- sidecar_ns(".bef_staging_dir")()
    expect_true(dir.exists(staging_dir))
    expect_equal(dirname(staging_dir), paths$cache_dir)
    expect_match(basename(staging_dir), "^\\.staging-")
    if (!identical(.Platform$OS.type, "windows")) {
      expect_equal(as.character(file.info(staging_dir)$mode), "700")
    }

    stage_paths <- sidecar_ns(".bef_stage_entry_paths")(key, staging_dir)
    expect_equal(basename(stage_paths$stan), basename(paths$stan))
    expect_equal(basename(stage_paths$exe), basename(paths$exe))
    expect_equal(basename(stage_paths$meta), basename(paths$meta))

    err <- tryCatch(sidecar_ns(".bef_cache_entry_paths")("not-a-key"), error = identity)
    expect_sidecar_invalid_args(err)
    expect_equal(err$arg, "key")
  })
})

test_that("sidecar builder emits required and optional provenance fields", {
  with_cache_sidecar_root({
    fixture <- write_sidecar_fixture()
    meta <- fixture$meta

    expect_named(
      meta,
      c(
        sidecar_ns(".bef_sidecar_required_fields")(),
        "host",
        "r_version",
        "compile_seconds"
      )
    )
    expect_equal(meta$cache_format_version, "v1")
    expect_equal(meta$cmdstan_version, "2.34.1")
    expect_equal(meta$cpp_options, list(stan_threads = TRUE))
    expect_equal(meta$stanc_options, list("O1"))
    expect_equal(meta$makevars_snapshot, sidecar_fixture_makevars())
    expect_equal(meta$os_major, "aarch64-apple-darwin22|Darwin-14")
    expect_equal(
      meta$binary_sha256,
      digest::digest(fixture$paths$exe, algo = "sha256", file = TRUE)
    )

    minimal <- sidecar_ns(".bef_build_sidecar")(
      fixture$key_info,
      fixture$paths$exe,
      compiled_at = "2026-05-11T12:34:56+0000",
      cmdstan_v_post = "2.34.1",
      host = NULL,
      r_version = NULL
    )
    expect_false("host" %in% names(minimal))
    expect_false("r_version" %in% names(minimal))
    expect_false("compile_seconds" %in% names(minimal))
    expect_named(minimal, sidecar_ns(".bef_sidecar_required_fields")())
  })
})

test_that("sidecar JSON round-trips and validates against files", {
  with_cache_sidecar_root({
    fixture <- write_sidecar_fixture()
    read <- sidecar_ns(".bef_read_sidecar")(fixture$paths$meta)

    expect_s3_class(read, "bef_sidecar_validation")
    expect_true(read$valid)
    expect_equal(read$meta$cache_format_version, "v1")
    expect_equal(read$meta$cpp_options, list(stan_threads = TRUE))
    expect_equal(read$meta$stanc_options, list("O1"))
    expect_equal(read$meta$compile_seconds, 207)

    valid <- sidecar_ns(".bef_validate_sidecar")(
      read$meta,
      key_info = fixture$key_info,
      paths = fixture$paths
    )
    expect_true(valid$valid)
    expect_equal(valid$paths$key, fixture$key_info$key)
  })
})

test_that("sidecar read reports parse and missing-file failures as invalid results", {
  with_cache_sidecar_root({
    paths <- sidecar_ns(".bef_cache_entry_paths")(sidecar_key_fixture()$key)

    missing <- sidecar_ns(".bef_read_sidecar")(paths$meta)
    expect_false(missing$valid)
    expect_equal(missing$reason, "missing")

    dir.create(dirname(paths$meta), recursive = TRUE, showWarnings = FALSE)
    writeLines("{not-json", paths$meta)
    parsed <- sidecar_ns(".bef_read_sidecar")(paths$meta)
    expect_false(parsed$valid)
    expect_equal(parsed$reason, "parse")
  })
})

test_that("sidecar schema rejects missing and malformed required fields", {
  with_cache_sidecar_root({
    fixture <- write_sidecar_fixture()
    validate <- sidecar_ns(".bef_validate_sidecar")

    for (field in sidecar_ns(".bef_sidecar_required_fields")()) {
      bad <- fixture$meta
      bad[[field]] <- NULL
      result <- validate(bad)
      expect_false(result$valid, info = field)
      expect_equal(result$reason, "schema", info = field)
      expect_true(field %in% result$missing_fields, info = field)
    }

    bad <- fixture$meta
    bad$binary_sha256 <- "abc"
    result <- validate(bad)
    expect_false(result$valid)
    expect_equal(result$field, "binary_sha256")

    bad <- fixture$meta
    bad$makevars_snapshot$makevars_bytes_sha256 <- "not-a-digest"
    result <- validate(bad)
    expect_false(result$valid)
    expect_equal(result$field, "makevars_snapshot.makevars_bytes_sha256")

    bad <- fixture$meta
    bad$compile_seconds <- -1
    result <- validate(bad)
    expect_false(result$valid)
    expect_equal(result$field, "compile_seconds")
  })
})

test_that("cache format mismatches raise the typed firewall condition", {
  with_cache_sidecar_root({
    fixture <- write_sidecar_fixture()
    bad <- fixture$meta
    bad$cache_format_version <- "v2"

    err <- tryCatch(
      sidecar_ns(".bef_validate_sidecar")(bad),
      error = identity
    )
    expect_cache_format_mismatch(err)
    expect_equal(err$cache_format_version, "v2")
    expect_equal(err$expected_cache_format_version, "v1")
  })
})

test_that("sidecar identity detects payload and cached source drift", {
  with_cache_sidecar_root({
    fixture <- write_sidecar_fixture()

    bad <- fixture$meta
    bad$cmdstanr_version <- "0.8.0"
    result <- sidecar_ns(".bef_validate_sidecar")(
      bad,
      key_info = fixture$key_info
    )
    expect_false(result$valid)
    expect_equal(result$reason, "payload_mismatch")
    expect_equal(result$field, "cmdstanr_version")

    writeBin(charToRaw("modified stan source\n"), fixture$paths$stan)
    result <- sidecar_ns(".bef_validate_sidecar")(
      fixture$meta,
      key_info = fixture$key_info,
      paths = fixture$paths
    )
    expect_false(result$valid)
    expect_equal(result$reason, "stan_file_sha256_mismatch")
  })
})

test_that("binary checksum mismatch warns with cache-corruption metadata", {
  with_cache_sidecar_root({
    fixture <- write_sidecar_fixture()
    writeBin(charToRaw("altered executable bytes\n"), fixture$paths$exe)

    warning <- NULL
    result <- withCallingHandlers(
      sidecar_ns(".bef_validate_sidecar")(
        fixture$meta,
        key_info = fixture$key_info,
        paths = fixture$paths
      ),
      bef_cache_corruption = function(w) {
        warning <<- w
        invokeRestart("muffleWarning")
      }
    )

    expect_false(result$valid)
    expect_equal(result$reason, "binary_sha256_mismatch")
    expect_s3_class(warning, "bef_cache_corruption")
    expect_s3_class(warning, "bef_cache_warning")
    expect_equal(warning$recorded_sha256, fixture$meta$binary_sha256)
    expect_equal(warning$key, fixture$key_info$key)
    expect_equal(warning$cache_dir, fixture$paths$cache_dir)
    expect_match(warning$computed_sha256, "^[0-9a-f]{64}$")
  })
})

test_that("binary digest hashes file contents and deletion removes entry artifacts", {
  with_cache_sidecar_root({
    fixture <- write_sidecar_fixture()

    expect_equal(
      sidecar_ns(".bef_binary_sha256")(fixture$paths$exe),
      digest::digest(fixture$paths$exe, algo = "sha256", file = TRUE)
    )
    expect_false(identical(
      sidecar_ns(".bef_binary_sha256")(fixture$paths$exe),
      digest::digest(fixture$paths$exe, algo = "sha256", serialize = FALSE)
    ))

    expect_true(all(file.exists(c(
      fixture$paths$stan,
      fixture$paths$exe,
      fixture$paths$meta
    ))))
    removed <- sidecar_ns(".bef_delete_cache_entry")(fixture$paths)
    expect_equal(removed, c(stan = 1L, exe = 1L, meta = 1L))
    expect_false(any(file.exists(c(
      fixture$paths$stan,
      fixture$paths$exe,
      fixture$paths$meta
    ))))
  })
})

test_that("staged artifacts promote in meta, source, executable order", {
  with_cache_sidecar_root({
    key <- sidecar_key_fixture()$key
    entry_paths <- sidecar_ns(".bef_cache_entry_paths")(key)
    staging_dir <- sidecar_ns(".bef_staging_dir")()
    stage_paths <- sidecar_ns(".bef_stage_entry_paths")(key, staging_dir)

    writeLines("meta", stage_paths$meta)
    writeLines("stan", stage_paths$stan)
    writeLines("exe", stage_paths$exe)

    promoted <- sidecar_ns(".bef_promote_staged_entry")(stage_paths, entry_paths)
    expect_identical(promoted, entry_paths)
    expect_equal(readLines(entry_paths$meta, warn = FALSE), "meta")
    expect_equal(readLines(entry_paths$stan, warn = FALSE), "stan")
    expect_equal(readLines(entry_paths$exe, warn = FALSE), "exe")
    expect_false(file.exists(stage_paths$meta))
    expect_false(file.exists(stage_paths$stan))
    expect_false(file.exists(stage_paths$exe))
  })
})
