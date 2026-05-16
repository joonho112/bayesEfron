cache_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

cache_fixture_stan <- function(contents = "parameters { real y; }\nmodel { y ~ normal(0, 1); }\n") {
  path <- tempfile(fileext = ".stan")
  writeBin(charToRaw(contents), path)
  path
}

cache_fixture_makevars <- function(flag = "-O2") {
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
      CXX17FLAGS = flag,
      CXX17PICFLAGS = "-fPIC",
      CXX17STD = "-std=gnu++17"
    )
  )
}

cache_key_fixture <- function(stan_file = cache_fixture_stan(),
                              cmdstan_version = "2.34.1",
                              cmdstanr_version = "0.7.1",
                              arch = "aarch64",
                              compiler = "clang++",
                              model_family = "RE",
                              cpp_options = list(stan_threads = TRUE),
                              stanc_options = list("O1"),
                              makevars_snapshot = cache_fixture_makevars(),
                              os_major = "aarch64-apple-darwin22|Darwin-14") {
  cache_ns(".bef_cache_key")(
    stan_file = stan_file,
    cmdstan_version = cmdstan_version,
    cmdstanr_version = cmdstanr_version,
    arch = arch,
    compiler = compiler,
    model_family = model_family,
    cpp_options = cpp_options,
    stanc_options = stanc_options,
    makevars_snapshot = makevars_snapshot,
    os_major = os_major
  )
}

expect_cache_key_invalid <- function(err) {
  expect_s3_class(err, "bef_invalid_args")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
  expect_equal(err$module, "cache-key")
  expect_equal(err$stage, 6L)
}

test_that("cache key returns the canonical 11-field payload and fixed SHA-256", {
  out <- cache_key_fixture()

  expect_match(out$key, "^[0-9a-f]{64}$")
  expect_equal(
    out$key,
    "4eea46747a8cc1087f958244a72382cae670a569a9fd4e91b7351c138f9dde6b"
  )
  expect_named(
    out$payload,
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
  )
  expect_equal(out$payload[["stan_file_sha256"]],
               "1588bb9e523be95316d671c9035e1cf1bf80506aee2aa6235816d5c83baf4e58")
  expect_equal(out$payload[["cmdstan_version"]], "2.34.1")
  expect_equal(out$payload[["cache_format_version"]], "v1")
  expect_equal(out$payload[["model_family"]], "RE")
  expect_equal(out$cmdstan_v_pre, "2.34.1")
  expect_type(out$payload_json, "character")
  expect_type(out$provenance$cpp_options_normalized, "character")
})

test_that("cache key is deterministic for a fixed injected payload", {
  stan_file <- cache_fixture_stan()
  first <- cache_key_fixture(stan_file = stan_file)
  second <- cache_key_fixture(stan_file = stan_file)

  expect_identical(second$key, first$key)
  expect_identical(second$payload, first$payload)
  expect_identical(second$payload_json, first$payload_json)
})

test_that("Stan source bytes and injected platform fields affect the key", {
  base <- cache_key_fixture()

  variants <- list(
    cache_key_fixture(stan_file = cache_fixture_stan("parameters { real z; }\n")),
    cache_key_fixture(cmdstan_version = "2.35.0"),
    cache_key_fixture(cmdstanr_version = "0.8.0"),
    cache_key_fixture(arch = "x86_64"),
    cache_key_fixture(compiler = "g++"),
    cache_key_fixture(makevars_snapshot = cache_fixture_makevars(flag = "-O0")),
    cache_key_fixture(os_major = "x86_64-pc-linux-gnu|Linux-22")
  )

  expect_true(all(vapply(variants, function(x) !identical(x$key, base$key), logical(1))))
})

test_that("compile option normalization is deterministic and value-sensitive", {
  stan_file <- cache_fixture_stan()

  null_opts <- cache_key_fixture(stan_file = stan_file, cpp_options = NULL)
  empty_opts <- cache_key_fixture(stan_file = stan_file, cpp_options = list())
  expect_identical(null_opts$key, empty_opts$key)
  expect_identical(
    null_opts$provenance$cpp_options_normalized,
    empty_opts$provenance$cpp_options_normalized
  )

  threads_true <- cache_key_fixture(
    stan_file = stan_file,
    cpp_options = list(stan_threads = TRUE)
  )
  threads_false <- cache_key_fixture(
    stan_file = stan_file,
    cpp_options = list(stan_threads = FALSE)
  )
  expect_false(identical(threads_true$key, threads_false$key))

  named_a <- cache_key_fixture(
    stan_file = stan_file,
    cpp_options = list(stan_threads = TRUE, STAN_NO_RANGE_CHECKS = FALSE)
  )
  named_b <- cache_key_fixture(
    stan_file = stan_file,
    cpp_options = list(STAN_NO_RANGE_CHECKS = FALSE, stan_threads = TRUE)
  )
  expect_identical(named_a$key, named_b$key)
  expect_identical(
    named_a$provenance$cpp_options_normalized,
    named_b$provenance$cpp_options_normalized
  )
})

test_that("stanc option normalization distinguishes optimization levels", {
  stan_file <- cache_fixture_stan()

  o1 <- cache_key_fixture(stan_file = stan_file, stanc_options = list("O1"))
  o0 <- cache_key_fixture(stan_file = stan_file, stanc_options = list("O0"))
  expect_false(identical(o1$key, o0$key))

  ordered <- cache_key_fixture(
    stan_file = stan_file,
    stanc_options = list("--allow-undefined", "O1")
  )
  reversed <- cache_key_fixture(
    stan_file = stan_file,
    stanc_options = list("O1", "--allow-undefined")
  )
  expect_identical(ordered$key, reversed$key)
  expect_identical(
    ordered$provenance$stanc_options_normalized,
    reversed$provenance$stanc_options_normalized
  )
})

test_that("named option permutations fuzz to the same cache key", {
  stan_file <- cache_fixture_stan()
  opts <- list(
    stan_threads = TRUE,
    STAN_NO_RANGE_CHECKS = FALSE,
    alpha = 1,
    beta = "x",
    gamma = list(inner = TRUE)
  )
  target <- cache_key_fixture(stan_file = stan_file, cpp_options = opts)$key

  set.seed(20260511)
  keys <- vapply(
    seq_len(100L),
    function(i) {
      permuted <- opts[sample(names(opts))]
      cache_key_fixture(stan_file = stan_file, cpp_options = permuted)$key
    },
    character(1)
  )

  expect_true(all(keys == target))
})

test_that("malformed cache-key inputs fail with typed argument errors", {
  err <- tryCatch(cache_key_fixture(stan_file = tempfile()), error = identity)
  expect_cache_key_invalid(err)
  expect_equal(err$arg, "stan_file")

  err <- tryCatch(
    cache_key_fixture(cpp_options = list(stan_threads = TRUE, FALSE)),
    error = identity
  )
  expect_cache_key_invalid(err)
  expect_equal(err$arg, "options")

  err <- tryCatch(
    cache_key_fixture(model_family = "HE"),
    error = identity
  )
  expect_cache_key_invalid(err)
  expect_equal(err$arg, "model_family")

  err <- tryCatch(
    cache_key_fixture(makevars_snapshot = list(env_vars = list())),
    error = identity
  )
  expect_cache_key_invalid(err)
  expect_equal(err$arg, "makevars_snapshot")
})
