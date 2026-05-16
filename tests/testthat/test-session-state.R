session_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

test_that("package session environments are initialized on load", {
  cache_env <- session_ns(".bef_cache_env")()
  msgs_env <- session_ns(".bef_session_msgs_env")()
  state_env <- session_ns(".bef_session_state")()

  expect_true(is.environment(cache_env))
  expect_true(is.environment(msgs_env))
  expect_true(is.environment(state_env))
  expect_identical(parent.env(cache_env), emptyenv())
  expect_identical(parent.env(msgs_env), emptyenv())
  expect_identical(parent.env(state_env), emptyenv())

  session_id <- session_ns(".bef_session_id")()
  expect_type(session_id, "character")
  expect_length(session_id, 1L)
  expect_true(nzchar(session_id))
  expect_identical(session_ns(".bef_session_id")(), session_id)

  expect_equal(state_env$cache_format_version, "v1")
  expect_type(state_env$cache_root, "character")
  expect_type(state_env$cache_dir, "character")
  expect_type(state_env$package_version, "character")
})

test_that("Tier 1 cache helpers round-trip and clear values", {
  session_ns(".bef_cache_clear_session")()

  expect_false(session_ns(".bef_cache_exists")("unit-key"))
  value <- list(model = "dummy")
  expect_invisible(session_ns(".bef_cache_set")("unit-key", value))
  expect_true(session_ns(".bef_cache_exists")("unit-key"))
  expect_identical(session_ns(".bef_cache_get")("unit-key"), value)
  expect_true("unit-key" %in% session_ns(".bef_cache_keys")())

  expect_invisible(session_ns(".bef_cache_remove")("unit-key"))
  expect_false(session_ns(".bef_cache_exists")("unit-key"))

  session_ns(".bef_cache_set")("a", 1)
  session_ns(".bef_cache_set")("b", 2)
  cleared <- session_ns(".bef_cache_clear_session")()
  expect_identical(cleared, 2L)
  expect_identical(session_ns(".bef_cache_keys")(), character())
})

test_that("cache path and compile-option defaults are stable", {
  expect_equal(session_ns(".bef_cache_format_version")(), "v1")
  expect_equal(
    session_ns(".bef_cache_root")(),
    normalizePath(
      tools::R_user_dir("bayesEfron", which = "cache"),
      winslash = "/",
      mustWork = FALSE
    )
  )

  cache_dir <- session_ns(".bef_cache_dir")(create = FALSE)
  expect_equal(basename(cache_dir), "v1")
  expect_equal(dirname(cache_dir), session_ns(".bef_cache_root")())

  lock_path <- session_ns(".bef_cache_lock_path")()
  expect_equal(basename(lock_path), ".lock")
  expect_equal(basename(dirname(lock_path)), "v1")

  expect_identical(session_ns(".bef_cache_lock_timeout")(), 600L)
  expect_null(session_ns(".bef_default_cpp_options")())
  expect_identical(session_ns(".bef_default_stanc_options")(), list("O1"))
})

test_that("BAYESEFRON_CACHE_ROOT overrides cache root for tests", {
  old <- Sys.getenv("BAYESEFRON_CACHE_ROOT", unset = NA_character_)
  on.exit({
    if (is.na(old)) {
      Sys.unsetenv("BAYESEFRON_CACHE_ROOT")
    } else {
      Sys.setenv(BAYESEFRON_CACHE_ROOT = old)
    }
  }, add = TRUE)

  root <- tempfile("bayesefron-cache-root-")
  Sys.setenv(BAYESEFRON_CACHE_ROOT = root)

  expected_root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  expect_equal(session_ns(".bef_cache_root")(), expected_root)
  expect_equal(
    session_ns(".bef_cache_dir")(create = FALSE),
    file.path(expected_root, "v1")
  )

  created <- session_ns(".bef_cache_dir")(create = TRUE)
  expect_true(dir.exists(created))
  expect_equal(
    created,
    normalizePath(file.path(expected_root, "v1"), winslash = "/", mustWork = TRUE)
  )
})

test_that("message session state can be reset independently of Tier 1 cache", {
  session_ns(".bef_cache_clear_session")()
  session_ns(".bef_reset_session_messages")()

  session_ns(".bef_cache_set")("unit-key", list(model = "dummy"))
  msgs <- session_ns(".bef_session_msgs_env")()
  msgs[["kl_target_experimental"]] <- TRUE

  expect_true(session_ns(".bef_cache_exists")("unit-key"))
  expect_true(isTRUE(msgs[["kl_target_experimental"]]))

  expect_invisible(session_ns(".bef_reset_session_messages")())

  expect_true(session_ns(".bef_cache_exists")("unit-key"))
  expect_false(isTRUE(msgs[["kl_target_experimental"]]))

  session_ns(".bef_cache_clear_session")()
})
