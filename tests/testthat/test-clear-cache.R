local_cache_with_files <- function(env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = env)
  withr::local_envvar(c(BAYESEFRON_CACHE_ROOT = root), .local_envir = env)
  writeLines("x", file.path(root, "efron_re-aaaaaaaaaaaa-cmdstan2.38.0-arm64"))
  writeLines("x", file.path(root, "efron_re-bbbbbbbbbbbb-cmdstan2.37.0-arm64"))
  writeLines("x", file.path(root, "unrelated.txt"))
  models <- getFromNamespace(".bef_models", "bayesEfron")
  assign("a model", "attached", envir = models)
  withr::defer(rm(list = ls(models, all.names = TRUE), envir = models), envir = env)
  root
}

attached_models <- function() {
  ls(getFromNamespace(".bef_models", "bayesEfron"), all.names = TRUE)
}

test_that("bayes_efron_clear_cache() removes the compiled models and nothing else", {
  root <- local_cache_with_files()

  expect_message(n <- bayes_efron_clear_cache(), "removed 2 cached files")
  expect_identical(n, 2L)
  expect_identical(list.files(root), "unrelated.txt")
  expect_length(attached_models(), 0L)

  expect_message(again <- bayes_efron_clear_cache(), "removed 0 cached files")
  expect_identical(again, 0L)
})

test_that("files left by the layout of earlier versions are removed too", {
  root <- local_cache_with_files()
  dir.create(file.path(root, "v1"))
  writeLines("x", file.path(root, "v1", paste(rep("a", 64L), collapse = "")))
  writeLines("{}", file.path(root, "v1", "meta.json"))

  n <- suppressMessages(bayes_efron_clear_cache())
  expect_identical(n, 4L)
  expect_false(dir.exists(file.path(root, "v1")))
})

test_that("the scope values of earlier versions still work, with a warning", {
  root <- local_cache_with_files()

  expect_warning(
    n <- suppressMessages(bayes_efron_clear_cache("lock_only")),
    "deprecated"
  )
  expect_identical(n, 0L)
  expect_length(list.files(root), 3L)
  expect_length(attached_models(), 1L)

  expect_warning(
    n <- suppressMessages(bayes_efron_clear_cache("session")),
    "deprecated"
  )
  expect_identical(n, 0L)
  expect_length(list.files(root), 3L)
  expect_length(attached_models(), 0L)

  expect_warning(
    n <- suppressMessages(bayes_efron_clear_cache("compiled_models")),
    "deprecated"
  )
  expect_identical(n, 2L)
  expect_identical(list.files(root), "unrelated.txt")
})

test_that("an unknown scope is an error, and an absent cache directory is not", {
  err <- tryCatch(bayes_efron_clear_cache("bogus"), error = identity)
  expect_s3_class(err, "bef_invalid_args")
  expect_identical(err$arg, "scope")

  withr::local_envvar(
    c(BAYESEFRON_CACHE_ROOT = file.path(withr::local_tempdir(), "absent"))
  )
  expect_identical(suppressMessages(bayes_efron_clear_cache()), 0L)
  expect_false(dir.exists(Sys.getenv("BAYESEFRON_CACHE_ROOT")))
})
