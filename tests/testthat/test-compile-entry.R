compile_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

new_compile_fake_model <- function(record = new.env(parent = emptyenv()),
                                   fail_sample = FALSE,
                                   chains_completed = NULL) {
  record$sample_calls <- list()
  structure(
    list(
      sample = function(...) {
        args <- list(...)
        record$sample_calls[[length(record$sample_calls) + 1L]] <- args
        record$output_dir_exists_at_sample <- dir.exists(args$output_dir)
        if (isTRUE(fail_sample)) {
          stop("fake smoke failure", call. = FALSE)
        }
        out <- list(args = args)
        if (!is.null(chains_completed)) {
          out$num_chains_completed <- function() chains_completed
        }
        structure(out, class = "fake_cmdstan_mcmc")
      }
    ),
    class = "fake_cmdstan_model"
  )
}

expect_compile_entry_invalid_args <- function(err, arg) {
  expect_s3_class(err, "bef_invalid_args")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
  expect_equal(err$module, "compile-entry")
  expect_equal(err$stage, 6L)
  expect_equal(err$arg, arg)
}

test_that("compile entry exposes the canonical four-argument signature", {
  fmls <- formals(bayes_efron_compile)

  expect_identical(
    names(fmls),
    c("model_family", "quiet", "force_recompile", "seed_for_check")
  )
  expect_identical(fmls$model_family, "RE")
  expect_true(identical(fmls$quiet, TRUE))
  expect_true(identical(fmls$force_recompile, FALSE))
  expect_identical(fmls$seed_for_check, 42L)
  expect_true("bayes_efron_compile" %in% getNamespaceExports("bayesEfron"))
})

test_that("compile entry validates arguments with typed errors", {
  fake_model <- new_compile_fake_model()
  model_fun <- function(...) fake_model
  smoke <- function(...) invisible(TRUE)
  entry <- compile_ns(".bef_compile_entry")

  err <- tryCatch(
    entry("HE", model_fun = model_fun, smoke_check_fun = smoke, check_installed = FALSE),
    error = identity
  )
  expect_compile_entry_invalid_args(err, "model_family")

  err <- tryCatch(
    entry(quiet = NA, model_fun = model_fun, smoke_check_fun = smoke, check_installed = FALSE),
    error = identity
  )
  expect_compile_entry_invalid_args(err, "quiet")

  err <- tryCatch(
    entry(
      force_recompile = c(TRUE, FALSE),
      model_fun = model_fun,
      smoke_check_fun = smoke,
      check_installed = FALSE
    ),
    error = identity
  )
  expect_compile_entry_invalid_args(err, "force_recompile")

  err <- tryCatch(
    entry(
      seed_for_check = NULL,
      model_fun = model_fun,
      smoke_check_fun = smoke,
      check_installed = FALSE
    ),
    error = identity
  )
  expect_compile_entry_invalid_args(err, "seed_for_check")
})

test_that("compile entry delegates to model cache and smoke check", {
  fake_model <- new_compile_fake_model()
  calls <- new.env(parent = emptyenv())
  calls$model <- NULL
  calls$smoke <- NULL

  model_fun <- function(model_name,
                        force_recompile,
                        cmdstan_model_fun,
                        check_installed,
                        ...) {
    calls$model <- list(
      model_name = model_name,
      force_recompile = force_recompile,
      cmdstan_model_fun = cmdstan_model_fun,
      check_installed = check_installed
    )
    fake_model
  }
  smoke <- function(model, model_family, seed_for_check, quiet) {
    calls$smoke <- list(
      model = model,
      model_family = model_family,
      seed_for_check = seed_for_check,
      quiet = quiet
    )
    invisible(TRUE)
  }

  out <- compile_ns(".bef_compile_entry")(
    model_family = "RE",
    quiet = FALSE,
    force_recompile = TRUE,
    seed_for_check = 123L,
    model_fun = model_fun,
    smoke_check_fun = smoke,
    check_installed = FALSE
  )

  expect_identical(out, fake_model)
  expect_equal(calls$model$model_name, "RE")
  expect_true(calls$model$force_recompile)
  expect_true(is.function(calls$model$cmdstan_model_fun))
  expect_false(calls$model$check_installed)
  expect_identical(calls$smoke$model, fake_model)
  expect_equal(calls$smoke$model_family, "RE")
  expect_identical(calls$smoke$seed_for_check, 123L)
  expect_false(calls$smoke$quiet)
})

test_that("quiet mode suppresses model and smoke-check messages", {
  fake_model <- new_compile_fake_model()
  noisy_model <- function(...) {
    message("compile message")
    cat("compile output\n")
    fake_model
  }
  noisy_smoke <- function(...) {
    message("smoke message")
    cat("smoke output\n")
    invisible(TRUE)
  }

  expect_silent(
    out <- compile_ns(".bef_compile_entry")(
      quiet = TRUE,
      model_fun = noisy_model,
      smoke_check_fun = noisy_smoke,
      check_installed = FALSE
    )
  )
  expect_identical(out, fake_model)

  messages <- character()
  output <- utils::capture.output(
    withCallingHandlers(
      compile_ns(".bef_compile_entry")(
        quiet = FALSE,
        model_fun = noisy_model,
        smoke_check_fun = noisy_smoke,
        check_installed = FALSE
      ),
      message = function(msg) {
        messages <<- c(messages, conditionMessage(msg))
        invokeRestart("muffleMessage")
      }
    )
  )
  expect_true(any(grepl("compile message", messages, fixed = TRUE)))
  expect_true(any(grepl("smoke message", messages, fixed = TRUE)))
  expect_true(any(grepl("compile output", output, fixed = TRUE)))
  expect_true(any(grepl("smoke output", output, fixed = TRUE)))
})

test_that("post-compile smoke check uses the synthetic five-site fixture", {
  record <- new.env(parent = emptyenv())
  model <- new_compile_fake_model(record)

  out <- compile_ns(".bef_compile_smoke_check")(
    model = model,
    model_family = "RE",
    seed_for_check = 99L,
    quiet = TRUE
  )

  expect_identical(out, model)
  expect_length(record$sample_calls, 1L)
  args <- record$sample_calls[[1L]]
  expect_named(
    args$data,
    c("K", "theta_hat", "sigma", "L", "grid", "M", "B")
  )
  expect_identical(args$data$K, 5L)
  expect_identical(args$data$L, 51L)
  expect_identical(args$data$M, 3L)
  expect_equal(dim(args$data$B), c(51L, 3L))
  expect_identical(args$chains, 1L)
  expect_identical(args$parallel_chains, 1L)
  expect_identical(args$iter_warmup, 2L)
  expect_identical(args$iter_sampling, 2L)
  expect_identical(args$seed, 99L)
  expect_identical(args$refresh, 0L)
  expect_equal(args$init, 0.5)
  expect_true(record$output_dir_exists_at_sample)
  expect_false(dir.exists(args$output_dir))
  expect_false(args$show_messages)
  expect_identical(args$diagnostics, c("divergences", "treedepth"))
})

test_that("post-compile smoke check failures are typed compile failures", {
  model_without_sample <- structure(list(), class = "fake_cmdstan_model")
  err <- tryCatch(
    compile_ns(".bef_compile_smoke_check")(
      model = model_without_sample,
      model_family = "RE",
      seed_for_check = 42L,
      quiet = TRUE
    ),
    error = identity
  )

  expect_s3_class(err, "bef_compile_failed")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
  expect_equal(err$model_family, "RE")
  expect_identical(err$seed_for_check, 42L)

  failing_model <- new_compile_fake_model(fail_sample = TRUE)
  err <- tryCatch(
    compile_ns(".bef_compile_smoke_check")(
      model = failing_model,
      model_family = "RE",
      seed_for_check = 42L,
      quiet = TRUE
    ),
    error = identity
  )

  expect_s3_class(err, "bef_compile_failed")
  expect_s3_class(err$parent, "simpleError")
  expect_equal(err$model_family, "RE")
  expect_identical(err$seed_for_check, 42L)

  incomplete_model <- new_compile_fake_model(chains_completed = 0L)
  err <- tryCatch(
    compile_ns(".bef_compile_smoke_check")(
      model = incomplete_model,
      model_family = "RE",
      seed_for_check = 42L,
      quiet = TRUE
    ),
    error = identity
  )

  expect_s3_class(err, "bef_compile_failed")
  expect_equal(err$model_family, "RE")
  expect_identical(err$seed_for_check, 42L)
  expect_identical(err$chains_completed, 0L)
})
