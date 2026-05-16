fit67_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

fit67_theta_hat <- function() {
  c(site_a = -0.45, site_b = -0.1, site_c = 0, site_d = 0.25, site_e = 0.55)
}

fit67_sigma <- function() {
  c(0.12, 0.18, 0.15, 0.22, 0.2)
}

fit67_context <- function(...,
                          theta_hat = fit67_theta_hat(),
                          sigma = fit67_sigma(),
                          seed = 987L,
                          chains = 2L,
                          iter_warmup = 7L,
                          iter_sampling = 11L,
                          adapt_delta = 0.95) {
  fit67_ns(".bef_fit_prepare_context")(
    call = quote(bayes_efron_fit(theta_hat, sigma)),
    theta_hat = theta_hat,
    sigma = sigma,
    ...,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    seed = seed,
    check_installed = FALSE
  )
}

fit67_clock <- function(elapsed = 2.5) {
  times <- as.POSIXct("2026-05-11 12:00:00", tz = "UTC") + c(0, elapsed)
  index <- 0L
  function() {
    index <<- index + 1L
    times[[min(index, length(times))]]
  }
}

fit67_fake_model <- function(record = new.env(parent = emptyenv()),
                             sample_error = NULL,
                             chains_completed = NULL,
                             interrupt = NULL) {
  record$sample_calls <- list()
  structure(
    list(
      sample = function(...) {
        args <- list(...)
        record$sample_calls[[length(record$sample_calls) + 1L]] <- args
        if (!is.null(interrupt)) {
          stop(interrupt)
        }
        if (!is.null(sample_error)) {
          stop(sample_error, call. = FALSE)
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

fit67_model_fun <- function(model, record = new.env(parent = emptyenv())) {
  force(model)
  force(record)
  function(...) {
    record$model_args <- list(...)
    model
  }
}

fit67_run <- function(context = fit67_context(),
                      model,
                      model_record = new.env(parent = emptyenv()),
                      elapsed = 2.5,
                      interactive = FALSE,
                      sample_fun = NULL) {
  fit67_ns(".bef_fit_run_stages_6_7")(
    context,
    model_fun = fit67_model_fun(model, model_record),
    sample_fun = sample_fun,
    now = fit67_clock(elapsed),
    interactive_fun = function() interactive
  )
}

expect_fit67_sampling_failed <- function(err) {
  expect_s3_class(err, "bef_sampling_failed")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
  expect_equal(err$stage, 7L)
}

test_that("Stage 6 retrieves the cache-backed model and Stage 7 samples with canonical args", {
  withr::local_envvar(c(BAYESEFRON_PARALLEL_CHAINS = NA, PARALLEL_CHAINS = NA))
  context <- fit67_context()
  sample_record <- new.env(parent = emptyenv())
  model_record <- new.env(parent = emptyenv())
  model <- fit67_fake_model(sample_record, chains_completed = 2L)

  out <- fit67_run(
    context = context,
    model = model,
    model_record = model_record,
    elapsed = 2.5,
    interactive = FALSE
  )

  expect_s3_class(out, "bef_fit_context")
  expect_identical(model_record$model_args$model_name, "RE")
  expect_false(model_record$model_args$check_installed)
  expect_identical(out$model, model)
  expect_s3_class(out$cmdstan_fit, "fake_cmdstan_mcmc")
  expect_equal(out$runtime_seconds, 2.5, tolerance = 1e-8)

  expect_length(sample_record$sample_calls, 1L)
  args <- sample_record$sample_calls[[1L]]
  expect_identical(args$data, context$stan_data)
  expect_identical(args$chains, 2L)
  expect_identical(args$parallel_chains, 2L)
  expect_identical(args$iter_warmup, 7L)
  expect_identical(args$iter_sampling, 11L)
  expect_identical(args$adapt_delta, 0.95)
  expect_identical(args$seed, 987L)
  expect_identical(args$refresh, 0L)
  expect_equal(args$init, 0.5)
})

test_that("Stage 7 can limit CmdStan parallel chains for distributed runs", {
  withr::local_envvar(c(BAYESEFRON_PARALLEL_CHAINS = "1", PARALLEL_CHAINS = NA))
  sample_record <- new.env(parent = emptyenv())
  model <- fit67_fake_model(sample_record, chains_completed = 4L)

  fit67_run(
    context = fit67_context(chains = 4L),
    model = model
  )

  expect_identical(sample_record$sample_calls[[1L]]$chains, 4L)
  expect_identical(sample_record$sample_calls[[1L]]$parallel_chains, 1L)
})

test_that("Stage 7 rejects invalid CmdStan parallel chain overrides", {
  withr::local_envvar(c(BAYESEFRON_PARALLEL_CHAINS = "5", PARALLEL_CHAINS = NA))
  sample_record <- new.env(parent = emptyenv())
  model <- fit67_fake_model(sample_record, chains_completed = 4L)

  err <- tryCatch(
    fit67_run(
      context = fit67_context(chains = 4L),
      model = model
    ),
    error = identity
  )

  expect_fit67_sampling_failed(err)
  expect_match(conditionMessage(err), "BAYESEFRON_PARALLEL_CHAINS", fixed = TRUE)
  expect_length(sample_record$sample_calls, 0L)
})

test_that("Stage 7 uses interactive refresh when the session is interactive", {
  withr::local_envvar(c(BAYESEFRON_PARALLEL_CHAINS = NA, PARALLEL_CHAINS = NA))
  sample_record <- new.env(parent = emptyenv())
  model <- fit67_fake_model(sample_record, chains_completed = 2L)

  fit67_run(model = model, interactive = TRUE)

  expect_identical(sample_record$sample_calls[[1L]]$refresh, 200L)
})

test_that("Stage 6 cache retrieval errors propagate without sampling wrap", {
  context <- fit67_context()
  called_sample <- FALSE

  err <- tryCatch(
    fit67_ns(".bef_fit_run_stages_6_7")(
      context,
      model_fun = function(...) {
        fit67_ns(".bef_abort_lock_timeout")(
          "fake lock timeout",
          stage = 6L
        )
      },
      sample_fun = function(...) {
        called_sample <<- TRUE
        list()
      },
      now = fit67_clock(),
      interactive_fun = function() FALSE
    ),
    error = identity
  )

  expect_s3_class(err, "bef_lock_timeout")
  expect_s3_class(err, "bef_cache_error")
  expect_s3_class(err, "bef_error")
  expect_false(called_sample)
})

test_that("Stage 6 wraps untyped model retrieval failures as compile failures", {
  context <- fit67_context()

  err <- tryCatch(
    fit67_ns(".bef_fit_run_stages_6_7")(
      context,
      model_fun = function(...) {
        stop("untyped model provider failure", call. = FALSE)
      },
      now = fit67_clock(),
      interactive_fun = function() FALSE
    ),
    error = identity
  )

  expect_s3_class(err, "bef_compile_failed")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
  expect_equal(err$stage, 6L)
  expect_equal(err$model_family, "RE")
  expect_s3_class(err$parent, "simpleError")
})

test_that("Stage 7 wraps sampler errors as typed sampling failures", {
  model <- fit67_fake_model(sample_error = "fake sample failure")

  err <- tryCatch(fit67_run(model = model), error = identity)

  expect_fit67_sampling_failed(err)
  expect_s3_class(err$parent, "simpleError")
})

test_that("Stage 7 propagates sampler warnings through the standard warning buffer", {
  sample_record <- new.env(parent = emptyenv())
  model <- fit67_fake_model(sample_record, chains_completed = 2L)

  expect_warning(
    out <- 
    fit67_run(
      model = model,
      sample_fun = function(...) {
        warning("sampler warning", call. = FALSE)
        model$sample(...)
      }
    ),
    "sampler warning"
  )

  expect_s3_class(out$cmdstan_fit, "fake_cmdstan_mcmc")
})

test_that("Stage 7 reports missing sample method as a typed sampling failure", {
  model <- structure(list(not_sample = TRUE), class = "fake_cmdstan_model")

  err <- tryCatch(fit67_run(model = model), error = identity)

  expect_fit67_sampling_failed(err)
  expect_equal(err$model_family, "RE")
})

test_that("Stage 7 propagates user interrupts without bef wrapping", {
  interruption <- structure(
    list(message = "user interrupt", call = NULL),
    class = c("interrupt", "condition")
  )
  model <- fit67_fake_model(interrupt = interruption)

  err <- tryCatch(
    fit67_run(model = model),
    interrupt = identity,
    error = identity
  )

  expect_s3_class(err, "interrupt")
  expect_false(inherits(err, "bef_error"))
})

test_that("Stage 7 detects partial chains", {
  model <- fit67_fake_model(chains_completed = 1L)

  err <- tryCatch(fit67_run(model = model), error = identity)

  expect_s3_class(err, "bef_sampling_partial")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
  expect_equal(err$stage, 7L)
  expect_equal(err$chains_requested, 2L)
  expect_equal(err$chains_completed, 1L)
})

test_that("Stage 7 rejects zero or invalid completed-chain counts", {
  zero <- tryCatch(
    fit67_run(model = fit67_fake_model(chains_completed = 0L)),
    error = identity
  )
  expect_fit67_sampling_failed(zero)
  expect_equal(zero$chains_completed, 0L)

  invalid <- tryCatch(
    fit67_run(model = fit67_fake_model(chains_completed = c(1L, 2L))),
    error = identity
  )
  expect_fit67_sampling_failed(invalid)

  missing <- tryCatch(
    fit67_run(model = fit67_fake_model(chains_completed = NA_integer_)),
    error = identity
  )
  expect_fit67_sampling_failed(missing)
})

test_that("Stage 7 wraps completed-chain inspection failures", {
  model <- structure(
    list(
      sample = function(...) {
        structure(
          list(
            num_chains_completed = function() {
              stop("chain count unavailable", call. = FALSE)
            }
          ),
          class = "fake_cmdstan_mcmc"
        )
      }
    ),
    class = "fake_cmdstan_model"
  )

  err <- tryCatch(fit67_run(model = model), error = identity)

  expect_fit67_sampling_failed(err)
  expect_s3_class(err$parent, "simpleError")
})

test_that("Stage 7 passes auto-generated and zero seeds to CmdStan explicitly", {
  auto_context <- fit67_ns(".bef_fit_prepare_context")(
    call = quote(bayes_efron_fit(theta_hat, sigma)),
    theta_hat = fit67_theta_hat(),
    sigma = fit67_sigma(),
    seed = NULL,
    check_installed = FALSE,
    now = function() as.POSIXct("2026-05-11 12:00:03", tz = "UTC")
  )
  auto_record <- new.env(parent = emptyenv())
  fit67_run(
    context = auto_context,
    model = fit67_fake_model(auto_record, chains_completed = 4L)
  )
  expect_identical(auto_record$sample_calls[[1L]]$seed, auto_context$effective_seed)

  zero_context <- fit67_context(seed = 0L)
  zero_record <- new.env(parent = emptyenv())
  fit67_run(
    context = zero_context,
    model = fit67_fake_model(zero_record, chains_completed = 2L)
  )
  expect_identical(zero_record$sample_calls[[1L]]$seed, 0L)
})
