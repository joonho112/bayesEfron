samp_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

samp_theta_hat <- function() {
  c(site_a = -0.45, site_b = -0.1, site_c = 0, site_d = 0.25, site_e = 0.55)
}

samp_sigma <- function() {
  c(0.12, 0.18, 0.15, 0.22, 0.2)
}

samp_context <- function(...,
                          theta_hat = samp_theta_hat(),
                          sigma = samp_sigma(),
                          seed = 987L,
                          chains = 2L,
                          iter_warmup = 7L,
                          iter_sampling = 11L,
                          adapt_delta = 0.95) {
  samp_ns(".bef_fit_prepare_context")(
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

samp_clock <- function(elapsed = 2.5) {
  times <- as.POSIXct("2026-05-11 12:00:00", tz = "UTC") + c(0, elapsed)
  index <- 0L
  function() {
    index <<- index + 1L
    times[[min(index, length(times))]]
  }
}

samp_fake_model <- function(record = new.env(parent = emptyenv()),
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

samp_model_fun <- function(model, record = new.env(parent = emptyenv())) {
  force(model)
  force(record)
  function(...) {
    record$model_args <- list(...)
    model
  }
}

samp_run <- function(context = samp_context(),
                      model,
                      model_record = new.env(parent = emptyenv()),
                      elapsed = 2.5,
                      interactive = FALSE,
                      sample_fun = NULL) {
  samp_ns(".bef_fit_run_sampler")(
    context,
    model_fun = samp_model_fun(model, model_record),
    sample_fun = sample_fun,
    now = samp_clock(elapsed),
    interactive_fun = function() interactive
  )
}

expect_samp_sampling_failed <- function(err) {
  expect_s3_class(err, "bef_sampling_failed")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
}

test_that("the model is retrieved and sampled with the arguments of the fit call", {
  context <- samp_context()
  sample_record <- new.env(parent = emptyenv())
  model_record <- new.env(parent = emptyenv())
  model <- samp_fake_model(sample_record, chains_completed = 2L)

  out <- samp_run(
    context = context,
    model = model,
    model_record = model_record,
    elapsed = 2.5,
    interactive = FALSE
  )

  expect_s3_class(out, "bef_fit_context")
  expect_identical(model_record$model_args, list(model_family = "RE"))
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
  expect_identical(args$max_treedepth, 10L)
  expect_identical(args$seed, 987L)
  expect_identical(args$refresh, 0L)
  expect_equal(args$init, 0.5)
  expect_identical(out$sampler_settings, args[names(args) != "data"])
})

test_that("sampling leaves R's random number generator as it found it", {
  withr::local_seed(20261001)
  before <- .Random.seed

  samp_run(model = samp_fake_model(chains_completed = 2L))
  expect_identical(.Random.seed, before)
})

test_that("max_treedepth is passed to the sampler", {
  sample_record <- new.env(parent = emptyenv())
  model <- samp_fake_model(sample_record, chains_completed = 2L)

  samp_run(context = samp_context(max_treedepth = 12L), model = model)
  expect_identical(sample_record$sample_calls[[1L]]$max_treedepth, 12L)
})

test_that("the number of parallel chains can be limited", {
  sample_record <- new.env(parent = emptyenv())
  model <- samp_fake_model(sample_record, chains_completed = 4L)

  samp_run(
    context = samp_context(chains = 4L, parallel_chains = 1L),
    model = model
  )

  expect_identical(sample_record$sample_calls[[1L]]$chains, 4L)
  expect_identical(sample_record$sample_calls[[1L]]$parallel_chains, 1L)
})

test_that("parallel_chains defaults to chains and cannot exceed it", {
  expect_identical(samp_context(chains = 3L)$args$parallel_chains, 3L)

  for (bad in list(0L, 5L, 1.5, "2")) {
    err <- tryCatch(
      samp_context(chains = 4L, parallel_chains = bad),
      error = identity
    )
    expect_s3_class(err, "bef_invalid_args")
    expect_identical(err$arg, "parallel_chains")
  }
})

test_that("progress is printed only in interactive sessions", {
  sample_record <- new.env(parent = emptyenv())
  model <- samp_fake_model(sample_record, chains_completed = 2L)

  samp_run(model = model, interactive = TRUE)

  expect_identical(sample_record$sample_calls[[1L]]$refresh, 200L)
})

test_that("package errors from model retrieval are passed on unchanged", {
  context <- samp_context()
  called_sample <- FALSE

  err <- tryCatch(
    samp_ns(".bef_fit_run_sampler")(
      context,
      model_fun = function(...) {
        samp_ns(".bef_abort_compile_failed")("fake compile failure")
      },
      sample_fun = function(...) {
        called_sample <<- TRUE
        list()
      },
      now = samp_clock(),
      interactive_fun = function() FALSE
    ),
    error = identity
  )

  expect_s3_class(err, "bef_compile_failed")
  expect_match(conditionMessage(err), "fake compile failure")
  expect_null(err$model_family)
  expect_false(called_sample)
})

test_that("other model retrieval errors become bef_compile_failed", {
  context <- samp_context()

  err <- tryCatch(
    samp_ns(".bef_fit_run_sampler")(
      context,
      model_fun = function(...) {
        stop("untyped model provider failure", call. = FALSE)
      },
      now = samp_clock(),
      interactive_fun = function() FALSE
    ),
    error = identity
  )

  expect_s3_class(err, "bef_compile_failed")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
  expect_equal(err$model_family, "RE")
  expect_s3_class(err$parent, "simpleError")
})

test_that("sampler errors become bef_sampling_failed", {
  model <- samp_fake_model(sample_error = "fake sample failure")

  err <- tryCatch(samp_run(model = model), error = identity)

  expect_samp_sampling_failed(err)
  expect_s3_class(err$parent, "simpleError")
})

test_that("sampler warnings are passed on", {
  sample_record <- new.env(parent = emptyenv())
  model <- samp_fake_model(sample_record, chains_completed = 2L)

  expect_warning(
    out <- 
    samp_run(
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

test_that("a model without a sample method is a sampling failure", {
  model <- structure(list(not_sample = TRUE), class = "fake_cmdstan_model")

  err <- tryCatch(samp_run(model = model), error = identity)

  expect_samp_sampling_failed(err)
  expect_equal(err$model_family, "RE")
})

test_that("user interrupts are passed on unchanged", {
  interruption <- structure(
    list(message = "user interrupt", call = NULL),
    class = c("interrupt", "condition")
  )
  model <- samp_fake_model(interrupt = interruption)

  err <- tryCatch(
    samp_run(model = model),
    interrupt = identity,
    error = identity
  )

  expect_s3_class(err, "interrupt")
  expect_false(inherits(err, "bef_error"))
})

test_that("fewer completed chains than requested is an error", {
  model <- samp_fake_model(chains_completed = 1L)

  err <- tryCatch(samp_run(model = model), error = identity)

  expect_s3_class(err, "bef_sampling_partial")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
  expect_equal(err$chains_requested, 2L)
  expect_equal(err$chains_completed, 1L)
})

test_that("zero or invalid completed-chain counts are errors", {
  zero <- tryCatch(
    samp_run(model = samp_fake_model(chains_completed = 0L)),
    error = identity
  )
  expect_samp_sampling_failed(zero)
  expect_equal(zero$chains_completed, 0L)

  invalid <- tryCatch(
    samp_run(model = samp_fake_model(chains_completed = c(1L, 2L))),
    error = identity
  )
  expect_samp_sampling_failed(invalid)

  missing <- tryCatch(
    samp_run(model = samp_fake_model(chains_completed = NA_integer_)),
    error = identity
  )
  expect_samp_sampling_failed(missing)
})

test_that("a failure while counting completed chains is a sampling failure", {
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

  err <- tryCatch(samp_run(model = model), error = identity)

  expect_samp_sampling_failed(err)
  expect_s3_class(err$parent, "simpleError")
})

test_that("generated and zero seeds are passed to CmdStan", {
  auto_context <- samp_ns(".bef_fit_prepare_context")(
    call = quote(bayes_efron_fit(theta_hat, sigma)),
    theta_hat = samp_theta_hat(),
    sigma = samp_sigma(),
    seed = NULL,
    check_installed = FALSE,
    now = function() as.POSIXct("2026-05-11 12:00:03", tz = "UTC")
  )
  auto_record <- new.env(parent = emptyenv())
  samp_run(
    context = auto_context,
    model = samp_fake_model(auto_record, chains_completed = 4L)
  )
  expect_identical(auto_record$sample_calls[[1L]]$seed, auto_context$effective_seed)

  zero_context <- samp_context(seed = 0L)
  zero_record <- new.env(parent = emptyenv())
  samp_run(
    context = zero_context,
    model = samp_fake_model(zero_record, chains_completed = 2L)
  )
  expect_identical(zero_record$sample_calls[[1L]]$seed, 0L)
})
