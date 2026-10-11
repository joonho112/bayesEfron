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
  expect_identical(model_record$model_args, list(model_family = "RE", program = "efron_mix"))
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

test_that("the program of this version is sampled with twelve significant digits", {
  sample_record <- new.env(parent = emptyenv())
  model_record <- new.env(parent = emptyenv())
  model <- samp_fake_model(sample_record, chains_completed = 2L)
  context <- samp_context()

  out <- samp_run(context = context, model = model, model_record = model_record)
  expect_identical(model_record$model_args, list(model_family = "RE", program = "efron_mix"))
  args <- sample_record$sample_calls[[1L]]
  expect_identical(args$sig_figs, 12L)
  expect_identical(args$data$tau_c, 0)
  expect_identical(out$sampler_settings, args[names(args) != "data"])

  # The call for the program of version 0.3.0 is the call of that version:
  # it has no such argument.
  samp_run(context = samp_context(settings = "0.3.0"), model = model,
           model_record = model_record)
  expect_identical(model_record$model_args, list(model_family = "RE", program = "efron_re"))
  expect_length(sample_record$sample_calls, 2L)
  old <- sample_record$sample_calls[[2L]]
  expect_identical(
    old[c("chains", "parallel_chains", "iter_warmup", "iter_sampling", "adapt_delta",
          "max_treedepth", "seed", "refresh", "init")],
    list(chains = 2L, parallel_chains = 2L, iter_warmup = 7L, iter_sampling = 11L,
         adapt_delta = 0.95, max_treedepth = 10L, seed = 987L, refresh = 0L, init = 0.5)
  )
  expect_false("sig_figs" %in% names(old))
  expect_false("tau_c" %in% names(old$data))
})

test_that("a model function that takes the family only is called without the program", {
  received <- NULL
  model <- samp_fake_model(chains_completed = 2L)
  samp_ns(".bef_fit_run_sampler")(
    samp_context(),
    model_fun = function(model_family) {
      received <<- model_family
      model
    },
    now = samp_clock(1),
    interactive_fun = function() FALSE
  )
  expect_identical(received, "RE")
})

test_that("refresh, init and show_messages reach the sampler, and thin does not", {
  sample_record <- new.env(parent = emptyenv())
  model <- samp_fake_model(sample_record, chains_completed = 2L)
  call_of <- function(..., interactive = FALSE) {
    samp_run(context = samp_context(...), model = model, interactive = interactive)
    sample_record$sample_calls[[length(sample_record$sample_calls)]]
  }

  # Outside an interactive session nothing is printed unless it is asked for.
  quiet <- call_of()
  expect_identical(quiet$refresh, 0L)
  expect_false(quiet$show_messages)
  expect_false(quiet$show_exceptions)
  expect_true("diagnostics" %in% names(quiet))
  expect_null(quiet$diagnostics)
  expect_identical(quiet$init, 0.5)
  expect_false("thin" %in% names(quiet))

  loud <- call_of(show_messages = TRUE, refresh = 50L)
  expect_identical(loud$refresh, 50L)
  expect_false(any(c("show_messages", "show_exceptions", "diagnostics") %in% names(loud)))

  # In an interactive session CmdStan reports its progress.
  shown <- call_of(interactive = TRUE)
  expect_identical(shown$refresh, 200L)
  expect_false("show_messages" %in% names(shown))
  hushed <- call_of(interactive = TRUE, show_messages = FALSE, refresh = 0L)
  expect_identical(hushed$refresh, 0L)
  expect_false(hushed$show_messages)

  starts <- list(list(alpha = rep(0, 6), lambda = 1), list(alpha = rep(0.1, 6), lambda = 2))
  expect_identical(call_of(init = starts)$init, starts)
  expect_identical(call_of(init = 0)$init, 0)
  fun <- function() list(alpha = rep(0, 6), lambda = 1)
  expect_identical(call_of(init = fun)$init, fun)

  # Thinning is done after sampling.
  thinned <- call_of(thin = 5L)
  expect_false("thin" %in% names(thinned))
  expect_identical(thinned$iter_sampling, 11L)

  for (bad in list(list(refresh = -1L), list(refresh = 1.5), list(refresh = "0"),
                   list(init = -1), list(init = c(0.5, 1)), list(init = "random"),
                   list(init = list(list(alpha = 0))), list(thin = 0L), list(thin = 12L),
                   list(thin = 2.5), list(show_messages = NA), list(show_messages = "yes"))) {
    err <- tryCatch(do.call(samp_context, bad), error = identity)
    expect_s3_class(err, "bef_invalid_args")
    expect_identical(err$arg, names(bad))
  }
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
