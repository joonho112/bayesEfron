fit15_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

fit15_theta_hat <- function() {
  c(site_a = -0.45, site_b = -0.1, site_c = 0, site_d = 0.25, site_e = 0.55)
}

fit15_sigma <- function() {
  c(0.12, 0.18, 0.15, 0.22, 0.2)
}

fit15_context <- function(...,
                          theta_hat = fit15_theta_hat(),
                          sigma = fit15_sigma(),
                          check_installed = FALSE,
                          now = function() as.POSIXct(
                            "2026-05-11 12:34:56",
                            tz = "UTC"
                          )) {
  fit15_ns(".bef_fit_prepare_context")(
    call = quote(bayes_efron_fit(theta_hat, sigma)),
    theta_hat = theta_hat,
    sigma = sigma,
    ...,
    check_installed = check_installed,
    now = now
  )
}

expect_fit15_invalid_args <- function(err) {
  expect_s3_class(err, "bef_invalid_args")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
}

test_that("fit entry exposes the canonical Stage 1-8 signature", {
  fmls <- formals(bayes_efron_fit)

  expect_identical(
    names(fmls),
    c(
      "theta_hat", "sigma", "...", "grid_method", "L", "expansion", "M",
      "theta_true", "bound_expansion", "model_family", "chains",
      "iter_warmup", "iter_sampling", "adapt_delta", "seed",
      "keep_cmdstan_fit"
    )
  )
  expect_identical(fmls$L, 101L)
  expect_identical(fmls$M, 6L)
  expect_identical(fmls$expansion, 0.5)
  expect_identical(fmls$model_family, "RE")
  expect_identical(fmls$chains, 4L)
  expect_identical(fmls$iter_warmup, 1000L)
  expect_identical(fmls$iter_sampling, 3000L)
  expect_identical(fmls$adapt_delta, 0.9)
  expect_null(fmls$seed)
  expect_false(fmls$keep_cmdstan_fit)
  expect_true("bayes_efron_fit" %in% getNamespaceExports("bayesEfron"))
})

test_that("Stage 1-5 context captures call, normalized args, data, grid, and Stan data", {
  ctx <- fit15_context(seed = 123L)

  expect_s3_class(ctx, "bef_fit_context")
  expect_equal(ctx$call, quote(bayes_efron_fit(theta_hat, sigma)))
  expect_equal(ctx$args$grid_method, "paper_realdata")
  expect_equal(ctx$args$model_family, "RE")
  expect_identical(ctx$args$seed, 123L)
  expect_identical(ctx$effective_seed, 123L)

  expect_s3_class(ctx$bef_data, "bef_data")
  expect_equal(ctx$bef_data$names, names(fit15_theta_hat()))
  expect_equal(ctx$bef_data$source, "list")

  expect_named(ctx$grid, c("grid", "B", "M", "L", "expansion", "kappa", "grid_method", "attribution"))
  expect_identical(ctx$grid$L, 101L)
  expect_identical(ctx$grid$M, 6L)
  expect_equal(ctx$grid$grid_method, "paper_realdata")

  expect_named(ctx$stan_data, c("K", "theta_hat", "sigma", "L", "grid", "M", "B"))
  expect_identical(ctx$stan_data$K, 5L)
  expect_identical(ctx$stan_data$L, 101L)
  expect_identical(ctx$stan_data$M, 6L)
  expect_equal(ctx$stan_data$theta_hat, unname(fit15_theta_hat()))
  expect_equal(ctx$stan_data$sigma, fit15_sigma())
  expect_equal(dim(ctx$stan_data$B), c(101L, 6L))
})

test_that("Stage 1-5 auto-generates reproducible effective seed from time", {
  fixed_time <- as.POSIXct("2026-05-11 12:34:56", tz = "UTC")
  ctx <- fit15_context(seed = NULL, now = function() fixed_time)

  expect_null(ctx$args$seed)
  expect_identical(
    ctx$effective_seed,
    as.integer(as.numeric(fixed_time)) %% .Machine$integer.max
  )
})

test_that("Stage 2 validation runs before the cmdstanr Stage 3 check", {
  called <- FALSE
  err <- tryCatch(
    fit15_ns(".bef_fit_prepare_context")(
      call = quote(bayes_efron_fit(theta_hat, sigma)),
      theta_hat = c(1, 2),
      sigma = c(0.1, 0),
      check_installed = TRUE,
      check_installed_fun = function() {
        called <<- TRUE
        stop("cmdstanr check should not run", call. = FALSE)
      }
    ),
    error = identity
  )

  expect_fit15_invalid_args(err)
  expect_equal(err$arg, "sigma")
  expect_equal(err$stage, 2L)
  expect_false(called)
})

test_that("Stage 3 check runs after valid Stage 2 validation", {
  called <- 0L

  ctx <- fit15_ns(".bef_fit_prepare_context")(
    call = quote(bayes_efron_fit(theta_hat, sigma)),
    theta_hat = fit15_theta_hat(),
    sigma = fit15_sigma(),
    check_installed = TRUE,
    check_installed_fun = function() {
      called <<- called + 1L
      invisible(TRUE)
    }
  )

  expect_s3_class(ctx, "bef_fit_context")
  expect_identical(called, 1L)
})

test_that("Stage 4 propagates oracle-required grid errors", {
  err <- tryCatch(
    fit15_context(grid_method = "paper_simulation"),
    error = identity
  )

  expect_s3_class(err, "bef_err_grid_oracle_required")
  expect_s3_class(err, "bef_grid_error")
  expect_s3_class(err, "bef_error")
})

test_that("Stage 4 bef_data validation enforces the class boundary", {
  err <- tryCatch(
    fit15_context(theta_hat = c(-0.1, 0, 0.1), sigma = rep(0.2, 3L)),
    error = identity
  )

  expect_fit15_invalid_args(err)
  expect_s3_class(err, "bayesEfron_validate_error")
})

test_that("Stage 5 prepares oracle and sensitivity Stan data explicitly", {
  theta_hat <- fit15_theta_hat()
  sigma <- fit15_sigma()
  theta_true <- theta_hat + 0.05

  sim <- fit15_context(
    theta_hat = theta_hat,
    sigma = sigma,
    grid_method = "paper_simulation",
    theta_true = theta_true
  )
  sens <- fit15_context(
    theta_hat = theta_hat,
    sigma = sigma,
    grid_method = "paper_sensitivity",
    theta_true = theta_true,
    bound_expansion = 0.25
  )

  expect_equal(sim$args$theta_true, theta_true)
  expect_equal(sim$grid$grid_method, "paper_simulation")
  expect_named(sim$stan_data, c("K", "theta_hat", "sigma", "L", "grid", "M", "B"))

  expect_equal(sens$args$bound_expansion, 0.25)
  expect_equal(sens$grid$grid_method, "paper_sensitivity")
  expect_named(sens$stan_data, c("K", "theta_hat", "sigma", "L", "grid", "M", "B"))
})

test_that("Stage 1-5 context remains ready for later pipeline stages", {
  ctx <- fit15_context(seed = 123L)

  expect_named(
    ctx,
    c("call", "args", "effective_seed", "grid", "bef_data", "stan_data")
  )
  expect_false("cmdstan_fit" %in% names(ctx))
  expect_false("runtime_seconds" %in% names(ctx))
})
