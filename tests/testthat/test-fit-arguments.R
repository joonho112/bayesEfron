prep_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

prep_theta_hat <- function() {
  c(site_a = -0.45, site_b = -0.1, site_c = 0, site_d = 0.25, site_e = 0.55)
}

prep_sigma <- function() {
  c(0.12, 0.18, 0.15, 0.22, 0.2)
}

prep_context <- function(...,
                          theta_hat = prep_theta_hat(),
                          sigma = prep_sigma(),
                          check_installed = FALSE,
                          now = function() as.POSIXct(
                            "2026-05-11 12:34:56",
                            tz = "UTC"
                          )) {
  prep_ns(".bef_fit_prepare_context")(
    call = quote(bayes_efron_fit(theta_hat, sigma)),
    theta_hat = theta_hat,
    sigma = sigma,
    ...,
    check_installed = check_installed,
    now = now
  )
}

expect_prep_invalid_args <- function(err) {
  expect_s3_class(err, "bef_invalid_args")
  expect_s3_class(err, "bef_pipeline_error")
  expect_s3_class(err, "bef_error")
}

test_that("the fit context holds the call, the normalized arguments, the data, the grid and the Stan data", {
  ctx <- prep_context(seed = 123L)

  expect_s3_class(ctx, "bef_fit_context")
  expect_equal(ctx$call, quote(bayes_efron_fit(theta_hat, sigma)))
  expect_equal(ctx$args$grid_method, "paper_realdata")
  expect_equal(ctx$args$model_family, "RE")
  expect_identical(ctx$args$seed, 123L)
  expect_identical(ctx$effective_seed, 123L)

  expect_s3_class(ctx$bef_data, "bef_data")
  expect_equal(ctx$bef_data$names, names(prep_theta_hat()))
  expect_equal(ctx$bef_data$source, "list")

  expect_named(ctx$grid, c("grid", "B", "M", "L", "expansion", "kappa", "grid_method", "attribution"))
  expect_identical(ctx$grid$L, 101L)
  expect_identical(ctx$grid$M, 6L)
  expect_equal(ctx$grid$grid_method, "paper_realdata")

  expect_named(ctx$stan_data, c("K", "theta_hat", "sigma", "L", "grid", "M", "B", "store_grid_quantities"))
  expect_identical(ctx$stan_data$K, 5L)
  expect_identical(ctx$stan_data$L, 101L)
  expect_identical(ctx$stan_data$M, 6L)
  expect_equal(ctx$stan_data$theta_hat, unname(prep_theta_hat()))
  expect_equal(ctx$stan_data$sigma, prep_sigma())
  expect_equal(dim(ctx$stan_data$B), c(101L, 6L))
})

test_that("a missing seed is generated from the clock", {
  fixed_time <- as.POSIXct("2026-05-11 12:34:56", tz = "UTC")
  ctx <- prep_context(seed = NULL, now = function() fixed_time)

  expect_null(ctx$args$seed)
  expect_identical(
    ctx$effective_seed,
    as.integer(as.numeric(fixed_time)) %% .Machine$integer.max
  )
})

test_that("at least one sampling iteration is required before sampling", {
  err <- tryCatch(prep_context(iter_sampling = 0L), error = identity)
  expect_prep_invalid_args(err)
  expect_identical(err$arg, "iter_sampling")
  expect_identical(prep_context(iter_sampling = 1L)$args$iter_sampling, 1L)
})

test_that("arguments are validated before cmdstanr is looked for", {
  called <- FALSE
  err <- tryCatch(
    prep_ns(".bef_fit_prepare_context")(
      call = quote(bayes_efron_fit(theta_hat, sigma)),
      theta_hat = c(1, 2, 3, 4, 5),
      sigma = c(0.1, 0.2, 0.1, 0.2, 0),
      check_installed = TRUE,
      check_installed_fun = function() {
        called <<- TRUE
        stop("cmdstanr check should not run", call. = FALSE)
      }
    ),
    error = identity
  )

  expect_prep_invalid_args(err)
  expect_equal(err$arg, "sigma")
  expect_false(called)
})

test_that("cmdstanr is looked for once the arguments are valid", {
  called <- 0L

  ctx <- prep_ns(".bef_fit_prepare_context")(
    call = quote(bayes_efron_fit(theta_hat, sigma)),
    theta_hat = prep_theta_hat(),
    sigma = prep_sigma(),
    check_installed = TRUE,
    check_installed_fun = function() {
      called <<- called + 1L
      invisible(TRUE)
    }
  )

  expect_s3_class(ctx, "bef_fit_context")
  expect_identical(called, 1L)
})

test_that("a grid method that needs theta_true fails without it", {
  err <- tryCatch(
    prep_context(grid_method = "paper_simulation"),
    error = identity
  )

  expect_s3_class(err, "bef_err_grid_oracle_required")
  expect_s3_class(err, "bef_grid_error")
  expect_s3_class(err, "bef_error")
})

test_that("fewer than five sites are rejected with the other argument checks", {
  for (K in 2:4) {
    err <- tryCatch(
      prep_context(theta_hat = seq(-0.1, 0.1, length.out = K), sigma = rep(0.2, K)),
      error = identity
    )
    expect_prep_invalid_args(err)
    expect_identical(err$arg, "theta_hat")
  }
  expect_s3_class(
    prep_context(theta_hat = seq(-0.2, 0.2, length.out = 5L), sigma = rep(0.2, 5L)),
    "bef_fit_context"
  )
})

test_that("Stan data are prepared for the simulation and sensitivity grid methods", {
  theta_hat <- prep_theta_hat()
  sigma <- prep_sigma()
  theta_true <- theta_hat + 0.05

  sim <- prep_context(
    theta_hat = theta_hat,
    sigma = sigma,
    grid_method = "paper_simulation",
    theta_true = theta_true
  )
  sens <- prep_context(
    theta_hat = theta_hat,
    sigma = sigma,
    grid_method = "paper_sensitivity",
    theta_true = theta_true,
    bound_expansion = 0.25
  )

  expect_equal(sim$args$theta_true, theta_true)
  expect_equal(sim$grid$grid_method, "paper_simulation")
  expect_named(sim$stan_data, c("K", "theta_hat", "sigma", "L", "grid", "M", "B", "store_grid_quantities"))

  expect_equal(sens$args$bound_expansion, 0.25)
  expect_equal(sens$grid$grid_method, "paper_sensitivity")
  expect_named(sens$stan_data, c("K", "theta_hat", "sigma", "L", "grid", "M", "B", "store_grid_quantities"))
})

test_that("the fit context carries everything that sampling needs", {
  ctx <- prep_context(seed = 123L)

  expect_named(
    ctx,
    c("call", "args", "effective_seed", "grid", "bef_data", "stan_data")
  )
  expect_false("cmdstan_fit" %in% names(ctx))
  expect_false("runtime_seconds" %in% names(ctx))
})

test_that("max_treedepth must be a whole number between 1 and 20", {
  for (bad in list(0L, 21L, 10.5, c(8L, 10L), "10", NA_integer_)) {
    err <- tryCatch(prep_context(max_treedepth = bad), error = identity)
    expect_prep_invalid_args(err)
    expect_identical(err$arg, "max_treedepth")
  }
  expect_identical(prep_context()$args$max_treedepth, 10L)
  expect_identical(prep_context(max_treedepth = 15)$args$max_treedepth, 15L)
})
