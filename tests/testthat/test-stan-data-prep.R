test_that("prepare_stan_data returns the exact locked Stan data shape", {
  theta_hat <- c(-1, 0, 2)
  sigma <- c(0.2, 0.3, 0.4)
  grid <- make_efron_grid(theta_hat, sigma)
  prepare_stan_data <- getFromNamespace("prepare_stan_data", "bayesEfron")

  stan_data <- prepare_stan_data(
    list(theta_hat = theta_hat, sigma = sigma),
    grid
  )

  expect_named(
    stan_data,
    c("K", "theta_hat", "sigma", "L", "grid", "M", "B")
  )
  expect_type(stan_data$K, "integer")
  expect_type(stan_data$L, "integer")
  expect_type(stan_data$M, "integer")
  expect_equal(stan_data$K, length(theta_hat))
  expect_equal(stan_data$theta_hat, theta_hat)
  expect_equal(stan_data$sigma, sigma)
  expect_equal(stan_data$L, grid$L)
  expect_equal(stan_data$grid, grid$grid)
  expect_equal(stan_data$M, grid$M)
  expect_equal(stan_data$B, unclass(grid$B), ignore_attr = TRUE)
  expect_equal(dim(stan_data$B), c(stan_data$L, stan_data$M))
})

test_that("prepare_stan_data rejects invalid inputs with typed conditions", {
  theta_hat <- c(-1, 0, 2)
  sigma <- c(0.2, 0.3, 0.4)
  grid <- make_efron_grid(theta_hat, sigma)
  prepare_stan_data <- getFromNamespace("prepare_stan_data", "bayesEfron")

  bad_grid <- grid
  bad_grid$B <- bad_grid$B[-1, , drop = FALSE]

  cases <- list(
    tryCatch(
      prepare_stan_data(list(theta_hat = theta_hat, sigma = c(0.2, 0, 0.4)), grid),
      error = identity
    ),
    tryCatch(
      prepare_stan_data(list(theta_hat = theta_hat, sigma = c(0.2, 0.3)), grid),
      error = identity
    ),
    tryCatch(
      prepare_stan_data(list(theta_hat = theta_hat, sigma = sigma), bad_grid),
      error = identity
    ),
    tryCatch(
      prepare_stan_data(
        list(theta_hat = theta_hat, sigma = sigma), grid,
        model_family = "HE"
      ),
      error = identity
    ),
    tryCatch(
      prepare_stan_data(
        list(theta_hat = theta_hat, sigma = sigma), grid,
        group = seq_along(theta_hat)
      ),
      error = identity
    ),
    tryCatch(
      prepare_stan_data(
        list(theta_hat = theta_hat, sigma = sigma), grid,
        rho = 0.1
      ),
      error = identity
    )
  )

  expect_true(all(vapply(cases, inherits, logical(1), "bef_invalid_args")))
  expect_true(all(vapply(cases, inherits, logical(1), "bef_pipeline_error")))
  expect_true(all(vapply(cases, inherits, logical(1), "bef_error")))
})
