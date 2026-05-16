test_that("make_efron_grid returns the locked grid schema and postconditions", {
  theta_hat <- c(-1, 0, 2, 3)
  sigma <- c(0.2, 0.3, 0.4, 0.5)
  theta_true <- c(-2, -0.5, 1, 2)

  grids <- list(
    make_efron_grid(theta_hat, sigma),
    make_efron_grid(
      theta_hat, sigma,
      grid_method = "paper_simulation", theta_true = theta_true
    ),
    make_efron_grid(
      theta_hat, sigma,
      grid_method = "paper_sensitivity", theta_true = theta_true
    ),
    suppressMessages(make_efron_grid(
      theta_hat, sigma,
      grid_method = "kl_target_experimental", kappa = 0.1
    ))
  )

  for (grid in grids) {
    expect_named(
      grid,
      c("grid", "B", "M", "L", "expansion", "kappa", "grid_method",
        "attribution")
    )
    expect_type(grid$grid, "double")
    expect_s3_class(grid$B, "matrix")
    expect_type(grid$M, "integer")
    expect_type(grid$L, "integer")
    expect_true(grid$L >= 51L)
    expect_true(grid$L <= 300L)
    expect_length(grid$grid, grid$L)
    expect_true(all(is.finite(grid$grid)))
    expect_true(all(diff(grid$grid) > 0))
    expect_equal(dim(grid$B), c(grid$L, grid$M))
    expect_true(all(is.finite(grid$B)))
    expect_equal(
      qr(cbind(1, grid$B), tol = sqrt(.Machine$double.eps))$rank,
      grid$M + 1L
    )
    expect_type(grid$attribution$formula, "character")
    expect_type(grid$attribution$source, "character")
    expect_true(nzchar(grid$attribution$formula))
    expect_true(nzchar(grid$attribution$source))
  }
})

test_that("grid methods keep their documented construction formulas", {
  theta_hat <- c(-1, 0, 2, 3)
  sigma <- c(0.2, 0.3, 0.4, 0.5)
  theta_true <- c(-2, -0.5, 1, 2)

  realdata <- make_efron_grid(theta_hat, sigma, L = 101L, expansion = 0.5)
  expect_equal(realdata$grid, seq(-3, 5, length.out = 101L))
  expect_equal(realdata$expansion, 0.5)
  expect_null(realdata$kappa)

  simulation <- make_efron_grid(
    theta_hat, sigma, grid_method = "paper_simulation",
    theta_true = theta_true, L = 101L
  )
  expect_equal(simulation$grid, seq(-2.5, 2.5, length.out = 101L))
  expect_true(is.na(simulation$expansion))
  expect_null(simulation$kappa)

  sensitivity <- make_efron_grid(
    theta_hat, sigma, grid_method = "paper_sensitivity",
    theta_true = theta_true, L = 101L, bound_expansion = 0.25
  )
  expected_sensitivity <- seq(
    min(theta_true) - diff(range(theta_true)) * 0.25,
    max(theta_true) + diff(range(theta_true)) * 0.25,
    length.out = 101L
  )
  expect_equal(sensitivity$grid, expected_sensitivity)
  expect_equal(sensitivity$expansion, 0.25)
  expect_null(sensitivity$kappa)

  kl <- suppressMessages(make_efron_grid(
    theta_hat, sigma, grid_method = "kl_target_experimental",
    kappa = 0.1, expansion = 0.5
  ))
  width <- diff(range(theta_hat))
  expected_L <- as.integer(min(
    max(51, ceiling(width * 2 / (2 * min(sigma) * sqrt(expm1(0.2)))) + 1),
    300
  ))
  expect_equal(kl$L, expected_L)
  expect_equal(kl$expansion, 0.5)
  expect_equal(kl$kappa, 0.1)
})

test_that("kl_target_experimental disclaimer fires once per session state", {
  reset_messages <- getFromNamespace(".bef_reset_session_messages", "bayesEfron")
  reset_messages()
  on.exit(reset_messages(), add = TRUE)

  theta_hat <- c(-1, 0, 2, 3)
  sigma <- c(0.2, 0.3, 0.4, 0.5)
  kl_call <- function() {
    make_efron_grid(
      theta_hat, sigma,
      grid_method = "kl_target_experimental",
      kappa = 0.02,
      expansion = 0.5
    )
  }

  expect_message(
    first <- kl_call(),
    "kl_target_experimental.*homoscedastic"
  )
  expect_silent(second <- kl_call())
  expect_equal(first$grid, second$grid)
  expect_equal(first$kappa, second$kappa)
})

test_that("grid validation failures use bayesEfron typed conditions", {
  err <- tryCatch(
    make_efron_grid(c(1, 1, 1), c(0.2, 0.3, 0.4)),
    error = identity
  )
  expect_s3_class(err, "bef_invalid_args")
  expect_s3_class(err, "bef_grid_error")
  expect_s3_class(err, "bef_error")

  err <- tryCatch(
    make_efron_grid(
      c(-1, 0, 2), c(0.2, 0.3, 0.4),
      grid_method = "paper_simulation"
    ),
    error = identity
  )
  expect_s3_class(err, "bef_err_grid_oracle_required")
  expect_s3_class(err, "bef_grid_error")
  expect_s3_class(err, "bef_error")
})
