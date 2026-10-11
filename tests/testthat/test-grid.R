test_that("make_efron_grid returns a grid with the documented elements", {
  theta_hat <- c(-1, 0, 2, 3)
  sigma <- c(0.2, 0.3, 0.4, 0.5)
  theta_true <- c(-2, -0.5, 1, 2)

  grids <- list(
    make_efron_grid(theta_hat, sigma, grid_method = "paper_realdata"),
    make_efron_grid(
      theta_hat, sigma,
      grid_method = "paper_simulation", theta_true = theta_true
    ),
    make_efron_grid(
      theta_hat, sigma,
      grid_method = "paper_sensitivity", theta_true = theta_true
    )
  )

  for (grid in grids) {
    expect_s3_class(grid, "bef_grid")
    expect_named(
      grid,
      c("grid", "B", "M", "L", "expansion", "grid_method",
        "attribution", "support", "spacing", "component", "width", "tau_c",
        "knots", "resolution", "data")
    )
    expect_identical(grid$data, list(theta_hat = theta_hat, sigma = sigma))
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

  realdata <- make_efron_grid(theta_hat, sigma, L = 101L, expansion = 0.5, grid_method = "paper_realdata")
  expect_equal(realdata$grid, seq(-3, 5, length.out = 101L))
  expect_equal(realdata$expansion, 0.5)

  simulation <- make_efron_grid(
    theta_hat, sigma, grid_method = "paper_simulation",
    theta_true = theta_true, L = 101L
  )
  expect_equal(simulation$grid, seq(-2.5, 2.5, length.out = 101L))
  expect_true(is.na(simulation$expansion))

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
})

test_that("the rule that version 0.4 removed names what takes its place", {
  theta_hat <- c(-1, 0, 2, 3, 4)
  sigma <- c(0.2, 0.3, 0.4, 0.5, 0.3)
  errors <- list(
    tryCatch(make_efron_grid(theta_hat, sigma, grid_method = "kl_target_experimental"),
             error = identity),
    tryCatch(make_efron_grid(theta_hat, sigma, grid_method = "paper_realdata", kappa = 0.1),
             error = identity),
    tryCatch(make_efron_grid(theta_hat, sigma, grid_method = "kl_target_experimental",
                             kappa = 0.1),
             error = identity),
    tryCatch(bayes_efron_fit(theta_hat, sigma, grid_method = "kl_target_experimental"),
             error = identity)
  )
  for (err in errors) {
    expect_s3_class(err, "bef_error")
    expect_s3_class(err, "bef_invalid_args")
    expect_match(conditionMessage(err), "removed in version 0.4", fixed = TRUE)
    expect_match(conditionMessage(err), "Leave `L` out", fixed = TRUE)
  }
  expect_false("kl_target_experimental" %in%
                 getFromNamespace(".bef_grid_methods", "bayesEfron")())
  expect_false("kappa" %in% names(formals(make_efron_grid)))
  expect_false("kappa" %in% names(make_efron_grid(theta_hat, sigma,
                                                  grid_method = "paper_realdata")))
  # Without `L` the grid is as fine as the smallest standard error asks for,
  # which is what the rule was for.
  fine <- make_efron_grid(theta_hat, c(0.02, 0.3, 0.4, 0.5, 0.3),
                          grid_method = "paper_realdata")
  expect_gt(fine$L, 101L)
  expect_gte(fine$resolution$sigma_min_over_spacing, 0.75 * (1 - 1e-9))
})

test_that("make_efron_grid gives a bef_grid_error for estimates that are all equal and for a missing theta_true", {
  err <- tryCatch(
    make_efron_grid(c(1, 1, 1), c(0.2, 0.3, 0.4), grid_method = "paper_realdata"),
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

# The three grid methods that follow Lee and Sui (2025), rebuilt from stored
# inputs and compared with stored grids. The tolerance leaves room for
# last-bit differences in splines::ns() between platforms.
grid_fixture <- function(name) {
  readRDS(test_path("_fixtures", "grid", name))
}

expect_grid_equal <- function(actual, expected) {
  expect_equal(actual$grid, expected$grid, tolerance = 1e-12)
  expect_equal(
    unclass(actual$B), unclass(expected$B),
    tolerance = 1e-12, ignore_attr = TRUE
  )
  expect_identical(actual$L, expected$L)
  expect_identical(actual$M, expected$M)
  expect_identical(actual$grid_method, expected$grid_method)
}

expect_grid_differs <- function(actual, expected) {
  expect_false(isTRUE(all.equal(actual$grid, expected$grid, tolerance = 1e-12)))
}

test_that("paper_realdata reproduces the stored grid", {
  expected <- grid_fixture("paper_realdata_grid.rds")
  inputs <- grid_fixture("paper_realdata_inputs.rds")$data

  actual <- make_efron_grid(
    theta_hat = inputs$theta_hat,
    sigma = inputs$sigma,
    L = expected$L,
    M = expected$M,
    expansion = expected$expansion,
    grid_method = "paper_realdata"
  )
  expect_grid_equal(actual, expected)

  narrower <- make_efron_grid(
    theta_hat = inputs$theta_hat,
    sigma = inputs$sigma,
    L = expected$L,
    M = expected$M,
    expansion = 0.25,
    grid_method = "paper_realdata"
  )
  expect_grid_differs(narrower, expected)
})

test_that("paper_simulation reproduces the stored grid and needs theta_true", {
  expected <- grid_fixture("paper_simulation_grid.rds")
  inputs <- grid_fixture("paper_simulation_inputs.rds")$data

  actual <- make_efron_grid(
    theta_hat = inputs$theta_hat,
    sigma = inputs$sigma,
    L = expected$L,
    M = expected$M,
    grid_method = "paper_simulation",
    theta_true = inputs$theta_true
  )
  expect_grid_equal(actual, expected)

  expect_error(
    make_efron_grid(
      theta_hat = inputs$theta_hat,
      sigma = inputs$sigma,
      L = expected$L,
      M = expected$M,
      grid_method = "paper_simulation"
    ),
    class = "bef_err_grid_oracle_required"
  )

  from_estimates <- make_efron_grid(
    theta_hat = inputs$theta_hat,
    sigma = inputs$sigma,
    L = expected$L,
    M = expected$M,
    grid_method = "paper_simulation",
    theta_true = inputs$theta_hat
  )
  expect_grid_differs(from_estimates, expected)
})

test_that("paper_sensitivity reproduces the stored grid and defaults to bound_expansion = 0.5", {
  expected <- grid_fixture("paper_sensitivity_grid.rds")
  inputs <- grid_fixture("paper_sensitivity_inputs.rds")$data
  expect_identical(inputs$bound_expansion, 0.5)
  expect_identical(expected$expansion, inputs$bound_expansion)

  actual <- make_efron_grid(
    theta_hat = inputs$theta_hat,
    sigma = inputs$sigma,
    L = expected$L,
    M = expected$M,
    grid_method = "paper_sensitivity",
    theta_true = inputs$theta_true,
    bound_expansion = inputs$bound_expansion
  )
  expect_grid_equal(actual, expected)

  default_bound <- make_efron_grid(
    theta_hat = inputs$theta_hat,
    sigma = inputs$sigma,
    L = expected$L,
    M = expected$M,
    grid_method = "paper_sensitivity",
    theta_true = inputs$theta_true
  )
  expect_identical(default_bound, actual)

  expect_error(
    make_efron_grid(
      theta_hat = inputs$theta_hat,
      sigma = inputs$sigma,
      L = expected$L,
      M = expected$M,
      grid_method = "paper_sensitivity",
      bound_expansion = expected$expansion
    ),
    class = "bef_err_grid_oracle_required"
  )

  narrower <- make_efron_grid(
    theta_hat = inputs$theta_hat,
    sigma = inputs$sigma,
    L = expected$L,
    M = expected$M,
    grid_method = "paper_sensitivity",
    theta_true = inputs$theta_true,
    bound_expansion = 0.25
  )
  expect_grid_differs(narrower, expected)

  fixed_padding <- make_efron_grid(
    theta_hat = inputs$theta_hat,
    sigma = inputs$sigma,
    L = expected$L,
    M = expected$M,
    grid_method = "paper_simulation",
    theta_true = inputs$theta_true
  )
  expect_grid_differs(fixed_padding, expected)
})
