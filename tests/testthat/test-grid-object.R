# The grid object of version 0.4: the rules of the paper build what they
# built in 0.3.0, and the support, the component and the knots are parts of
# the object that can be set and read one by one.

grid_ref <- function() {
  readRDS(testthat::test_path("_fixtures", "grid-v030.rds"))
}

test_that("the rules of the paper give the grid and the basis of version 0.3.0", {
  ref <- grid_ref()
  d <- ref$data
  same <- function(grid, stored) {
    expect_identical(grid$grid, stored$grid)
    expect_identical(matrix(as.numeric(grid$B), nrow = grid$L), stored$B)
    expect_identical(grid$L, stored$L)
    expect_identical(grid$M, stored$M)
    expect_identical(grid$grid_method, stored$grid_method)
    expect_identical(grid$expansion, stored$expansion)
  }
  same(make_efron_grid(d$theta_hat, d$sigma, grid_method = "paper_realdata"), ref$realdata)
  same(make_efron_grid(d$theta_hat, d$sigma, L = 201L, expansion = 1, M = 8L,
                       grid_method = "paper_realdata"), ref$realdata_wide)
  same(make_efron_grid(d$theta_hat, d$sigma, grid_method = "paper_simulation",
                       theta_true = d$theta_true), ref$simulation)
  same(make_efron_grid(d$theta_hat, d$sigma, L = 151L, M = 5L,
                       grid_method = "paper_sensitivity", theta_true = d$theta_true,
                       bound_expansion = 0.25), ref$sensitivity)
  same(make_efron_grid(krw_discrimination$yi, krw_discrimination$sei,
                       grid_method = "paper_realdata"), ref$krw)
})

test_that("the basis of rule \"grid\" does not depend on where the grid lies", {
  d <- grid_ref()$data
  a <- make_efron_grid(d$theta_hat, d$sigma, L = 101L, grid_method = "paper_realdata")
  # The second grid is far too coarse for these data, which is not the
  # matter here.
  b <- suppressWarnings(
    make_efron_grid(d$theta_hat, d$sigma, support = c(100, 250), L = 101L),
    classes = "bef_warn_grid_resolution"
  )
  expect_equal(matrix(as.numeric(a$B), nrow = a$L), matrix(as.numeric(b$B), nrow = b$L),
               tolerance = 1e-10)
})

test_that("support sets the ends and takes the place of a rule", {
  d <- grid_ref()$data
  grid <- make_efron_grid(d$theta_hat, d$sigma, support = c(-2, 3), L = 51L)
  expect_identical(grid$support, c(-2, 3))
  expect_identical(range(grid$grid), c(-2, 3))
  expect_equal(grid$spacing, 0.1)
  expect_identical(grid$grid_method, "support")
  expect_true(is.na(grid$expansion))
  expect_identical(grid$attribution$formula, "ends given in `support`")

  expect_error(
    make_efron_grid(d$theta_hat, d$sigma, support = c(-2, 3), grid_method = "paper_realdata"),
    "not both", class = "bef_invalid_args"
  )
  for (bad in list(c(1, 1), c(2, 1), 1, c(0, NA), c(-Inf, 1), "a")) {
    expect_error(make_efron_grid(d$theta_hat, d$sigma, support = bad),
                 "two finite numbers", class = "bef_invalid_args")
  }
})

test_that("the component and its width are parts of the grid", {
  d <- grid_ref()$data
  point <- make_efron_grid(d$theta_hat, d$sigma, grid_method = "paper_realdata")
  expect_identical(point$component, "point")
  expect_identical(point$tau_c, 0)
  expect_true(is.na(point$width))

  normal <- make_efron_grid(d$theta_hat, d$sigma, component = "normal", grid_method = "paper_realdata")
  expect_identical(normal$component, "normal")
  expect_equal(normal$width, 2 / 3)
  expect_equal(normal$tau_c, 2 / 3 * normal$spacing)
  expect_identical(normal$grid, point$grid)

  half <- make_efron_grid(d$theta_hat, d$sigma, component = "normal", width = 0.5, grid_method = "paper_realdata")
  expect_equal(half$tau_c, 0.5 * half$spacing)

  expect_error(make_efron_grid(d$theta_hat, d$sigma, component = "t", grid_method = "paper_realdata"),
               "`component` must be", class = "bef_invalid_args")
  expect_error(make_efron_grid(d$theta_hat, d$sigma, component = "normal", width = 0, grid_method = "paper_realdata"),
               "`width`", class = "bef_invalid_args")
  expect_error(make_efron_grid(d$theta_hat, d$sigma, knots = "quantile", grid_method = "paper_realdata"),
               "`knots` must be", class = "bef_invalid_args")
})

test_that("the knots and the resolution are recorded", {
  grid <- make_efron_grid(krw_discrimination$yi, krw_discrimination$sei, grid_method = "paper_realdata")
  expect_identical(grid$knots$rule, "grid")
  expect_length(grid$knots$interior, grid$M - 1L)
  expect_equal(grid$knots$boundary, grid$support)
  expect_true(all(grid$knots$interior > grid$support[1L] &
                    grid$knots$interior < grid$support[2L]))
  expect_equal(diff(grid$knots$interior), rep(diff(grid$support) / grid$M, grid$M - 2L),
               tolerance = 1e-8)

  r <- grid$resolution
  expect_named(r, c("sigma_min_over_spacing", "scale_over_spacing", "share_outside_estimates",
                    "width_over_sigma", "L_chosen", "L_needed", "adequate"))
  expect_true(r$L_chosen)
  expect_true(r$adequate)
  expect_equal(r$width_over_sigma, diff(grid$support) / stats::median(krw_discrimination$sei))
  expect_equal(r$sigma_min_over_spacing, min(krw_discrimination$sei) / grid$spacing)
  expect_equal(r$share_outside_estimates,
               mean(grid$grid < min(krw_discrimination$yi) | grid$grid > max(krw_discrimination$yi)))
  # For the callback data, 2.21 and 7.89 spacings.
  expect_equal(round(r$sigma_min_over_spacing, 2), 2.21)
  expect_equal(round(r$scale_over_spacing, 2), 7.89)
})

test_that("for every number of spline functions the basis has full rank and its knots lie inside the grid", {
  for (M in 3:10) {
    for (L in c(51L, 101L, 401L)) {
      grid <- suppressWarnings(
        make_efron_grid(krw_discrimination, M = M, L = L, grid_method = "paper_realdata"),
        classes = "bef_warn_grid_resolution"
      )
      expect_identical(dim(grid$B), c(L, M))
      expect_identical(qr(cbind(1, grid$B))$rank, M + 1L)
      expect_identical(grid$knots$rule, "grid")
      expect_length(grid$knots$interior, M - 1L)
      expect_true(all(grid$knots$interior > grid$support[1L] &
                        grid$knots$interior < grid$support[2L]))
      expect_equal(grid$knots$boundary, grid$support)
      # The knots of rule "grid": quantiles of the grid points, which are
      # equally spaced, so the knots are too.
      expect_equal(diff(c(grid$support[1L], grid$knots$interior, grid$support[2L])),
                   rep(diff(grid$support) / M, M), tolerance = 1e-8)
    }
  }
  for (bad in list(2L, 11L, 6.5, "6", NA_integer_)) {
    expect_error(make_efron_grid(krw_discrimination, M = bad, grid_method = "paper_realdata"), class = "bef_grid_error")
  }
  expect_identical(make_efron_grid(krw_discrimination, grid_method = "paper_realdata")$M, 6L)
})

test_that("the robust scale falls back on the sampling variance", {
  expect_equal(.bef_robust_scale(c(0, 0.01, -0.01, 0.02), rep(1, 4)), 1)
  y <- c(-3, -1, 0, 1, 3)
  expect_equal(.bef_robust_scale(y, rep(0.1, 5)), sqrt(stats::mad(y)^2 - 0.01))
})

test_that("settings must be named in full and a grid of points is not taken", {
  d <- grid_ref()$data
  expect_error(make_efron_grid(d$theta_hat, d$sigma, 101L, grid_method = "paper_realdata"),
               "must be given by name", class = "bef_invalid_args")
  expect_error(make_efron_grid(d$theta_hat, d$sigma, grid = seq(-3, 3, 0.1), grid_method = "paper_realdata"),
               "`grid` is not an argument", class = "bef_invalid_args")
  # `grid_m` would match `grid_method` by its beginning if the settings could
  # be abbreviated.
  expect_error(make_efron_grid(d$theta_hat, d$sigma, grid_m = "paper_realdata", grid_method = "paper_realdata"),
               "Unused arguments: `grid_m`", class = "bef_invalid_args")
})

test_that("a grid prints in a few lines", {
  grid <- make_efron_grid(krw_discrimination$yi, krw_discrimination$sei, grid_method = "paper_realdata")
  lines <- format(grid)
  expect_lte(length(lines), 10L)
  expect_match(lines[[1L]], "^<bef_grid> 101 points from ")
  expect_true(any(grepl("point masses", lines, fixed = TRUE)))
  expect_output(shown <- withVisible(print(grid)), "<bef_grid> 101 points")
  expect_false(shown$visible)
  expect_identical(shown$value, grid)

  normal <- make_efron_grid(krw_discrimination$yi, krw_discrimination$sei,
                            support = c(-0.05, 0.15), component = "normal")
  expect_true(any(grepl("ends given in `support`", format(normal), fixed = TRUE)))
  expect_true(any(grepl("normal, standard deviation", format(normal), fixed = TRUE)))
})
