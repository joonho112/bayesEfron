# The number of grid points and the warning for a grid that is too coarse or
# too narrow for the data.

resolution_ns <- function(name) getFromNamespace(name, "bayesEfron")

resolution_warning <- function(expr) {
  caught <- NULL
  value <- withCallingHandlers(
    expr,
    bef_warn_grid_resolution = function(w) {
      caught <<- w
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warning = caught)
}

test_that("the three example data sets get 101 points and no warning", {
  for (data in list(krw_discrimination, raudenbush1985)) {
    expect_no_warning(grid <- make_efron_grid(data, grid_method = "paper_realdata"))
    expect_identical(grid$L, 101L)
    expect_true(grid$resolution$L_chosen)
    expect_true(grid$resolution$adequate)
    expect_lte(grid$resolution$L_needed, 101L)
    expect_identical(grid$grid, make_efron_grid(data, L = 101L, grid_method = "paper_realdata")$grid)
  }
  # The simulated sites on the range of 0.3.0, the coarsest of the three:
  # the smallest standard error is 0.98 spacings.
  expect_no_warning(grid <- make_efron_grid(twin_towers$yi, twin_towers$sei, grid_method = "paper_realdata"))
  expect_identical(grid$L, 101L)
  expect_equal(round(grid$resolution$sigma_min_over_spacing, 2), 0.98)
  expect_gte(grid$resolution$scale_over_spacing, 4)
})

test_that("without `L` the grid gets as many points as the standard errors and the spread ask for", {
  y <- krw_discrimination$yi
  s <- krw_discrimination$sei
  base <- make_efron_grid(y, s, grid_method = "paper_realdata")
  width <- diff(base$support)

  # Standard errors a fifth as large: 0.44 spacings at 101 points.
  expect_no_warning(fine <- make_efron_grid(y, s / 5, grid_method = "paper_realdata"))
  expect_gt(fine$L, 101L)
  expect_lte(fine$L, 401L)
  expect_true(fine$resolution$L_chosen)
  expect_true(fine$resolution$adequate)
  expect_identical(fine$L, fine$resolution$L_needed)
  expect_gte(fine$resolution$sigma_min_over_spacing, 0.75 * (1 - 1e-9))
  expect_gte(fine$resolution$scale_over_spacing, 4 * (1 - 1e-9))
  # It is the smallest such number.
  expect_lt(min(s / 5) / (width / (fine$L - 2L)), 0.75)
  expect_identical(fine$support, base$support)

  # The number that is needed is the larger of the two that the conditions
  # give.
  need <- resolution_ns(".bef_resolution_need")(base$support, y, s / 5, "point", NA_real_)
  scale <- resolution_ns(".bef_robust_scale")(y, s / 5)
  expect_identical(
    need$L_needed,
    as.integer(max(ceiling(width * 0.75 / min(s / 5)), ceiling(width * 4 / scale))) + 1L
  )
})

test_that("a standard error of a hundredth of the spacing is reported", {
  # An estimate midway between two grid points whose standard error is the
  # spacing divided by 100.
  base <- make_efron_grid(krw_discrimination, grid_method = "paper_realdata")
  y <- c(krw_discrimination$yi, base$grid[40L] + base$spacing / 2)
  s <- c(krw_discrimination$sei, base$spacing / 100)

  out <- resolution_warning(make_efron_grid(y, s, grid_method = "paper_realdata"))
  expect_s3_class(out$warning, "bef_warn_grid_resolution")
  expect_s3_class(out$warning, "bef_warning")
  expect_identical(out$value$L, 401L)
  expect_false(out$value$resolution$adequate)
  expect_gt(out$value$resolution$L_needed, 401L)
  expect_match(conditionMessage(out$warning), "smallest standard error is 0.04 spacings")
  expect_match(conditionMessage(out$warning), "No `L` up to 401 is enough", fixed = TRUE)
  expect_match(conditionMessage(out$warning), "component = \"normal\"", fixed = TRUE)
  expect_equal(out$warning$sigma_min_over_spacing, out$value$resolution$sigma_min_over_spacing)
  expect_identical(out$warning$L, 401L)

  # With a number of points that was given, the warning names one that is
  # enough when there is one.
  out <- resolution_warning(make_efron_grid(krw_discrimination$yi, krw_discrimination$sei / 5,
                                            L = 101L, grid_method = "paper_realdata"))
  expect_s3_class(out$warning, "bef_warn_grid_resolution")
  expect_identical(out$value$L, 101L)
  expect_false(out$value$resolution$L_chosen)
  needed <- out$value$resolution$L_needed
  expect_match(conditionMessage(out$warning),
               sprintf("`L = %d` is enough, and `L = NULL` chooses it.", needed), fixed = TRUE)
  expect_no_warning(make_efron_grid(krw_discrimination$yi, krw_discrimination$sei / 5, L = needed, grid_method = "paper_realdata"))

  # Normal components do not need the standard errors to exceed the spacing.
  expect_no_warning(
    normal <- make_efron_grid(y, s, component = "normal", grid_method = "paper_realdata")
  )
  expect_identical(normal$L, 101L)
  expect_true(normal$resolution$adequate)
})

test_that("a grid far narrower than the standard errors is reported", {
  # Ten estimates that differ by 1e-12.
  out <- resolution_warning(make_efron_grid(0.2 + 1e-12 * 1:10, rep(0.1, 10), grid_method = "paper_realdata"))
  expect_s3_class(out$warning, "bef_warn_grid_resolution")
  expect_equal(diff(out$value$support), 1.8e-11, tolerance = 1e-3)
  expect_lt(out$value$resolution$width_over_sigma, 1e-9)
  expect_match(conditionMessage(out$warning), "less than the median standard error")
  expect_match(conditionMessage(out$warning), "Check the estimates and the standard errors")

  # The rule "pilot" does not fit on such a range, and reports the same.
  testthat::local_mocked_bindings(
    .bef_pilot_fit = function(...) stop("the pilot fit should not run"),
    .package = "bayesEfron"
  )
  out <- resolution_warning(
    make_efron_grid(0.2 + 1e-12 * 1:10, rep(0.1, 10), grid_method = "pilot", seed = 3L)
  )
  expect_s3_class(out$warning, "bef_warn_grid_resolution")
  expect_false(out$value$pilot$narrowed)
  expect_identical(out$value$pilot$chains, 0L)
  expect_identical(out$value$grid,
                   suppressWarnings(make_efron_grid(0.2 + 1e-12 * 1:10, rep(0.1, 10), grid_method = "paper_realdata"))$grid)

  # A range given by the caller is held to the same standard.
  expect_warning(
    make_efron_grid(krw_discrimination, support = c(0.01, 0.012)),
    class = "bef_warn_grid_resolution"
  )
})

test_that("the spread of the effects sets the number of points when it is small next to the range", {
  # One estimate far from the others doubles the range of the grid and more.
  set.seed(8)
  y <- c(stats::rnorm(60, 0, 0.1), 6)
  s <- c(rep(0.2, 60), 3)
  expect_no_warning(grid <- make_efron_grid(y, s, grid_method = "paper_realdata"))
  expect_gt(grid$L, 101L)
  expect_gte(grid$resolution$scale_over_spacing, 4 * (1 - 1e-9))
  expect_true(grid$resolution$adequate)

  out <- resolution_warning(make_efron_grid(y, s, L = 101L, grid_method = "paper_realdata"))
  expect_match(conditionMessage(out$warning), "spread of the effects is about")
  expect_false(grepl("smallest standard error", conditionMessage(out$warning)))
})

test_that("normal components are held to their width next to the spread of the effects", {
  y <- krw_discrimination$yi
  s <- krw_discrimination$sei
  expect_no_warning(grid <- make_efron_grid(y, s, component = "normal", grid_method = "paper_realdata"))
  expect_identical(grid$L, 101L)
  # tau_c over the scale: the width over the scale in spacings.
  expect_lte(grid$tau_c / resolution_ns(".bef_robust_scale")(y, s), 1 / 6)

  # Wider components need a finer grid.
  expect_no_warning(wide <- make_efron_grid(y, s, component = "normal", width = 2, grid_method = "paper_realdata"))
  expect_gt(wide$L, 101L)
  expect_lte(wide$tau_c / resolution_ns(".bef_robust_scale")(y, s), (1 / 6) * (1 + 1e-9))

  out <- resolution_warning(make_efron_grid(y, s, component = "normal", width = 2, L = 101L, grid_method = "paper_realdata"))
  expect_s3_class(out$warning, "bef_warn_grid_resolution")
  expect_match(conditionMessage(out$warning), "A normal component is 0.25 times as wide")
  expect_match(conditionMessage(out$warning), sprintf("`L = %d` is enough", wide$L), fixed = TRUE)
})

test_that("the number of grid points must be between 51 and 401", {
  y <- krw_discrimination$yi
  s <- krw_discrimination$sei
  for (bad in list(50L, 402L, 101.5, c(101L, 201L), "101", NA_integer_)) {
    expect_error(make_efron_grid(y, s, L = bad, grid_method = "paper_realdata"), class = "bef_grid_error")
  }
  expect_identical(suppressWarnings(make_efron_grid(y, s, L = 51L, grid_method = "paper_realdata"))$L, 51L)
  expect_identical(make_efron_grid(y, s, L = 401L, grid_method = "paper_realdata")$L, 401L)
  expect_false(make_efron_grid(y, s, L = 401L, grid_method = "paper_realdata")$resolution$L_chosen)
})

test_that("a fit chooses the number of grid points in the same way", {
  prepare <- function(...) {
    resolution_ns(".bef_fit_prepare_context")(
      call = quote(bayes_efron_fit(theta_hat, sigma)), check_installed = FALSE, ...
    )
  }
  y <- krw_discrimination$yi
  s <- krw_discrimination$sei
  expect_no_warning(context <- prepare(theta_hat = y, sigma = s))
  expect_identical(context$stan_data$L, 101L)
  expect_no_warning(context <- prepare(theta_hat = y, sigma = s / 5))
  expect_identical(context$stan_data$L, make_efron_grid(y, s / 5, grid_method = "paper_realdata")$L)
  expect_gt(context$stan_data$L, 101L)
  expect_warning(prepare(theta_hat = y, sigma = s / 5, L = 101L),
                 class = "bef_warn_grid_resolution")
  expect_identical(suppressWarnings(prepare(theta_hat = y, sigma = s / 5, L = 101L))$stan_data$L, 101L)
})

test_that("a grid that is too coarse prints as such", {
  lines <- format(suppressWarnings(
    make_efron_grid(krw_discrimination$yi, krw_discrimination$sei / 5, L = 101L, grid_method = "paper_realdata")
  ))
  expect_match(lines[[6L]], "(too coarse)", fixed = TRUE)
  expect_false(grepl("too coarse", format(make_efron_grid(krw_discrimination, grid_method = "paper_realdata"))[[6L]]))
})
