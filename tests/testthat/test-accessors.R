# g_summary(), site_summary() and exceedance_prob().

accessor_fit <- function() {
  readRDS(test_path("_fixtures", "short_fit.rds"))
}

test_that("g_summary() returns means inside their intervals for every type", {
  fit <- accessor_fit()
  grid <- fit$metadata$data_list$grid
  spacing <- diff(grid)[1L]

  mass <- g_summary(fit, type = "mass")
  expect_named(mass, c("value", "mean", "lower", "upper"))
  expect_identical(attr(mass, "level"), 0.9)
  expect_equal(mass$value, grid)
  expect_equal(sum(mass$mean), 1, tolerance = 1e-8)
  expect_true(all(mass$lower <= mass$mean + 1e-12 & mass$mean <= mass$upper + 1e-12))

  density <- g_summary(fit)
  expect_equal(sum(density$mean * spacing), 1, tolerance = 1e-8)
  expect_equal(density$mean, mass$mean / spacing)
  # At other values the density is that of the cell that holds the value,
  # and zero beyond the cells.
  off <- g_summary(fit, at = c(grid[3L] + 0.2 * spacing, grid[1L] - spacing))
  expect_equal(off$mean, c(density$mean[3L], 0))

  cdf <- g_summary(fit, type = "cdf", at = c(grid, max(grid) + spacing))
  expect_true(all(diff(cdf$mean) >= -1e-12))
  expect_equal(cdf$mean[nrow(cdf)], 1, tolerance = 1e-8)
  step <- g_summary(fit, type = "cdf", reading = "step")
  expect_equal(step$mean, cumsum(mass$mean), tolerance = 1e-8)

  quantiles <- g_summary(fit, type = "quantile")
  expect_identical(quantiles$prob, c(0.05, 0.25, 0.5, 0.75, 0.95))
  expect_true(all(diff(quantiles$mean) > 0))
  expect_true(all(quantiles$lower <= quantiles$mean & quantiles$mean <= quantiles$upper))

  wider <- g_summary(fit, type = "mass", level = 0.99)
  expect_true(all(wider$lower <= mass$lower & wider$upper >= mass$upper))
})

test_that("the probability of an interval is the difference of the distribution function", {
  fit <- accessor_fit()
  engine <- .bef_engine(fit)
  for (reading in c("interpolated", "step")) {
    cdf <- .bef_g_cdf_draws(engine, exp(engine$log_g), c(0.1, 0.45), reading)
    prob <- g_summary(fit, type = "prob", lower = 0.1, upper = 0.45, reading = reading)
    expect_named(prob, c("from", "to", "mean", "lower", "upper"))
    expect_equal(prob$mean, mean(cdf[, 2L] - cdf[, 1L]))
    above <- g_summary(fit, type = "prob", lower = 0.1, reading = reading)
    below <- g_summary(fit, type = "prob", upper = 0.1, reading = reading)
    expect_equal(above$mean + below$mean, 1, tolerance = 1e-10)
    expect_identical(above$to, Inf)
  }
  expect_error(g_summary(fit, type = "prob", lower = 1, upper = 0), class = "bef_invalid_args")
})

test_that("scale changes the units of effects and densities and leaves probabilities", {
  fit <- accessor_fit()
  expect_equal(g_summary(fit, scale = 100)$value, 100 * g_summary(fit)$value)
  expect_equal(g_summary(fit, scale = 100)$mean, g_summary(fit)$mean / 100)
  expect_equal(g_summary(fit, type = "cdf", scale = 100)$mean, g_summary(fit, type = "cdf")$mean)
  expect_equal(g_summary(fit, type = "quantile", scale = 100)$mean,
               100 * g_summary(fit, type = "quantile")$mean)
  expect_equal(g_summary(fit, type = "prob", lower = 0, scale = 100)$mean,
               g_summary(fit, type = "prob", lower = 0)$mean)
  expect_error(g_summary(fit, scale = 0), class = "bef_invalid_args")
  expect_error(g_summary(fit, probs = c(0, 0.5), type = "quantile"), class = "bef_invalid_args")
  expect_error(g_summary(fit, type = "curve"), class = "bef_invalid_args")
})

test_that("the moments of g are those that confint() returns", {
  fit <- accessor_fit()
  moments <- g_summary(fit, type = "moments", level = 0.8)
  intervals <- confint(fit, type = "g", level = 0.8)
  expect_identical(moments$parameter, c("mean_g", "var_g", "sd_g"))
  expect_equal(moments$mean, intervals$mean)
  expect_equal(moments$lower, intervals$lower)
  expect_equal(moments$upper, intervals$upper)
  scaled <- g_summary(fit, type = "moments", scale = 10)
  expect_equal(scaled$mean, g_summary(fit, type = "moments")$mean * c(10, 100, 10))
})

test_that("for the KRW fit the share of firms with a positive gap is 0.87", {
  # The mass of g above zero, with the probabilities left on the grid
  # points: the value that the documentation of version 0.3.0 computed by
  # hand from the draws. With the fit of that version the interval was 0.799
  # to 0.936; the stored fit of this version has other draws.
  share <- g_summary(bayesEfron::krw_fit, type = "prob", lower = 0, reading = "step")
  expect_equal(round(c(share$mean, share$lower, share$upper), 3), c(0.870, 0.801, 0.933))
  interpolated <- g_summary(bayesEfron::krw_fit, type = "prob", lower = 0)
  expect_lt(abs(interpolated$mean - share$mean), 0.01)
})

test_that("site_summary() extends the site table and its parts add up", {
  fit <- accessor_fit()
  table <- as.data.frame(fit)
  sites <- site_summary(fit)
  expect_named(sites, c("site", "label", "estimate", "se", "mean", "sd", "lower", "upper",
                        "mode", "sd_given_g", "sd_about_g"))
  for (column in c("site", "label", "estimate", "se", "mean", "sd", "lower", "upper", "mode")) {
    expect_equal(sites[[column]], table[[column]])
  }
  expect_equal(sites$sd^2, sites$sd_given_g^2 + sites$sd_about_g^2)
  expect_identical(attr(sites, "level"), 0.9)

  full <- site_summary(fit, level = 0.95, inner_level = 0.5, threshold = 0.1)
  expect_true(all(c("inner_lower", "inner_upper", "prob_above") %in% names(full)))
  expect_true(all(full$lower <= full$inner_lower & full$inner_upper <= full$upper))
  expect_true(all(full$lower <= sites$lower & full$upper >= sites$upper))
  expect_equal(full$prob_above, unname(exceedance_prob(fit, threshold = 0.1)))
  expect_identical(attr(full, "threshold"), 0.1)

  scaled <- site_summary(fit, scale = 100, threshold = 0.1)
  expect_equal(scaled$mean, 100 * sites$mean)
  expect_equal(scaled$sd_about_g, 100 * sites$sd_about_g)
  expect_equal(scaled$prob_above, full$prob_above)

  picked <- site_summary(fit, sites = c(4, 2))
  expect_identical(picked$site, c(4L, 2L))
  expect_equal(picked$mean, sites$mean[c(4L, 2L)])
  expect_identical(site_summary(fit, sites = "3")$site, 3L)
  expect_error(site_summary(fit, sites = 9), "`sites`", class = "bef_invalid_args")
  expect_error(site_summary(fit, threshold = c(0, 1)), class = "bef_invalid_args")
})

test_that("exceedance probabilities fall as the threshold rises and the two directions sum to one", {
  fit <- accessor_fit()
  low <- exceedance_prob(fit, threshold = 0)
  high <- exceedance_prob(fit, threshold = 0.3)
  expect_named(low, as.character(1:5))
  expect_true(all(high <= low))
  expect_true(all(low >= 0 & low <= 1))
  for (reading in c("interpolated", "step")) {
    expect_equal(
      exceedance_prob(fit, 0.2, reading = reading) +
        exceedance_prob(fit, 0.2, direction = "below", reading = reading),
      stats::setNames(rep(1, 5L), 1:5)
    )
  }
  # As steps, the probability is the mass on the grid points above the
  # threshold; a grid point equal to it is not above.
  engine <- .bef_engine(fit)
  weights <- .bef_site_weights(engine)
  point <- engine$grid[30L]
  expect_equal(unname(exceedance_prob(fit, point, reading = "step")),
               rowSums(weights[, engine$grid > point, drop = FALSE]))
  expect_equal(unname(exceedance_prob(fit, 0, sites = c(5, 1))), unname(low[c(5L, 1L)]))
})

test_that("the stored firm probabilities are reproduced by exceedance_prob()", {
  stored <- readRDS(system.file("examples", "lee_sui_examples.rds", package = "bayesEfron"))
  expected <- stored$cases$krw$positive_probability
  skip_if(is.null(expected), "The stored comparison holds no site probabilities.")
  # The stored values leave the probabilities on the grid points.
  expect_equal(unname(exceedance_prob(bayesEfron::krw_fit, 0, reading = "step")),
               as.numeric(expected), tolerance = 1e-8)
})
