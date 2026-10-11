# as_normalmix(): the posterior mean of g as a mixture that ebnm takes, and
# the posterior distributions that ebnm computes from it.

normalmix_ns <- function(name) getFromNamespace(name, "bayesEfron")

# The posterior means and standard deviations of the sites under one g, from
# the engine of the package.
normalmix_plug_in <- function(fit) {
  engine <- normalmix_ns(".bef_engine")(fit)
  g <- colMeans(exp(engine$log_g))
  one <- normalmix_ns(".bef_engine_from_parts")(
    log_g = matrix(log(g / sum(g)), nrow = 1L), grid = engine$grid,
    theta_hat = engine$theta_hat, sigma = engine$sigma, tau_c = engine$tau_c
  )
  moments <- normalmix_ns(".bef_site_moments")(one)
  list(mean = moments$mean, sd = sqrt(moments$var))
}

test_that("the object has the form of a normal mixture without ashr", {
  fit <- plot_fit(12L)
  g <- as_normalmix(fit)
  expect_s3_class(g, "normalmix", exact = TRUE)
  expect_named(unclass(g), c("pi", "mean", "sd"))
  grid <- fit$metadata$data_list$grid
  expect_identical(g$mean, as.numeric(grid))
  expect_equal(g$pi, g_summary(fit, type = "mass")$mean / sum(g_summary(fit, type = "mass")$mean))
  expect_equal(sum(g$pi), 1)
  expect_true(all(g$pi >= 0))
  expect_identical(g$sd, rep(0, length(grid)))
  expect_identical(attr(g, "row.names"), seq_along(grid))

  normal <- plot_fit(12L, component = "normal")
  tau_c <- normal$metadata$data_list$tau_c
  expect_gt(tau_c, 0)
  expect_identical(as_normalmix(normal)$sd, rep(tau_c, length(grid)))

  expect_error(as_normalmix(list()), class = "bef_invalid_fit")
  expect_error(as_normalmix(fit, level = 0.9), "no argument", class = "bef_invalid_args")
})

test_that("ebnm computes the plug-in posterior distributions from it", {
  skip_if_not_installed("ebnm")
  for (component in c("point", "normal")) {
    fit <- plot_fit(150L, component = component)
    data <- fit$metadata$data_list
    plug_in <- ebnm::ebnm_npmle(unname(data$theta_hat), data$sigma, g_init = as_normalmix(fit),
                                fix_g = TRUE)
    expected <- normalmix_plug_in(fit)
    expect_equal(unname(stats::coef(plug_in)), expected$mean, tolerance = 1e-8)
    expect_equal(stats::fitted(plug_in)$sd, expected$sd, tolerance = 1e-8)
    # The prior that ebnm reports, which it builds with ashr, is the object
    # that it was given: the same class, elements and attributes.
    expect_equal(plug_in$fitted_g, as_normalmix(fit))
  }
})

test_that("the fit differs from the plug-in by the uncertainty about g", {
  skip_if_not_installed("ebnm")
  fit <- plot_fit(150L)
  data <- fit$metadata$data_list
  plug_in <- ebnm::ebnm_npmle(unname(data$theta_hat), data$sigma, g_init = as_normalmix(fit),
                              fix_g = TRUE)
  sites <- site_summary(fit)
  # The means are close; the standard deviations of the fit are not smaller
  # than what one g gives on average.
  expect_lt(max(abs(sites$mean - stats::coef(plug_in)) / data$sigma), 0.25)
  expect_gt(mean(sites$sd), mean(stats::fitted(plug_in)$sd) * 0.999)
  expect_true(all(sites$sd_about_g >= 0))
})
