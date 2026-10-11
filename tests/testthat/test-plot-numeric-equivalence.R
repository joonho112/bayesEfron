test_that("the base caterpillar plot has sorted means and nested intervals, and drawing does not change its data", {
  fit <- plot_test_fit()
  payload <- plot_test_ns(".bef_plot_payload_bef_fit_re")(
    fit,
    type = "caterpillar",
    level = 0.8,
    sort_by = "mean"
  )

  expect_equal(payload$type, "caterpillar")
  expect_equal(payload$reference, fit$metadata$mean_g_summary$mean)
  expect_named(
    payload$data,
    c(
      "site", "lower", "upper", "point", "inner_lower", "inner_upper",
      "sd", "sigma", "position"
    )
  )
  expect_equal(payload$data$point, sort(payload$data$point))
  expect_true(all(payload$data$lower <= payload$data$inner_lower))
  expect_true(all(payload$data$inner_lower <= payload$data$inner_upper))
  expect_true(all(payload$data$inner_upper <= payload$data$upper))

  before <- payload
  plot_test_with_pdf_device({
    expect_null(plot_test_ns(".bef_plot_base")(payload))
  })
  expect_identical(payload, before)

})

test_that("the base plots and the ggplot2 plots draw the same numbers", {
  fit <- bayesEfron::raudenbush_fit
  payload <- plot_test_ns(".bef_plot_payload_bef_fit_re")
  sites <- payload(fit, type = "caterpillar", level = 0.8, sort_by = "sigma")$data
  drawn <- plot_sites(fit, level = 0.8, sort_by = "sigma")
  # The same sites, intervals and means; plot_sites() puts the most precise
  # estimate at the top and the base plot at the bottom.
  expect_identical(rev(drawn$data$site), sites$site)
  expect_equal(rev(drawn$data$lower), sites$lower)
  expect_equal(rev(drawn$data$upper), sites$upper)
  expect_equal(rev(drawn$data$mean), sites$point)
  expect_equal(rev(drawn$data$inner_lower), sites$inner_lower)
  expect_equal(plot_layer_data(drawn, 1L)$xintercept,
               payload(fit, type = "caterpillar", level = 0.8, sort_by = "sigma")$reference)

  mass <- payload(fit, type = "g", level = 0.8, sort_by = "mean")$data
  drawn <- plot_g(fit, type = "mass", level = 0.8, xlim = range(mass$grid))
  expect_equal(plot_layer_data(drawn, 2L)$y, mass$point)
  expect_equal(plot_layer_data(drawn, 1L)$ymin, mass$lower)
  expect_equal(plot_layer_data(drawn, 1L)$ymax, mass$upper)
})

test_that("for a fit with draws of g, the base plot of g has a mean inside an interval at each grid point", {
  fit <- plot_test_fit(L = 13L)
  payload <- plot_test_ns(".bef_plot_payload_bef_fit_re")(
    fit,
    type = "g",
    level = 0.8,
    sort_by = "mean"
  )

  expect_equal(payload$type, "g")
  expect_named(payload$data, c("kind", "grid", "lower", "upper", "point"))
  expect_equal(payload$data$kind, rep("density", fit$metadata$data_list$L))
  expect_equal(payload$data$grid, fit$metadata$data_list$grid)
  expect_true(all(payload$data$lower <= payload$data$point))
  expect_true(all(payload$data$point <= payload$data$upper))
})

test_that("type = \"g\" draws the density band for a fit made with the default settings", {
  # The example fit was run with store_grid_quantities = FALSE, so it carries
  # log_g but not g.
  fit <- readRDS(test_path("_fixtures", "short_fit.rds"))
  expect_false(any(grepl("^g\\[", posterior::variables(fit$draws))))

  payload <- plot_test_ns(".bef_plot_payload_bef_fit_re")(
    fit,
    type = "g",
    level = 0.9,
    sort_by = "mean"
  )
  L <- fit$metadata$data_list$L
  expect_equal(payload$data$kind, rep("density", L))
  expect_equal(payload$data$grid, fit$metadata$data_list$grid)
  expect_true(all(payload$data$lower <= payload$data$point))
  expect_true(all(payload$data$point <= payload$data$upper))
  expect_equal(sum(payload$data$point), 1, tolerance = 1e-8)
})

test_that("for a fit with neither g nor log_g, the base plot of g has the intervals of the moments of g", {
  fit <- plot_test_fit()
  fit$draws <- fit$draws[, , c("mean_g", "var_g"), drop = FALSE]

  payload <- plot_test_ns(".bef_plot_payload_bef_fit_re")(
    fit,
    type = "g",
    level = 0.8,
    sort_by = "mean"
  )

  expect_named(payload$data, c("site", "lower", "upper", "point", "kind", "position"))
  expect_equal(payload$data$kind, rep("moment", 3L))
  expect_equal(payload$data$site, c("mean_g", "var_g", "sd_g"))
})

test_that("plot() of a fit returns ggplot objects without the attribute of version 0.3", {
  fit <- plot_test_fit()

  gg <- plot(fit, type = "caterpillar", level = 0.8)
  expect_s3_class(gg, "ggplot")
  expect_null(attr(gg, "bef_plot_payload", exact = TRUE))

  g_plot <- plot(fit, type = "g", level = 0.8)
  expect_s3_class(g_plot, "ggplot")
  expect_null(attr(g_plot, "bef_plot_payload", exact = TRUE))
})
