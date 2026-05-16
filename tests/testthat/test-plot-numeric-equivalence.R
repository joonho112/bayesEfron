test_that("plot payload contains backend-independent caterpillar numeric content", {
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

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- plot_test_ns(".bef_plot_ggplot2")(payload)
    expect_s3_class(gg, "ggplot")
    expect_identical(attr(gg, "bef_plot_payload", exact = TRUE), payload)
  }
})

test_that("plot payload uses grid-level g density draws when available", {
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

test_that("plot payload falls back to prior-moment g summaries without g draws", {
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

test_that("plot.bef_fit_re returns ggplot objects when ggplot2 backend is active", {
  skip_if_not_installed("ggplot2")
  fit <- plot_test_fit()
  withr::local_envvar(c(BAYESEFRON_NO_GGPLOT2 = NA))

  gg <- plot(fit, type = "caterpillar", level = 0.8)
  expect_s3_class(gg, "ggplot")
  expect_equal(attr(gg, "bef_plot_payload", exact = TRUE)$type, "caterpillar")

  g_plot <- plot(fit, type = "g", level = 0.8)
  expect_s3_class(g_plot, "ggplot")
  expect_equal(attr(g_plot, "bef_plot_payload", exact = TRUE)$type, "g")
})
