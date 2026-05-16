test_that("plot.bef_fit_re base backend draws all plot types and returns NULL", {
  fit <- plot_test_fit()
  withr::local_envvar(c(BAYESEFRON_NO_GGPLOT2 = "1"))

  plot_test_with_pdf_device({
    expect_null(plot(fit, type = "caterpillar", level = 0.8, sort_by = "mean"))
    expect_null(plot(fit, type = "g", level = 0.8))
    expect_null(plot(fit, type = "sensitivity", level = 0.8, sort_by = "sigma"))
    expect_null(plot(fit, type = "diagnostic", level = 0.8))
  })
})

test_that("BAYESEFRON_NO_GGPLOT2 forces the base backend", {
  fit <- plot_test_fit()
  use_ggplot2 <- plot_test_ns(".bef_plot_use_ggplot2")
  expect_equal(use_ggplot2(), requireNamespace("ggplot2", quietly = TRUE))

  withr::local_envvar(c(BAYESEFRON_NO_GGPLOT2 = "1"))
  expect_false(use_ggplot2())

  plot_test_with_pdf_device({
    expect_null(plot(fit, type = "caterpillar"))
  })
})
