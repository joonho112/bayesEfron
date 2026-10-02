test_that("the base backend draws every plot type and returns NULL", {
  fit <- plot_test_fit()

  plot_test_with_pdf_device({
    expect_null(plot(fit, type = "caterpillar", level = 0.8, sort_by = "mean", backend = "base"))
    expect_null(plot(fit, type = "g", level = 0.8, backend = "base"))
    expect_null(plot(fit, type = "diagnostic", level = 0.8, backend = "base"))
  })
})

test_that("backend picks ggplot2 when it is installed and can be set to base", {
  use_ggplot2 <- plot_test_ns(".bef_plot_use_ggplot2")
  expect_equal(use_ggplot2(), requireNamespace("ggplot2", quietly = TRUE))
  expect_false(use_ggplot2("base", explicit = TRUE))

  err <- tryCatch(plot(plot_test_fit(), backend = "lattice"), error = identity)
  expect_s3_class(err, "bef_invalid_args")
  expect_identical(err$arg, "backend")
})

test_that("the sensitivity type is no longer accepted", {
  err <- tryCatch(plot(plot_test_fit(), type = "sensitivity"), error = identity)
  expect_s3_class(err, "bef_invalid_args")
  expect_identical(err$arg, "type")
})
