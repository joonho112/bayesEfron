test_that("the base backend draws every plot type and returns NULL", {
  fit <- plot_test_fit()

  plot_test_with_pdf_device({
    expect_null(plot(fit, type = "caterpillar", level = 0.8, sort_by = "mean", backend = "base"))
    expect_null(plot(fit, type = "g", level = 0.8, backend = "base"))
    expect_null(plot(fit, type = "diagnostic", level = 0.8, backend = "base"))
  })
})

test_that("the base backend is deprecated and draws the four plots of version 0.3 only", {
  fit <- plot_test_fit()
  shown <- plot_test_ns(".bayesEfron_msgs_emitted")
  shown[["plot_backend_base"]] <- NULL
  withr::defer(shown[["plot_backend_base"]] <- TRUE)
  plot_test_with_pdf_device({
    expect_warning(plot(fit, type = "g", backend = "base"), "deprecated",
                   class = "bef_deprecated")
    # The warning is given once in a session.
    expect_no_warning(plot(fit, type = "g", backend = "base"))
    expect_null(plot(fit, type = "sites", backend = "base"))
    expect_null(plot(fit, type = "diagnostics", backend = "base"))
  })
  expect_error(plot(fit, type = "shrinkage", backend = "base"), "draws only the types",
               class = "bef_invalid_args")
  expect_error(plot(fit, type = "g", scale = 100, backend = "base"),
               "only further arguments are `level` and `sort_by`", class = "bef_invalid_args")

  err <- tryCatch(plot(plot_test_fit(), backend = "lattice"), error = identity)
  expect_s3_class(err, "bef_invalid_args")
  expect_identical(err$arg, "backend")
  # `backend = "ggplot2"`, which version 0.3 took, is what is drawn anyway.
  expect_s3_class(plot(bayesEfron::raudenbush_fit, backend = "ggplot2"), "ggplot")
})

test_that("the sensitivity type is no longer accepted", {
  err <- tryCatch(plot(plot_test_fit(), type = "sensitivity"), error = identity)
  expect_s3_class(err, "bef_invalid_args")
  expect_identical(err$arg, "type")
})
