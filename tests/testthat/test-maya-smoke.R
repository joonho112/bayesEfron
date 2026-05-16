maya_example_path <- function() {
  installed <- system.file(
    "examples",
    "example-bayes-efron-fit.R",
    package = "bayesEfron"
  )
  if (nzchar(installed) && file.exists(installed)) {
    return(installed)
  }
  root <- bef_test_source_root(
    required = file.path("inst", "examples", "example-bayes-efron-fit.R")
  )
  file.path(root, "inst", "examples", "example-bayes-efron-fit.R")
}

test_that("Maya smoke example runs end-to-end under live gate", {
  skip_if_not(
    live_cmdstan_run_requested(),
    "Maya smoke path is gated by BAYESEFRON_RUN_LIVE=1 or BAYESEFRON_RUN_FULL_LIVE=1."
  )
  live_cmdstan_skip_if_unavailable()
  skip_if_not_installed("metafor")

  with_live_cmdstan_cache_root({
    withr::local_envvar(c(BAYESEFRON_NO_GGPLOT2 = "1"))
    grDevices::pdf(tempfile(fileext = ".pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)

    env <- new.env(parent = globalenv())
    sys.source(maya_example_path(), envir = env)
    result <- env$maya_result

    expect_type(result, "list")
    expect_s3_class(result$data, "bef_data")
    expect_identical(result$data$source, "metafor::escalc")
    expect_equal(length(result$data$theta_hat), 10L)
    expect_true(all(is.finite(result$data$theta_hat)))
    expect_true(all(is.finite(result$data$sigma)))
    expect_true(all(result$data$sigma > 0))
    expect_s3_class(result$fit, "bef_fit_re")
    expect_s3_class(result$fit, "bef_fit")
    expect_false("cmdstan_fit" %in% names(result$fit))
    expect_equal(result$fit$metadata$model_family, "RE")
    expect_equal(result$fit$metadata$data_list$K, 10L)
    expect_s3_class(result$summary, "summary.bef_fit_re")
    expect_true("theta_summary" %in% names(result$summary))
    expect_s3_class(result$diagnostic, "bef_diagnostic")
    expect_named(result$theta_ci, c("site", "lower", "upper", "point"))
    expect_equal(nrow(result$theta_ci), 10L)
    expect_true(all(is.finite(result$theta_ci$lower)))
    expect_true(all(is.finite(result$theta_ci$upper)))
    expect_true(all(is.finite(result$theta_ci$point)))
    expect_true(all(result$theta_ci$lower <= result$theta_ci$upper))
    expect_null(result$plot)
  })
})
