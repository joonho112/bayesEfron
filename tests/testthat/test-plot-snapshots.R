# Pictures of the plots, compared with stored ones. A new version of
# ggplot2 can change a picture without changing what it shows, so the
# comparison is made where BAYESEFRON_SNAPSHOTS is "true": in one job of
# the checks and in the test runs before a release, and never on CRAN.
# After a change that is meant, review the new pictures with
# testthat::snapshot_review().

snapshot_plot <- function(title, plot) {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not(identical(Sys.getenv("BAYESEFRON_SNAPSHOTS"), "true"),
              "BAYESEFRON_SNAPSHOTS is not \"true\"")
  vdiffr::expect_doppelganger(title, plot)
}

test_that("plot_g() looks as stored", {
  fit <- plot_fit(12L, labels = TRUE)
  snapshot_plot("g density with rug and normal", plot_g(fit, observed = "rug", normal = TRUE))
  snapshot_plot("g distribution function", plot_g(fit, type = "cdf", inner_level = 0.5))
})

test_that("plot_sites() looks as stored", {
  snapshot_plot("sites with labels", plot_sites(plot_fit(12L, labels = TRUE), highlight = 3))
  snapshot_plot("sites without labels", suppressMessages(plot_sites(plot_fit(150L))))
})

test_that("plot_shrinkage() looks as stored", {
  fit <- plot_fit(12L, labels = TRUE)
  snapshot_plot("shrinkage scatter", plot_shrinkage(fit))
  snapshot_plot("shrinkage arrows", plot_shrinkage(fit, style = "arrows"))
})

test_that("plot_exceedance() looks as stored", {
  snapshot_plot("exceedance with cutoff",
                plot_exceedance(plot_fit(12L, labels = TRUE), cutoff = 0.9))
})

test_that("plot_diagnostics() looks as stored", {
  snapshot_plot("diagnostics", plot_diagnostics(plot_fit(12L)))
})

test_that("plot_loo() looks as stored", {
  skip_if_not_installed("loo")
  fit <- plot_fit(12L, labels = TRUE)
  result <- suppressWarnings(loo::loo(fit))
  result$diagnostics$pareto_k <- seq(-0.1, 0.9, length.out = 12L)
  snapshot_plot("pareto k", plot_loo(fit, loo = result))
})

test_that("plot() of a grid and of a comparison look as stored", {
  data <- plot_fit(12L)$metadata$data_list
  grid <- make_efron_grid(unname(data$theta_hat), data$sigma, grid_method = "paper_realdata")
  snapshot_plot("grid", plot(grid))
  skip_if_not_installed("loo")
  snapshot_plot("comparison", plot(plot_comparison()))
})

test_that("plot_recovery() looks as stored", {
  fit <- plot_fit(12L, labels = TRUE)
  truth <- withr::with_seed(5L, {
    table <- site_summary(fit)
    table$mean + stats::rnorm(12L, 0, 1.2 * table$sd)
  })
  snapshot_plot("recovery intervals", plot_recovery(fit, truth))
  many <- plot_fit(150L)
  truth <- withr::with_seed(5L, {
    table <- site_summary(many)
    table$mean + stats::rnorm(150L, 0, table$sd)
  })
  snapshot_plot("recovery calibration", plot_recovery(many, truth, type = "calibration"))
})
