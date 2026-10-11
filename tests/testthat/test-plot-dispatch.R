# plot() of a fit hands the fit to the function of the type, and the calls
# of version 0.3 still run.

test_that("each type is drawn by its function", {
  fit <- plot_fit(12L, labels = TRUE)
  same <- function(a, b, layer = 1L) {
    expect_equal(plot_layer_data(a, layer), plot_layer_data(b, layer))
    expect_identical(plot_caption(a), plot_caption(b))
  }
  same(plot(fit), plot_sites(fit), 2L)
  same(plot(fit, type = "sites", sort_by = "sigma", level = 0.8),
       plot_sites(fit, sort_by = "sigma", level = 0.8), 2L)
  same(plot(fit, type = "g", scale = 10, normal = TRUE), plot_g(fit, scale = 10, normal = TRUE), 2L)
  same(plot(fit, type = "shrinkage", style = "arrows"), plot_shrinkage(fit, style = "arrows"), 2L)
  same(plot(fit, type = "exceedance", cutoff = 0.9), plot_exceedance(fit, cutoff = 0.9), 3L)
  same(plot(fit, type = "diagnostics", what = "rhat"), plot_diagnostics(fit, what = "rhat"), 3L)
  same(plot(fit, type = "site_posterior", sites = 2:3), plot_site_posterior(fit, sites = 2:3), 2L)
  same(plot(fit, type = "ppc"), plot_ppc(fit), 3L)
  truth <- seq(-1, 1, length.out = 12L)
  same(plot(fit, type = "recovery", truth = truth), plot_recovery(fit, truth), 1L)
  skip_if_not_installed("loo")
  same(suppressWarnings(plot(fit, type = "loo")), suppressWarnings(plot_loo(fit)), 2L)
  skip_if_not_installed("bayesplot")
  expect_s3_class(plot(fit, type = "mcmc", pars = "lambda"), "ggplot")
})

test_that("the grid of a fit is drawn with the data of the fit", {
  fit <- plot_fit(12L)
  plot <- plot(fit, type = "grid", scale = 10)
  expect_identical(unname(plot_geoms(plot)),
                   c("GeomVline", "GeomLine", "GeomRug", "GeomSegment", "GeomPoint"))
  expect_equal(plot_layer_data(plot, 5L)$x, 10 * sort(unname(fit$metadata$data_list$theta_hat)))
  expect_equal(plot_layer_data(plot, 3L)$x, 10 * fit$metadata$data_list$grid)

  # A fit of version 0.3 holds no grid object.
  old <- readRDS(test_path("_fixtures", "short_fit.rds"))
  attr(old$metadata, "grid") <- NULL
  expect_error(suppressWarnings(plot(old, type = "grid")), "make_efron_grid",
               class = "bef_invalid_fit")
})

test_that("the names and the arguments of version 0.3 are still taken", {
  fit <- plot_fit(12L, labels = TRUE)
  expect_equal(plot_layer_data(plot(fit, type = "caterpillar"), 2L),
               plot_layer_data(plot_sites(fit), 2L))
  expect_equal(plot_layer_data(plot(fit, type = "diagnostic"), 3L),
               plot_layer_data(plot_diagnostics(fit), 3L))
  # `level` and `sort_by` were taken for every type and ignored where they
  # had no meaning.
  expect_s3_class(plot(fit, type = "diagnostic", level = 0.8, sort_by = "mean"), "ggplot")
  expect_s3_class(plot(fit, type = "g", level = 0.8, sort_by = "sigma"), "ggplot")
  expect_equal(plot_layer_data(plot(fit, type = "g", level = 0.8, sort_by = "sigma"), 1L),
               plot_layer_data(plot_g(fit, level = 0.8), 1L))
  # Under its new name a type takes the arguments of its function only.
  expect_error(plot(fit, type = "diagnostics", level = 0.8), "no argument",
               class = "bef_invalid_args")
  expect_error(plot(fit, type = "shrinkage", level = 0.8), "no argument",
               class = "bef_invalid_args")

  # type = "g" now draws the density; the probabilities are another type of
  # plot_g().
  expect_identical(plot(fit, type = "g")$labels$y, "Density")
  expect_error(plot(fit, type = "sensitivity"), "must be one of", class = "bef_invalid_args")
  expect_error(plot(fit, type = c("g", "sites")), class = "bef_invalid_args")
})

test_that("the eight plot() calls of the vignettes of version 0.3 run", {
  fit <- bayesEfron::krw_fit
  studies <- bayesEfron::raudenbush_fit
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  withr::defer({
    grDevices::dev.off()
    unlink(path)
  })
  # vignettes/bayesEfron.Rmd
  expect_null(plot(fit, type = "g", backend = "base"))
  expect_null(plot(fit, type = "caterpillar", backend = "base"))
  # vignettes/diagnostics.Rmd
  expect_s3_class(plot(fit, type = "caterpillar", sort_by = "sigma", level = 0.95), "ggplot")
  relabeled <- plot(fit, type = "g") +
    ggplot2::labs(title = NULL, x = "Callback probability difference",
                  y = "Probability of the grid point")
  expect_s3_class(relabeled, "ggplot")
  expect_s3_class(plot(fit, type = "diagnostic"), "ggplot")
  # vignettes/teacher-expectancy.Rmd
  expect_null(plot(studies, type = "g", backend = "base"))
  expect_null(plot(studies, type = "caterpillar", backend = "base"))
  skip_if_not_installed("loo")
  expect_null(plot(fit, type = "loo", backend = "base"))
})
