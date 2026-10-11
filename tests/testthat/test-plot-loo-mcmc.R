# plot_loo(): the Pareto k values and the threshold are those of the loo
# package. plot_mcmc(): the draws go to bayesplot.

test_that("the threshold of the Pareto k plot is that of loo", {
  skip_if_not_installed("loo")
  threshold <- getFromNamespace(".bef_pareto_k_threshold", "bayesEfron")
  expect_equal(threshold(800), 1 - 1 / log10(800))
  expect_equal(round(threshold(800), 3L), 0.656)
  expect_identical(threshold(12000), 0.7)
  if (exists("ps_khat_threshold", envir = asNamespace("loo"), inherits = FALSE)) {
    of_loo <- get("ps_khat_threshold", envir = asNamespace("loo"))
    for (draws in c(200, 800, 2200, 12000)) {
      expect_equal(threshold(draws), of_loo(draws))
    }
  }
})

test_that("the points are the Pareto k values of loo()", {
  skip_if_not_installed("loo")
  fit <- plot_fit(12L, labels = TRUE)
  result <- suppressWarnings(loo::loo(fit))
  plot <- suppressWarnings(plot_loo(fit))
  expect_s3_class(plot, "ggplot")
  expect_null(plot$labels$title)
  expect_identical(plot$labels$y, "Pareto k")
  expect_equal(plot$data$pareto_k, as.numeric(result$diagnostics$pareto_k))
  line <- plot_layer_data(plot, 1L)$yintercept
  expect_equal(line, 1 - 1 / log10(200))
  # The sites that the plot marks are those that loo reports.
  expect_identical(which(plot$data$mark == "above"), unname(loo::pareto_k_ids(result)))
  expect_match(plot_caption(plot), "not reliable with 200 draws")
})

test_that("a result of loo() can be given, and the sites above the line are labeled", {
  skip_if_not_installed("loo")
  fit <- plot_fit(12L, labels = TRUE)
  result <- suppressWarnings(loo::loo(fit))
  result$diagnostics$pareto_k[c(3L, 8L)] <- c(0.9, 1.2)
  plot <- plot_loo(fit, loo = result)
  expect_identical(unname(plot_geoms(plot)), c("GeomHline", "GeomPoint", "GeomText"))
  expect_identical(plot_layer_data(plot, 3L)$label, c("Site 003", "Site 008"))
  expect_identical(which(plot$data$mark == "above"), c(3L, 8L))
  expect_match(plot_caption(plot), "2 of the 12 sites are above it \\(orange triangles\\)")

  unlabeled <- plot_loo(fit, loo = result, label = FALSE)
  expect_identical(unname(plot_geoms(unlabeled)), c("GeomHline", "GeomPoint"))
  result$diagnostics$pareto_k[] <- 0.1
  none <- plot_loo(fit, loo = result)
  expect_identical(unname(plot_geoms(none)), c("GeomHline", "GeomPoint"))
  expect_match(plot_caption(none), "No site is above it")

  other <- suppressWarnings(loo::loo(plot_fit(150L)))
  expect_error(plot_loo(fit, loo = other), "each of the 12 sites", class = "bef_invalid_args")
  expect_error(plot_loo(fit, loo = list()), class = "bef_invalid_args")
  expect_error(plot_loo(fit, label = NA), class = "bef_invalid_args")
  expect_error(plot_loo(fit, threshold = 0.7), "no argument", class = "bef_invalid_args")
})

test_that("`pars` names parameters, with or without an index", {
  variables <- getFromNamespace(".bef_plot_mcmc_variables", "bayesEfron")
  all <- c("lp__", "alpha[1]", "alpha[2]", "lambda", "log_g[1]", "log_g[2]", "mean_g", "sd_g")
  expect_identical(variables(c("alpha", "lambda"), all), c("alpha[1]", "alpha[2]", "lambda"))
  expect_identical(variables(c("sd_g", "alpha[2]", "alpha"), all),
                   c("sd_g", "alpha[2]", "alpha[1]"))
  expect_error(variables("theta", all), "names no parameter of the fit: \"theta\"",
               class = "bef_invalid_args")
  expect_error(variables(character(), all), class = "bef_invalid_args")
  expect_error(variables(1:2, all), class = "bef_invalid_args")
})

test_that("the plots of bayesplot are drawn from a fit", {
  skip_if_not_installed("bayesplot")
  fit <- plot_fit(12L)
  trace <- plot_mcmc(fit)
  expect_s3_class(trace, "ggplot")
  # The spline coefficients, lambda, and the mean and the SD of g.
  expect_identical(levels(trace$data$parameter),
                   c(sprintf("alpha[%d]", 1:6), "lambda", "mean_g", "sd_g"))
  expect_s3_class(plot_mcmc(fit, type = "rank", pars = "lambda"), "ggplot")
  expect_s3_class(plot_mcmc(fit, type = "dens", pars = c("mean_g", "sd_g")), "ggplot")
  expect_s3_class(plot_mcmc(fit, type = "energy"), "ggplot")
  # The grid of plots is laid out on a device when it is made.
  grDevices::pdf(withr::local_tempfile(fileext = ".pdf"))
  withr::defer(grDevices::dev.off())
  expect_s3_class(plot_mcmc(fit, type = "pairs", pars = c("lambda", "sd_g")), "bayesplot_grid")

  # A thinned fit does not hold neighboring iterations, which the energy plot
  # compares; the other plots are drawn from the draws that it keeps.
  expect_error(plot_mcmc(bayesEfron::krw_fit, type = "energy"),
               "keeps one iteration in 15", class = "bef_invalid_args")
  expect_s3_class(plot_mcmc(bayesEfron::krw_fit, type = "dens", pars = "sd_g"), "ggplot")

  expect_error(plot_mcmc(fit, type = "violin"), class = "bef_invalid_args")
  expect_error(plot_mcmc(fit, pars = "theta"), class = "bef_invalid_args")
  expect_error(plot_mcmc(fit, regex_pars = "alpha"), "no argument", class = "bef_invalid_args")
  expect_error(plot_mcmc(1), class = "bef_invalid_fit")
})
