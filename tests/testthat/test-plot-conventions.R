# What the plot functions have in common: the order of their arguments, a
# ggplot object without a title and with a caption, and no silent
# arguments. R/plot-utils.R states the conventions.

plot_function_names <- function() {
  sort(grep("^plot_", getNamespaceExports("bayesEfron"), value = TRUE))
}

# One call of every plot function and method on a small fit.
plot_convention_plots <- function() {
  fit <- plot_fit(12L, labels = TRUE)
  many <- plot_fit(150L)
  truth <- seq(-1, 1, length.out = 12L)
  data <- fit$metadata$data_list
  plots <- list(
    plot_g = plot_g(fit, observed = "rug", normal = TRUE),
    plot_g_cdf = plot_g(fit, type = "cdf", inner_level = 0.5),
    plot_sites = plot_sites(fit),
    plot_sites_many = suppressMessages(plot_sites(many)),
    plot_shrinkage = plot_shrinkage(fit),
    plot_shrinkage_arrows = plot_shrinkage(fit, style = "arrows", normal = TRUE),
    plot_exceedance = plot_exceedance(fit, cutoff = 0.9),
    plot_exceedance_many = suppressMessages(plot_exceedance(many)),
    plot_diagnostics = plot_diagnostics(fit),
    plot_recovery = plot_recovery(fit, truth),
    plot_recovery_calibration = plot_recovery(fit, truth, type = "calibration"),
    plot_recovery_scatter = plot_recovery(fit, truth, type = "scatter"),
    plot_ppc = plot_ppc(fit),
    plot_ppc_hist = plot_ppc(many, type = "pit_hist"),
    plot_ppc_marginal = plot_ppc(fit, type = "marginal"),
    plot_site_posterior = plot_site_posterior(fit, sites = 1:2),
    plot_grid = plot(make_efron_grid(unname(data$theta_hat), data$sigma,
                                     grid_method = "paper_realdata")),
    plot_fit_grid = plot(fit, type = "grid")
  )
  if (requireNamespace("loo", quietly = TRUE)) {
    plots$plot_loo <- suppressWarnings(plot_loo(fit))
    plots$plot_comparison <- plot(plot_comparison())
  }
  plots
}

test_that("the arguments of every plot function come in the same order", {
  expect_identical(
    plot_function_names(),
    c("plot_diagnostics", "plot_exceedance", "plot_g", "plot_loo", "plot_mcmc", "plot_ppc",
      "plot_recovery", "plot_shrinkage", "plot_site_posterior", "plot_sites")
  )
  for (name in plot_function_names()) {
    args <- names(formals(get(name, asNamespace("bayesEfron"))))
    expect_identical(args[[1L]], "x", info = name)
    expect_identical(args[[length(args)]], "...", info = name)
    if ("scale" %in% args) {
      expect_identical(which(args == "scale"), length(args) - 1L, info = name)
    }
    if ("max_sites" %in% args) {
      expect_identical(which(args == "max_sites"), which(args == "scale") - 1L, info = name)
    }
    if (all(c("level", "inner_level") %in% args)) {
      expect_identical(which(args == "inner_level"), which(args == "level") + 1L, info = name)
    }
    # What is drawn comes before the levels.
    if (all(c("type", "level") %in% args)) {
      expect_lt(which(args == "type"), which(args == "level"))
    }
  }
  for (name in c("plot.bef_grid", "plot.bef_grid_comparison")) {
    args <- names(formals(get(name, asNamespace("bayesEfron"))))
    expect_identical(utils::tail(args, 2L), c("scale", "..."), info = name)
  }
})

test_that("the defaults that the functions share are the same", {
  default <- function(name, arg) eval(formals(get(name, asNamespace("bayesEfron")))[[arg]])
  expect_identical(default("plot_g", "level"), 0.9)
  expect_identical(default("plot_sites", "level"), 0.9)
  expect_identical(default("plot_recovery", "level"), 0.9)
  expect_identical(default("plot_sites", "max_sites"), 100)
  expect_identical(default("plot_exceedance", "max_sites"), 100)
  for (name in plot_function_names()) {
    if ("scale" %in% names(formals(get(name, asNamespace("bayesEfron"))))) {
      expect_identical(default(name, "scale"), 1, info = name)
    }
  }
})

test_that("every plot is a ggplot object without a title and with a caption", {
  plots <- plot_convention_plots()
  for (name in names(plots)) {
    plot <- plots[[name]]
    expect_s3_class(plot, "ggplot")
    expect_null(plot$labels$title, info = name)
    expect_null(plot$labels$subtitle, info = name)
    caption <- plot$labels$caption
    expect_true(is.character(caption) && nzchar(caption), info = name)
    # Whole sentences, wrapped to the width of a plot.
    expect_match(caption, "[.]$", info = name)
    expect_lte(max(nchar(strsplit(caption, "\n", fixed = TRUE)[[1L]])), 100L)
    expect_no_match(caption, "  ", fixed = TRUE, info = name)
    # It can be built, which finds a wrong aesthetic or scale.
    expect_no_error(ggplot2::ggplot_build(plot))
    expect_null(attr(plot, "bef_plot_payload", exact = TRUE), info = name)
  }
})

test_that("an argument that a function does not have is an error", {
  fit <- plot_fit(12L)
  truth <- seq(-1, 1, length.out = 12L)
  required <- list(plot_recovery = list(truth = truth), plot_site_posterior = list(sites = 1L))
  for (name in plot_function_names()) {
    fun <- get(name, asNamespace("bayesEfron"))
    args <- c(list(fit), required[[name]])
    expect_error(do.call(fun, c(args, list(not_an_argument = 1))), "no argument",
                 class = "bef_invalid_args", info = name)
    expect_error(do.call(fun, c(list(list()), required[[name]])), class = "bef_invalid_fit",
                 info = name)
    if ("scale" %in% names(formals(fun))) {
      expect_error(do.call(fun, c(args, list(scale = 0))), class = "bef_invalid_args",
                   info = name)
      expect_error(do.call(fun, c(args, list(scale = c(1, 100)))), class = "bef_invalid_args",
                   info = name)
    }
  }
})

test_that("`scale` multiplies the values of the effects and leaves probabilities alone", {
  fit <- plot_fit(12L, labels = TRUE)
  sites <- plot_layer_data(plot_sites(fit, scale = 100), 2L)
  expect_equal(sites$x, 100 * plot_layer_data(plot_sites(fit), 2L)$x)
  density <- plot_layer_data(plot_g(fit, scale = 100), 2L)
  expect_equal(sum(density$y) * diff(density$x[1:2]), 1, tolerance = 1e-8)
  expect_equal(plot_layer_data(plot_exceedance(fit, scale = 100), 2L)$x,
               plot_layer_data(plot_exceedance(fit), 2L)$x)
})
