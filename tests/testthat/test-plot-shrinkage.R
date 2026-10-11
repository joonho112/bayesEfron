# plot_shrinkage(): the estimates and the posterior means are those of
# site_summary(), and the normal model beside them is the linear one.

test_that("the scatter sets the posterior means against the estimates", {
  fit <- plot_fit(12L, labels = TRUE)
  table <- site_summary(fit)
  plot <- plot_shrinkage(fit)
  expect_s3_class(plot, "ggplot")
  expect_identical(unname(plot_geoms(plot)), c("GeomAbline", "GeomPoint"))
  expect_null(plot$labels$title)
  expect_identical(plot$labels$x, "Estimate")
  expect_identical(plot$labels$y, "Posterior mean")

  points <- plot_layer_data(plot, 2L)
  expect_equal(points$x, table$estimate)
  expect_equal(points$y, table$mean)
  expect_equal(points$y, unname(coef(fit)))
  expect_equal(plot$data$color, table$se)
  expect_identical(plot$scales$get_scales("colour")$name, "Standard error")
  # Both axes have the same limits, those of the estimates and the means.
  expect_equal(plot$coordinates$limits$x, range(table$estimate, table$mean))
  expect_identical(plot$coordinates$limits$y, plot$coordinates$limits$x)
  expect_match(plot_caption(plot), "Dashed line: the two are equal")
  expect_match(plot_caption(plot), "Color: standard error of the estimate")

  scaled <- plot_layer_data(plot_shrinkage(fit, scale = 100), 2L)
  expect_equal(scaled$x, 100 * points$x)
  expect_equal(scaled$y, 100 * points$y)

  plain <- plot_shrinkage(fit, color_by = "none")
  expect_null(plain$scales$get_scales("colour"))
  expect_no_match(plot_caption(plain), "Color")
})

test_that("the arrows run from the estimates to the posterior means", {
  fit <- plot_fit(12L, labels = TRUE)
  table <- site_summary(fit)
  plot <- plot_shrinkage(fit, style = "arrows")
  expect_identical(unname(plot_geoms(plot)),
                   c("GeomVline", "GeomSegment", "GeomPoint", "GeomText"))
  expect_identical(plot$labels$y, "Standard error")
  arrows <- plot_layer_data(plot, 2L)
  expect_equal(arrows$x, table$estimate)
  expect_equal(arrows$xend, table$mean)
  expect_equal(arrows$y, table$se)
  expect_equal(arrows$yend, table$se)
  expect_s3_class(plot$layers[[2L]]$geom_params$arrow, "arrow")
  expect_equal(plot_layer_data(plot, 1L)$xintercept, mean(fit$posterior$mean_g))
  expect_identical(plot_layer_data(plot, 4L)$label, table$label)
  # The standard error is the vertical axis and not the color.
  expect_null(plot$scales$get_scales("colour"))
  expect_match(plot_caption(plot), "Each arrow runs from the estimate of a site \\(open circle\\)")

  # Many sites: lines without heads, a point at the posterior mean, no
  # labels.
  many <- plot_fit(150L)
  lines <- plot_shrinkage(many, style = "arrows")
  expect_identical(unname(plot_geoms(lines)), c("GeomVline", "GeomSegment", "GeomPoint"))
  expect_null(lines$layers[[2L]]$geom_params$arrow)
  expect_equal(plot_layer_data(lines, 3L)$x, site_summary(many)$mean)
  expect_match(plot_caption(lines), "which a point marks")
})

test_that("many points are transparent", {
  alpha <- getFromNamespace(".bef_plot_alpha", "bayesEfron")
  expect_identical(alpha(12L), 1)
  expect_identical(alpha(200L), 0.5)
  expect_identical(alpha(1500L), 0.15)
  expect_equal(unique(plot_layer_data(plot_shrinkage(plot_fit(150L)), 2L)$alpha), 100 / 150)
})

test_that("the normal model is drawn in a second panel", {
  fit <- plot_fit(12L, labels = TRUE)
  table <- site_summary(fit)
  reml <- getFromNamespace(".bef_reml_normal", "bayesEfron")(table$estimate, table$se)
  plot <- plot_shrinkage(fit, normal = TRUE)
  expect_identical(levels(plot$data$model), c("This fit", "Normal (REML)"))
  points <- plot_layer_data(plot, 2L)
  expect_identical(nrow(points), 24L)
  second <- points[points$PANEL == 2L, ]
  share <- reml$sd^2 / (reml$sd^2 + table$se^2)
  expect_equal(second$x, table$estimate)
  expect_equal(second$y, reml$mean + share * (table$estimate - reml$mean))
  expect_equal(points$y[points$PANEL == 1L], table$mean)
  expect_match(plot_caption(plot), "estimated by REML and held there")

  given <- plot_shrinkage(fit, normal = c(mean = 0.5, sd = 0), scale = 10)
  expect_identical(levels(given$data$model), c("This fit", "Normal"))
  points <- plot_layer_data(given, 2L)
  # Without variation of the effects every site is moved to the mean.
  expect_equal(points$y[points$PANEL == 2L], rep(5, 12L))
  expect_match(plot_caption(given), "mean and standard deviation given")

  arrows <- plot_shrinkage(fit, style = "arrows", normal = TRUE)
  centers <- plot_layer_data(arrows, 1L)
  expect_equal(centers$xintercept, c(mean(fit$posterior$mean_g), reml$mean))
  expect_error(plot_shrinkage(fit, normal = c(1, 2)), class = "bef_invalid_args")
})

test_that("true effects are the color, and crosses in the arrows style", {
  fit <- plot_fit(12L, labels = TRUE)
  truth <- seq(-1, 1, length.out = 12L)
  plot <- plot_shrinkage(fit, truth = truth)
  expect_equal(plot$data$color, truth)
  expect_identical(plot$scales$get_scales("colour")$name, "True effect")
  expect_match(plot_caption(plot), "Color: true effect")

  # A choice of the color is kept.
  by_sigma <- plot_shrinkage(fit, truth = truth, color_by = "sigma")
  expect_equal(by_sigma$data$color, site_summary(fit)$se)
  expect_error(plot_shrinkage(fit, color_by = "truth"), "needs the true effects",
               class = "bef_invalid_args")
  expect_error(plot_shrinkage(fit, truth = truth[-1L]), class = "bef_invalid_args")

  arrows <- plot_shrinkage(fit, style = "arrows", truth = truth, scale = 10)
  expect_identical(unname(plot_geoms(arrows)),
                   c("GeomVline", "GeomSegment", "GeomPoint", "GeomPoint", "GeomText"))
  expect_equal(plot_layer_data(arrows, 4L)$x, 10 * truth)
  expect_match(plot_caption(arrows), "Crosses: the true effects")
})

test_that("highlighted sites are marked and labelled", {
  fit <- plot_fit(150L, labels = TRUE)
  labels <- site_summary(fit)$label
  colors <- getFromNamespace(".bef_plot_colors", "bayesEfron")()

  plot <- suppressMessages(plot_shrinkage(fit, highlight = labels[c(4L, 90L)]))
  expect_identical(unname(plot_geoms(plot)),
                   c("GeomAbline", "GeomPoint", "GeomPoint", "GeomText"))
  marked <- plot_layer_data(plot, 3L)
  expect_equal(marked$x, site_summary(fit)$estimate[c(4L, 90L)])
  expect_identical(unique(marked$colour), colors[["highlight"]])
  expect_identical(plot_layer_data(plot, 4L)$label, labels[c(4L, 90L)])
  expect_match(plot_caption(plot), "Orange: the highlighted sites")

  arrows <- plot_shrinkage(fit, style = "arrows", highlight = c(4, 90))
  expect_identical(unname(plot_geoms(arrows)),
                   c("GeomVline", "GeomSegment", "GeomPoint", "GeomSegment", "GeomText"))
  expect_identical(plot_layer_data(arrows, 5L)$label, labels[c(4L, 90L)])
  expect_error(plot_shrinkage(fit, highlight = "no such site"), class = "bef_invalid_args")
})

test_that("wrong arguments are errors", {
  fit <- plot_fit(12L)
  expect_error(plot_shrinkage(fit, style = "lines"), class = "bef_invalid_args")
  expect_error(plot_shrinkage(fit, color_by = "se"), class = "bef_invalid_args")
  expect_error(plot_shrinkage(fit, scale = 0), class = "bef_invalid_args")
  expect_error(plot_shrinkage(fit, alpha = 0.5), "no argument", class = "bef_invalid_args")
  expect_error(plot_shrinkage(1:3), class = "bef_invalid_fit")
})
