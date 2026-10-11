# plot() of a grid: the basis functions, the knots and the grid points are
# those of the grid, and the estimates are those that it was built from.

grid_data <- function(K = 12L) {
  withr::with_seed(11L, {
    sigma <- stats::runif(K, 0.2, 0.6)
    list(theta_hat = stats::rnorm(K, 0.5, 1) + stats::rnorm(K, 0, sigma), sigma = sigma)
  })
}

test_that("the basis, the knots and the grid points are those of the grid", {
  data <- grid_data()
  grid <- make_efron_grid(data$theta_hat, data$sigma, grid_method = "paper_realdata")
  plot <- plot(grid)
  expect_s3_class(plot, "ggplot")
  expect_identical(unname(plot_geoms(plot)),
                   c("GeomVline", "GeomLine", "GeomRug", "GeomSegment", "GeomPoint"))
  expect_null(plot$labels$title)
  expect_identical(plot$labels$x, "Effect")

  knots <- plot_layer_data(plot, 1L)
  # In both panels: the two ends and the interior knots.
  expect_equal(sort(unique(knots$xintercept)),
               sort(c(grid$knots$boundary, grid$knots$interior)))
  expect_identical(nrow(knots), 2L * (2L + length(grid$knots$interior)))
  expect_equal(grid$knots$boundary, grid$support)

  curves <- plot_layer_data(plot, 2L)
  expect_identical(nrow(curves), grid$L * grid$M)
  expect_identical(length(unique(curves$group)), grid$M)
  first <- curves[curves$group == 1L, ]
  expect_equal(first$x, grid$grid)
  expect_equal(first$y, as.numeric(grid$B[, 1L]))
  expect_equal(plot_layer_data(plot, 3L)$x, grid$grid)

  # The estimates in increasing order, with two standard errors on each side.
  by_estimate <- order(data$theta_hat)
  lines <- plot_layer_data(plot, 4L)
  expect_equal(lines$x, (data$theta_hat - 2 * data$sigma)[by_estimate])
  expect_equal(lines$xend, (data$theta_hat + 2 * data$sigma)[by_estimate])
  expect_equal(plot_layer_data(plot, 5L)$x, sort(data$theta_hat))
  expect_match(plot_caption(plot), "Curves: the 6 basis functions of the natural cubic spline")
  expect_match(plot_caption(plot), "Ticks: the 101 grid points")
  expect_match(plot_caption(plot), "the 12 estimates, in increasing order")
  expect_no_match(plot_caption(plot), "outside the grid")

  scaled <- plot(grid, scale = 100)
  expect_equal(plot_layer_data(scaled, 3L)$x, 100 * grid$grid)
  expect_equal(plot_layer_data(scaled, 5L)$x, 100 * sort(data$theta_hat))
  expect_equal(sort(unique(plot_layer_data(scaled, 1L)$xintercept)),
               100 * sort(c(grid$knots$boundary, grid$knots$interior)))
})

test_that("estimates outside the grid are counted", {
  data <- grid_data()
  narrow <- suppressWarnings(
    make_efron_grid(data$theta_hat, data$sigma, support = stats::quantile(data$theta_hat, c(0.3, 0.8)),
                    L = 51L)
  )
  outside <- sum(data$theta_hat < narrow$support[1L] | data$theta_hat > narrow$support[2L])
  expect_gt(outside, 1L)
  expect_match(plot_caption(plot(narrow)), sprintf("%d estimates lie outside the grid", outside))
})

test_that("many estimates are a histogram, and a grid without estimates is drawn alone", {
  data <- grid_data(150L)
  grid <- make_efron_grid(data$theta_hat, data$sigma, grid_method = "paper_realdata")
  plot <- plot(grid)
  expect_identical(unname(plot_geoms(plot)), c("GeomVline", "GeomLine", "GeomRug", "GeomRect"))
  bars <- plot_layer_data(plot, 4L)
  expect_identical(sum(bars$ymax), 150)
  expect_match(plot_caption(plot), "Bars: the 150 estimates")

  # The grid that a fit stores holds no estimates.
  fit <- plot_fit(12L)
  stored <- attr(fit$metadata, "grid", exact = TRUE)
  expect_null(stored$data)
  alone <- plot(stored)
  expect_identical(unname(plot_geoms(alone)), c("GeomVline", "GeomLine", "GeomRug"))
  expect_match(plot_caption(alone), "The grid holds no estimates")
  # They can be given, as two vectors or as a data frame.
  theta_hat <- unname(fit$metadata$data_list$theta_hat)
  sigma <- fit$metadata$data_list$sigma
  given <- plot(stored, theta_hat, sigma)
  expect_equal(plot_layer_data(given, 5L)$x, sort(theta_hat))
  from_frame <- plot(stored, data.frame(yi = theta_hat, sei = sigma))
  expect_equal(plot_layer_data(from_frame, 5L)$x, sort(theta_hat))
})

test_that("wrong arguments are errors", {
  data <- grid_data()
  grid <- make_efron_grid(data$theta_hat, data$sigma, grid_method = "paper_realdata")
  expect_error(plot(grid, sigma = data$sigma), "needs the estimates", class = "bef_invalid_args")
  expect_error(plot(grid, data$theta_hat, data$sigma[-1L]), class = "bef_invalid_args")
  expect_error(plot(grid, scale = -1), class = "bef_invalid_args")
  expect_error(plot(grid, main = "Grid"), "no argument", class = "bef_invalid_args")
})

# ---- plot() of a comparison --------------------------------------------------

test_that("a comparison is drawn row by row from what it holds", {
  skip_if_not_installed("loo")
  comparison <- plot_comparison()
  plot <- plot(comparison)
  expect_s3_class(plot, "ggplot")
  expect_identical(unname(plot_geoms(plot)),
                   c("GeomBlank", "GeomRibbon", "GeomLine", "GeomSegment", "GeomPoint"))
  expect_null(plot$labels$title)
  # The first specification is the top row, named by what differs.
  expect_identical(rev(plot_axis_labels(plot)),
                   c("expansion 0.25", "expansion 0.50", "expansion 1.00"))

  g <- attr(comparison, "g", exact = TRUE)
  top <- max(vapply(g, function(curve) max(curve$upper), numeric(1L)))
  curves <- plot_layer_data(plot, 3L)
  first <- curves[curves$group == 1L, ]
  expect_equal(first$x, g[["1"]]$value)
  expect_equal(first$y, 3 + 0.85 * g[["1"]]$mean / top)
  last <- curves[curves$group == 3L, ]
  expect_equal(last$y, 1 + 0.85 * g[["3"]]$mean / top)
  bands <- plot_layer_data(plot, 2L)
  expect_equal(bands$ymax[bands$group == 2L], 2 + 0.85 * g[["2"]]$upper / top)
  # No band reaches into the row above.
  expect_true(all(bands$ymax - floor(bands$ymin + 1e-9) <= 0.85 + 1e-9))

  by_spec <- comparison[order(comparison$spec), ]
  points <- plot_layer_data(plot, 5L)
  expect_equal(points$x[points$PANEL == 2L], by_spec$mean_g)
  expect_equal(points$x[points$PANEL == 3L], by_spec$sd_g)
  expect_equal(points$y[points$PANEL == 2L], c(3, 2, 1))
  kept <- attr(comparison, "sites", exact = TRUE)
  expect_equal(points$x[points$PANEL == 4L], kept$mean[kept$site == 2L])
  expect_equal(points$x[points$PANEL == 4L], by_spec[["site_Site 002"]])
  lines <- plot_layer_data(plot, 4L)
  expect_equal(lines$x[lines$PANEL == 5L], kept$lower[kept$site == 9L])
  expect_equal(lines$xend[lines$PANEL == 5L], kept$upper[kept$site == 9L])
  panels <- ggplot2::ggplot_build(plot)$layout$layout$panel
  expect_identical(as.character(panels),
                   c("Density of g", "Mean of g", "SD of g", "Site Site 002", "Site Site 009"))
  expect_match(plot_caption(plot), "pointwise 90% posterior interval, all on one vertical scale")
  expect_no_match(plot_caption(plot), "failed")

  scaled <- plot(comparison, scale = 10)
  expect_equal(plot_layer_data(scaled, 5L)$x, 10 * points$x)
  expect_equal(plot_layer_data(scaled, 3L)$x, 10 * curves$x)
  expect_equal(plot_layer_data(scaled, 3L)$y, curves$y)
})

test_that("the stored comparison of three ranges shows what the range moves", {
  comparison <- readRDS(system.file("examples", "raudenbush_ranges.rds", package = "bayesEfron"))
  expect_s3_class(comparison, "bef_grid_comparison")
  by_range <- comparison[order(comparison$expansion), ]
  expect_identical(by_range$expansion, c(0.25, 0.5, 1))
  # Study 4, the largest estimate, is moved further as the range widens,
  # and the standard deviation of g changes little.
  expect_true(all(diff(by_range$site_4) < 0))
  expect_gt(by_range$site_4[[1L]] - by_range$site_4[[3L]], 0.2)
  expect_lt(diff(range(by_range$sd_g)) / mean(by_range$sd_g), 0.1)

  plot <- plot(comparison)
  expect_identical(plot_axis_labels(plot),
                   rev(c("expansion 0.25", "expansion 0.50", "expansion 1.00")))
  points <- plot_layer_data(plot, 5L)
  expect_equal(points$x[points$PANEL == 3L], by_range$sd_g)
  expect_equal(points$x[points$PANEL == 4L], by_range$site_4)
  panels <- ggplot2::ggplot_build(plot)$layout$layout$panel
  expect_identical(as.character(panels),
                   c("Density of g", "Mean of g", "SD of g", "Site 4", "Site 8"))
})

test_that("a specification whose fit failed keeps its row", {
  skip_if_not_installed("loo")
  comparison <- plot_comparison(fail = 0.5)
  plot <- plot(comparison)
  expect_identical(rev(plot_axis_labels(plot)),
                   c("expansion 0.25", "expansion 0.50 (fit failed)", "expansion 1.00"))
  expect_identical(sort(unique(plot_layer_data(plot, 3L)$group)), 1:2)
  expect_false(any(plot_layer_data(plot, 5L)$y == 2))
  expect_match(plot_caption(plot), "A row without marks is a specification whose fit failed")
})

test_that("specifications are named by the settings that differ", {
  labels <- getFromNamespace(".bef_comparison_labels", "bayesEfron")
  specs <- data.frame(
    spec = 1:4, L = c(51L, 101L, 51L, 101L), M = c(4L, 4L, 6L, 6L), expansion = 0.5,
    grid_method = "paper_realdata", component = "point", lower = -1, upper = 2,
    stringsAsFactors = FALSE
  )
  expect_identical(labels(specs), c("L = 51, M = 4", "L = 101, M = 4", "L = 51, M = 6",
                                    "L = 101, M = 6"))
  one <- specs[1L, ]
  expect_identical(labels(one), "L = 51, M = 4")
  ranges <- transform(specs[1:2, ], L = 101L, grid_method = "support", expansion = NA_real_,
                      lower = c(-1, -2), upper = c(2, 3))
  expect_identical(labels(ranges), c("-1 to 2", "-2 to 3"))
  mixed <- transform(specs[1:2, ], L = 101L, grid_method = c("paper_realdata", "pilot"),
                     component = c("point", "normal"))
  expect_identical(labels(mixed), c("rule \"paper_realdata\", point masses",
                                    "rule \"pilot\", normal components"))
})

test_that("a comparison without the estimates of g cannot be drawn", {
  old <- readRDS(test_path("_fixtures", "grid_comparison_v030.rds"))$raudenbush
  expect_s3_class(old, "bef_grid_comparison")
  expect_error(plot(old), "Run `compare_efron_grids\\(\\)` again", class = "bef_invalid_args")
  skip_if_not_installed("loo")
  comparison <- plot_comparison()
  expect_error(plot(comparison, level = 0.5), "no argument", class = "bef_invalid_args")
  expect_error(plot(comparison, scale = "a"), class = "bef_invalid_args")
})
