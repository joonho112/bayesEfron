# plot_g(): the plot draws the numbers of g_summary(), and what is laid over
# the estimate of g is what was asked for.

# log_g of a fit as a plain matrix.
gd_log_g <- function(fit) {
  getFromNamespace(".bef_engine", "bayesEfron")(fit)$log_g
}

test_that("the density plot draws the posterior mean and the band of g_summary()", {
  fit <- plot_fit(12L)
  plot <- plot_g(fit)
  expect_s3_class(plot, "ggplot")
  expect_identical(unname(plot_geoms(plot)), c("GeomRibbon", "GeomLine"))
  expect_null(plot$labels$title)
  expect_identical(plot$labels$x, "Effect")
  expect_identical(plot$labels$y, "Density")
  expect_match(plot_caption(plot), "pointwise 90% posterior interval")

  summary <- g_summary(fit)
  band <- plot_layer_data(plot, 1L)
  line <- plot_layer_data(plot, 2L)
  expect_equal(band$x, summary$value)
  expect_equal(band$ymin, summary$lower)
  expect_equal(band$ymax, summary$upper)
  expect_equal(line$y, summary$mean)
  # The density integrates to one over the grid.
  expect_equal(sum(line$y) * diff(summary$value[1:2]), 1, tolerance = 1e-8)

  # In other units the values are multiplied and the density is divided.
  scaled <- plot_layer_data(plot_g(fit, scale = 100), 2L)
  expect_equal(scaled$x, 100 * line$x)
  expect_equal(scaled$y, line$y / 100)

  wide <- plot_g(fit, level = 0.5)
  expect_equal(plot_layer_data(wide, 1L)$ymax, g_summary(fit, level = 0.5)$upper)
  expect_match(plot_caption(wide), "pointwise 50% posterior interval")
})

test_that("the distribution function and the probabilities are those of g_summary()", {
  fit <- plot_fit(12L)
  grid <- fit$metadata$data_list$grid
  spacing <- diff(grid[1:2])

  cdf <- plot_g(fit, type = "cdf")
  expect_identical(cdf$labels$y, "Cumulative probability")
  line <- plot_layer_data(cdf, 2L)
  # Evaluated at the edges of the cells, between which it is linear.
  edges <- c(grid[1L] - spacing / 2, grid + spacing / 2)
  expect_equal(line$x, edges)
  expect_equal(line$y, g_summary(fit, type = "cdf", at = edges)$mean)
  expect_equal(line$y[[1L]], 0)
  expect_equal(line$y[[length(line$y)]], 1)
  expect_true(all(diff(line$y) >= 0))

  mass <- plot_g(fit, type = "mass")
  expect_identical(mass$labels$y, "Probability of a grid point")
  drawn <- plot_layer_data(mass, 2L)
  # The numbers that plot(type = "g") drew in version 0.3.0: the mean and
  # the quantiles of the probabilities over the draws.
  g <- unname(exp(gd_log_g(fit)))
  expect_equal(drawn$x, grid)
  expect_equal(drawn$y, colMeans(g), tolerance = 1e-8)
  expect_equal(plot_layer_data(mass, 1L)$ymax,
               apply(g, 2L, stats::quantile, probs = 0.95, names = FALSE), tolerance = 1e-8)
  expect_equal(sum(drawn$y), 1, tolerance = 1e-8)
})

test_that("a second band, single draws and a rug are added as layers", {
  fit <- plot_fit(12L)
  two <- plot_g(fit, inner_level = 0.5)
  expect_identical(unname(plot_geoms(two)), c("GeomRibbon", "GeomRibbon", "GeomLine"))
  expect_equal(plot_layer_data(two, 2L)$ymax, g_summary(fit, level = 0.5)$upper)
  expect_match(plot_caption(two), "with the 50% interval inside")
  expect_error(plot_g(fit, inner_level = 0.9), "smaller than `level`", class = "bef_invalid_args")

  drawn <- plot_g(fit, draws = 7, observed = "rug")
  expect_identical(unname(plot_geoms(drawn)), c("GeomRibbon", "GeomLine", "GeomLine", "GeomRug"))
  single <- plot_layer_data(drawn, 2L)
  expect_length(unique(single$group), 7L)
  # Each draw is a density: it integrates to one.
  spacing <- diff(fit$metadata$data_list$grid[1:2])
  expect_equal(as.numeric(tapply(single$y, single$group, sum)) * spacing, rep(1, 7L),
               tolerance = 1e-8)
  expect_equal(sort(plot_layer_data(drawn, 4L)$x), sort(fit$metadata$data_list$theta_hat))
  expect_match(plot_caption(drawn), "Thin lines: 7 posterior draws of g")
  expect_match(plot_caption(drawn), "Marks on the axis: the estimates")
  # No more draws than the fit holds.
  expect_length(unique(plot_layer_data(plot_g(fit, draws = 5000), 2L)$group), 200L)
})

test_that("the estimates, a normal distribution and the truth are laid over the density", {
  fit <- plot_fit(150L)
  data <- fit$metadata$data_list

  with_hist <- plot_g(fit, observed = "histogram", normal = TRUE)
  expect_identical(unname(plot_geoms(with_hist)), c("GeomRect", "GeomRibbon", "GeomLine"))
  bars <- plot_layer_data(with_hist, 1L)
  expect_equal(sum((bars$xmax - bars$xmin) * bars$ymax), 1, tolerance = 1e-8)
  expect_equal(range(c(bars$xmin, bars$xmax)), range(data$theta_hat))
  lines <- plot_layer_data(with_hist, 3L)
  expect_length(unique(lines$group), 2L)
  reml <- getFromNamespace(".bef_reml_normal", "bayesEfron")(data$theta_hat, data$sigma)
  normal <- lines[lines$group == 2L, ]
  expect_equal(normal$y, stats::dnorm(normal$x, reml$mean, reml$sd), tolerance = 1e-10)
  expect_match(plot_caption(with_hist), "estimated by REML")
  expect_match(plot_caption(with_hist), "not an estimate of g")

  given <- plot_g(fit, normal = c(mean = 0.3, sd = 1.2), scale = 10)
  curve <- plot_layer_data(given, 2L)
  curve <- curve[curve$group == 2L, ]
  expect_equal(curve$y, stats::dnorm(curve$x / 10, 0.3, 1.2) / 10, tolerance = 1e-10)
  expect_match(plot_caption(given), "the normal distribution that was given")

  truth <- c(stats::qnorm(seq(0.01, 0.99, length.out = 60), -1, 0.3),
             stats::qnorm(seq(0.01, 0.99, length.out = 90), 1, 0.5))
  both <- plot_g(fit, truth = truth, observed = "histogram")
  expect_identical(unname(plot_geoms(both)), c("GeomRect", "GeomRect", "GeomRibbon", "GeomLine"))
  expect_equal(range(c(plot_layer_data(both, 1L)$xmin, plot_layer_data(both, 1L)$xmax)),
               range(truth))
  expect_match(plot_caption(both), "Gray bars: the true effects. Outlined bars: the estimates")

  density <- function(x) 0.4 * stats::dnorm(x, -1, 0.3) + 0.6 * stats::dnorm(x, 1, 0.5)
  with_function <- plot_g(fit, truth = density)
  curve <- plot_layer_data(with_function, 2L)
  curve <- curve[curve$group == 2L, ]
  expect_equal(curve$y, density(curve$x), tolerance = 1e-10)

  # With the distribution function, samples are step lines that end at one.
  steps <- plot_g(fit, type = "cdf", truth = truth, observed = "histogram", normal = TRUE)
  expect_identical(unname(plot_geoms(steps)), c("GeomRibbon", "GeomLine"))
  lines <- plot_layer_data(steps, 2L)
  expect_length(unique(lines$group), 4L)
  expect_equal(as.numeric(tapply(lines$y, lines$group, max)), rep(1, 4L), tolerance = 1e-6)

  reml_metafor <- function() {
    model <- metafor::rma(yi = data$theta_hat, vi = data$sigma^2, method = "REML")
    c(as.numeric(stats::coef(model)), sqrt(model$tau2))
  }
  skip_if_not_installed("metafor")
  expect_equal(c(reml$mean, reml$sd), reml_metafor(), tolerance = 1e-5)
})

test_that("other estimates of the distribution are drawn as lines with a legend", {
  fit <- plot_fit(12L)
  other <- plot_fit(12L, seed = 7L)
  grid <- fit$metadata$data_list$grid
  outside <- data.frame(x = grid, density = stats::dnorm(grid, 0.5, 1))
  plot <- plot_g(fit, compare = list("Another fit" = other, "Elsewhere" = outside), scale = 2)
  lines <- plot_layer_data(plot, 2L)
  expect_length(unique(lines$group), 3L)
  expect_equal(lines$y[lines$group == 2L], g_summary(other, scale = 2)$mean)
  expect_equal(lines$y[lines$group == 3L], outside$density / 2)
  expect_equal(lines$x[lines$group == 3L], 2 * grid)
  expect_false(identical(plot$theme$legend.position, "none"))
  expect_match(plot_caption(plot), "named in the legend, without bands")
  expect_identical(plot_g(fit)$theme$legend.position, "none")

  # A density becomes a distribution function; a distribution function with
  # jumps has no density.
  as_cdf <- plot_layer_data(plot_g(fit, type = "cdf", compare = list(Elsewhere = outside)), 2L)
  elsewhere <- as_cdf[as_cdf$group == 2L, ]
  expect_equal(elsewhere$y[[nrow(elsewhere)]], 1)
  expect_true(all(diff(elsewhere$y) >= 0))
  jumps <- data.frame(x = c(-1, 0, 1), cdf = c(0.2, 0.7, 1))
  expect_s3_class(plot_g(fit, type = "cdf", compare = list(NPMLE = jumps)), "ggplot")
  expect_error(plot_g(fit, compare = list(NPMLE = jumps)), "use `type = \"cdf\"`",
               class = "bef_invalid_args")
  for (bad in list(list(other), list(a = other, a = other), other, list(a = 1:3),
                   list(a = data.frame(y = 1:3)))) {
    expect_error(plot_g(fit, compare = bad), class = "bef_invalid_args")
  }
})

test_that("the axis shows where g has its mass when the grid is much wider", {
  limits <- getFromNamespace(".bef_plot_g_limits", "bayesEfron")
  engine_of <- function(sd) {
    grid <- seq(-10, 10, length.out = 101L)
    g <- stats::dnorm(grid, 1, sd)
    getFromNamespace(".bef_engine_from_parts", "bayesEfron")(
      log_g = matrix(log(g / sum(g)), nrow = 3L, ncol = 101L, byrow = TRUE),
      grid = grid, theta_hat = c(0, 1, 2), sigma = c(1, 1, 1)
    )
  }
  # g fills the grid: the whole grid with its outer cells.
  expect_equal(limits(engine_of(5), 1), c(-10.1, 10.1))
  # g takes up a small part: that part with a margin, in the units asked for.
  narrow <- limits(engine_of(0.5), 1)
  expect_gt(narrow[[1L]], -3)
  expect_lt(narrow[[2L]], 5)
  expect_lt(narrow[[1L]], 1 - 3 * 0.5)
  expect_gt(narrow[[2L]], 1 + 3 * 0.5)
  expect_equal(limits(engine_of(0.5), 100), 100 * narrow)

  fit <- plot_fit(12L)
  expect_equal(plot_g(fit, xlim = c(-1, 2))$coordinates$limits$x, c(-1, 2))
  for (bad in list(c(2, 1), 1, c(0, NA), "a")) {
    expect_error(plot_g(fit, xlim = bad), "`xlim`", class = "bef_invalid_args")
  }
})

test_that("for every type the caption counts the draws of g that are drawn, and the band and the mean stay the same", {
  fit <- plot_fit(12L, n_draws = 100L)
  for (type in c("density", "cdf", "mass")) {
    base <- plot_g(fit, type = type)
    for (requested in c(7L, 1000L)) {
      drawn <- plot_g(fit, type = type, draws = requested)
      single <- plot_layer_data(drawn, 2L)
      expected <- min(requested, posterior::ndraws(fit$draws))
      expect_length(unique(single$group), expected)
      expect_match(plot_caption(drawn), sprintf("Thin lines: %d posterior draws of g", expected))
      expect_equal(plot_layer_data(drawn, 1L), plot_layer_data(base, 1L))
      expect_equal(plot_layer_data(drawn, 3L)$y, plot_layer_data(base, 2L)$y)
    }
  }
})

test_that("normal components give a smooth density beyond the ends of the grid", {
  fit <- plot_fit(12L, component = "normal")
  data <- fit$metadata$data_list
  line <- plot_layer_data(plot_g(fit), 2L)
  expect_lt(min(line$x), data$grid[[1L]])
  expect_gt(max(line$x), data$grid[[data$L]])
  expect_equal(line$y, g_summary(fit, at = line$x)$mean)
  expect_equal(sum(line$y) * diff(line$x[1:2]), 1, tolerance = 1e-3)
})

test_that("arguments that do not fit the plot are refused", {
  fit <- plot_fit(12L)
  expect_error(plot_g(fit, type = "mass", normal = TRUE), "type = \"density\"",
               class = "bef_invalid_args")
  expect_error(plot_g(fit, type = "mass", observed = "histogram"), class = "bef_invalid_args")
  expect_s3_class(plot_g(fit, type = "mass", observed = "rug", draws = 3), "ggplot")
  expect_error(plot_g(fit, colour = "red"), "has no argument `colour`", class = "bef_invalid_args")
  expect_error(plot_g(fit, type = "histogram"), class = "bef_invalid_args")
  expect_error(plot_g(fit, level = 1.2), class = "bef_invalid_args")
  expect_error(plot_g(fit, scale = 0), class = "bef_invalid_args")
  expect_error(plot_g(fit, draws = -1), class = "bef_invalid_args")
  expect_error(plot_g(fit, normal = c(0, 1)), "c\\(mean = , sd = \\)", class = "bef_invalid_args")
  expect_error(plot_g(fit, truth = "a"), class = "bef_invalid_args")
  expect_error(plot_g(fit, truth = function(x) 1), "one number for each effect value",
               class = "bef_invalid_args")
  expect_error(plot_g(list()), class = "bef_invalid_fit")
})

test_that("the two figures of the README are one call each", {
  krw <- plot_g(krw_fit, observed = "histogram", normal = TRUE, scale = 100)
  expect_identical(unname(plot_geoms(krw)), c("GeomRect", "GeomRibbon", "GeomLine"))
  # The vertical axis is a density per percentage point.
  line <- plot_layer_data(krw, 3L)
  line <- line[line$group == 1L, ]
  expect_equal(sum(line$y) * diff(line$x[1:2]), 1, tolerance = 1e-8)
  expect_equal(diff(line$x[1:2]), 100 * diff(krw_fit$metadata$data_list$grid[1:2]))

  twin <- plot_g(twin_towers_fit, truth = twin_towers$theta_true, normal = TRUE)
  expect_identical(unname(plot_geoms(twin)), c("GeomRect", "GeomRibbon", "GeomLine"))
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(print(krw))
  expect_silent(print(twin))
})
