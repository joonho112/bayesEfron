# plot_recovery(): the intervals are those of confint(), and the shares are
# counted from them.

# True effects for the sites of a fit: its posterior means moved, so that
# some intervals hold them and some do not.
recovery_truth <- function(fit) {
  table <- site_summary(fit)
  withr::with_seed(5L, table$mean + stats::rnorm(nrow(table), 0, 1.2 * table$sd))
}

test_that("the intervals are those of confint(), and the misses are marked", {
  fit <- plot_fit(12L, labels = TRUE)
  truth <- recovery_truth(fit)
  limits <- confint(fit, level = 0.9)
  inside <- truth >= limits$lower & truth <= limits$upper
  expect_true(any(inside) && any(!inside))

  plot <- plot_recovery(fit, truth)
  expect_s3_class(plot, "ggplot")
  expect_identical(unname(plot_geoms(plot)), c("GeomSegment", "GeomPoint", "GeomPoint"))
  expect_null(plot$labels$title)
  by_mean <- order(coef(fit))
  lines <- plot_layer_data(plot, 1L)
  expect_equal(lines$x, limits$lower[by_mean])
  expect_equal(lines$xend, limits$upper[by_mean])
  expect_equal(plot_layer_data(plot, 2L)$x, unname(coef(fit)[by_mean]))
  marks <- plot_layer_data(plot, 3L)
  expect_equal(marks$x, truth[by_mean])
  # A true effect outside its interval is a triangle, inside a cross.
  expect_identical(marks$shape == 24, !inside[by_mean])
  expect_identical(plot$data$mark, ifelse(inside[by_mean], "inside", "outside"))
  expect_match(plot_caption(plot), sprintf("outside it \\(%d of 12\\)", sum(!inside)))
  expect_match(plot_caption(plot), "they are not the coverage over repeated samples")
  expect_no_match(plot_caption(plot), "chosen at random")

  wide <- plot_recovery(fit, truth, level = 0.5, scale = 10)
  half <- confint(fit, level = 0.5)
  expect_equal(plot_layer_data(wide, 1L)$x, 10 * half$lower[by_mean])
  expect_equal(plot_layer_data(wide, 3L)$x, 10 * truth[by_mean])
  expect_match(plot_caption(wide), "50% posterior intervals of the 12 sites")
})

test_that("`n` sites are chosen at random, and `seed` fixes the choice", {
  fit <- plot_fit(150L)
  truth <- recovery_truth(fit)
  limits <- confint(fit)
  share <- mean(truth >= limits$lower & truth <= limits$upper)

  withr::local_seed(99L)
  before <- stats::runif(1L)
  withr::local_seed(99L)
  first <- plot_recovery(fit, truth, n = 40, seed = 7)
  # The state of R's generator is as it was.
  expect_identical(stats::runif(1L), before)
  second <- plot_recovery(fit, truth, n = 40, seed = 7)
  other <- plot_recovery(fit, truth, n = 40, seed = 8)
  expect_identical(nrow(first$data), 40L)
  expect_identical(first$data$site, second$data$site)
  expect_false(identical(sort(first$data$site), sort(other$data$site)))
  expect_identical(anyDuplicated(first$data$site), 0L)
  expect_match(plot_caption(first), "40 of the 150 sites, chosen at random")
  # The share in the caption is that of all the sites.
  expect_match(plot_caption(first),
               sprintf("Of all 150 sites, %s have their true effect inside",
                       paste0(format(100 * round(share, 3L)), "%")), fixed = TRUE)

  all_sites <- plot_recovery(fit, truth, n = 150)
  expect_identical(sort(all_sites$data$site), 1:150)
  sampler <- getFromNamespace(".bef_sample_sites", "bayesEfron")
  expect_identical(sampler(10L, 20L), 1:10)
  expect_length(sampler(100L, 20L, seed = 1L), 20L)
})

test_that("the calibration counts the true effects inside the intervals of each level", {
  fit <- plot_fit(150L)
  truth <- recovery_truth(fit)
  plot <- plot_recovery(fit, truth, type = "calibration")
  expect_identical(unname(plot_geoms(plot)), c("GeomAbline", "GeomLine", "GeomPoint"))
  expect_identical(plot$labels$x, "Level of the posterior interval")
  points <- plot_layer_data(plot, 3L)
  expect_equal(points$x, seq(0.05, 0.95, by = 0.05))
  for (level in c(0.5, 0.9)) {
    limits <- confint(fit, level = level)
    expect_equal(points$y[which.min(abs(points$x - level))],
                 mean(truth >= limits$lower & truth <= limits$upper))
  }
  expect_true(all(diff(points$y) >= 0))
  expect_match(plot_caption(plot), "they are not the coverage over repeated samples")

  # With normal components the intervals are those of that model.
  normal <- plot_fit(12L, component = "normal")
  truth <- recovery_truth(normal)
  limits <- confint(normal, level = 0.8)
  shares <- plot_layer_data(plot_recovery(normal, truth, type = "calibration"), 3L)
  expect_equal(shares$y[which.min(abs(shares$x - 0.8))],
               mean(truth >= limits$lower & truth <= limits$upper))
})

test_that("the stored fit of the simulated sites is drawn, and its share is the count", {
  truth <- bayesEfron::twin_towers$theta_true
  fit <- bayesEfron::twin_towers_fit
  points <- plot_layer_data(plot_recovery(fit, truth, type = "calibration"), 3L)
  limits <- confint(fit, level = 0.9)
  expect_equal(points$y[which.min(abs(points$x - 0.9))],
               mean(truth >= limits$lower & truth <= limits$upper))
})

test_that("the scatter sets the posterior means against the true effects", {
  fit <- plot_fit(12L)
  truth <- recovery_truth(fit)
  table <- site_summary(fit)
  plot <- plot_recovery(fit, truth, type = "scatter")
  expect_identical(unname(plot_geoms(plot)), c("GeomAbline", "GeomPoint"))
  points <- plot_layer_data(plot, 2L)
  expect_equal(points$x, truth)
  expect_equal(points$y, table$mean)
  expect_identical(plot$labels$x, "True effect")
  rmse <- function(values) format(signif(sqrt(mean((values - truth)^2)), 3L))
  expect_match(plot_caption(plot),
               sprintf("%s for the posterior means and %s for the estimates",
                       rmse(table$mean), rmse(table$estimate)), fixed = TRUE)
  expect_equal(plot_layer_data(plot_recovery(fit, truth, type = "scatter", scale = 100), 2L)$y,
               100 * table$mean)
})

test_that("wrong arguments are errors", {
  fit <- plot_fit(12L)
  truth <- recovery_truth(fit)
  expect_error(plot_recovery(fit), "true effect of every site", class = "bef_invalid_args")
  expect_error(plot_recovery(fit, truth[-1L]), class = "bef_invalid_args")
  expect_error(plot_recovery(fit, c(truth[-1L], NA)), class = "bef_invalid_args")
  expect_error(plot_recovery(fit, truth, type = "coverage"), class = "bef_invalid_args")
  expect_error(plot_recovery(fit, truth, n = 0), class = "bef_invalid_args")
  expect_error(plot_recovery(fit, truth, level = 1), class = "bef_invalid_args")
  expect_error(plot_recovery(fit, truth, seed = "a"), class = "bef_invalid_args")
  expect_error(plot_recovery(fit, truth, sites = 1:3), "no argument", class = "bef_invalid_args")
})
