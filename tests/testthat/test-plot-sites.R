# plot_sites(): the intervals are those of confint(), the sites are those
# that were asked for, and many sites are drawn without labels.

test_that("the intervals and the points are those of the accessors", {
  fit <- plot_fit(12L, labels = TRUE)
  plot <- plot_sites(fit)
  expect_s3_class(plot, "ggplot")
  expect_identical(
    unname(plot_geoms(plot)),
    c("GeomVline", "GeomSegment", "GeomSegment", "GeomPoint", "GeomPoint")
  )
  expect_null(plot$labels$title)
  expect_identical(plot$labels$x, "Effect")

  limits <- confint(fit, level = 0.9)
  by_mean <- order(coef(fit))
  outer <- plot_layer_data(plot, 2L)
  expect_equal(outer$x, limits$lower[by_mean])
  expect_equal(outer$xend, limits$upper[by_mean])
  expect_equal(outer$y, seq_len(12L))
  inner <- plot_layer_data(plot, 3L)
  half <- confint(fit, level = 0.5)
  expect_equal(inner$x, half$lower[by_mean])
  expect_equal(inner$xend, half$upper[by_mean])
  expect_equal(plot_layer_data(plot, 4L)$x, unname(fit$metadata$data_list$theta_hat[by_mean]))
  expect_equal(plot_layer_data(plot, 5L)$x, unname(coef(fit)[by_mean]))
  expect_equal(plot_layer_data(plot, 1L)$xintercept, mean(fit$posterior$mean_g))

  # The labels of the fit are on the axis, the largest mean at the top.
  expect_identical(plot_axis_labels(plot), names(coef(fit))[by_mean])
  expect_match(plot_caption(plot), "Thick and thin lines: 50% and 90% posterior intervals")
  expect_match(plot_caption(plot), "Open circles: the estimates")
  expect_match(plot_caption(plot), "Dashed line: posterior mean of the mean of g")
  expect_match(plot_caption(plot), "ordered by posterior mean, the largest at the top")

  # Other levels, and other units.
  wide <- plot_sites(fit, level = 0.95, inner_level = NULL, scale = 100)
  expect_identical(unname(plot_geoms(wide)), c("GeomVline", "GeomSegment", "GeomPoint", "GeomPoint"))
  expect_equal(plot_layer_data(wide, 2L)$x, 100 * confint(fit, level = 0.95)$lower[by_mean])
  expect_equal(plot_layer_data(wide, 1L)$xintercept, 100 * mean(fit$posterior$mean_g))
  expect_match(plot_caption(wide), "Lines: 95% posterior intervals")
})

test_that("the order of the sites follows `sort_by`", {
  fit <- plot_fit(12L, labels = TRUE)
  table <- site_summary(fit, threshold = 0)
  drawn <- function(...) plot_sites(fit, ...)$data$site

  expect_identical(drawn(), order(table$mean))
  expect_identical(drawn(sort_by = "observed"), order(table$estimate))
  # The most precise estimate is at the top, and so is the first site.
  expect_identical(drawn(sort_by = "sigma"), order(-table$se))
  expect_identical(drawn(sort_by = "none"), 12:1)
  expect_identical(drawn(sort_by = "prob"), order(table$prob_above, -table$site))
  expect_match(plot_caption(plot_sites(fit, sort_by = "sigma")),
               "ordered by standard error, the smallest at the top")
  expect_no_match(plot_caption(plot_sites(fit, sort_by = "none")), "ordered")
})

test_that("sites are chosen by position, by label, or as the extremes", {
  fit <- plot_fit(12L, labels = TRUE)
  labels <- names(coef(fit))

  by_position <- plot_sites(fit, sites = c(3, 7, 1), sort_by = "none")
  expect_identical(by_position$data$site, c(1L, 7L, 3L))
  expect_identical(plot_axis_labels(by_position), labels[c(1L, 7L, 3L)])
  by_label <- plot_sites(fit, sites = labels[c(3L, 7L, 1L)], sort_by = "none")
  expect_identical(by_label$data$site, c(1L, 7L, 3L))

  extremes <- plot_sites(fit, n_extreme = 3)
  by_mean <- order(coef(fit))
  expect_identical(extremes$data$site, by_mean[c(1:3, 10:12)])
  expect_identical(as.character(extremes$data$panel), rep(c("Lowest 3", "Highest 3"), each = 3L))
  # The panel of the highest sites is the upper one.
  expect_identical(levels(extremes$data$panel), c("Highest 3", "Lowest 3"))
  expect_identical(plot_axis_labels(extremes, panel = 1L), labels[by_mean[10:12]])
  expect_identical(plot_axis_labels(extremes, panel = 2L), labels[by_mean[1:3]])

  # In the order of the data, the first of the chosen sites is at the top.
  in_order <- plot_sites(fit, n_extreme = 3, sort_by = "none")
  expect_identical(in_order$data$site,
                   c(sort(by_mean[1:3], decreasing = TRUE), sort(by_mean[10:12], decreasing = TRUE)))

  expect_error(plot_sites(fit, sites = 1:3, n_extreme = 2), class = "bef_invalid_args")
  expect_error(plot_sites(fit, n_extreme = 7), "at most 6", class = "bef_invalid_args")
  expect_error(plot_sites(fit, sites = "no such site"), "names no site",
               class = "bef_invalid_args")
  expect_error(plot_sites(fit, sites = 13), class = "bef_invalid_args")
})

test_that("many sites are drawn as lines without labels, with a message", {
  fit <- plot_fit(150L)
  expect_message(plot <- plot_sites(fit, max_sites = 100), class = "bef_plot_many_sites")
  expect_length(plot_axis_labels(plot), 0L)
  # The intervals, the reference line in front of them, and the means.
  expect_identical(unname(plot_geoms(plot)), c("GeomSegment", "GeomVline", "GeomPoint"))
  expect_identical(nrow(plot_layer_data(plot, 1L)), 150L)
  expect_equal(plot_layer_data(plot, 1L)$x, confint(fit)$lower[order(coef(fit))])
  expect_match(plot_caption(plot), "Each line is the 90% posterior interval of one site")
  expect_no_match(plot_caption(plot), "Open circles")

  # Up to `max_sites` sites keep their labels, and a choice of sites does.
  expect_no_message(labelled <- plot_sites(fit, max_sites = 150))
  expect_length(plot_axis_labels(labelled), 150L)
  expect_no_message(chosen <- plot_sites(fit, n_extreme = 10))
  expect_length(plot_axis_labels(chosen, panel = 1L), 10L)
  expect_identical(chosen$data$site, order(coef(fit))[c(1:10, 141:150)])
  expect_message(plot_sites(fit, n_extreme = 60), "120 sites", class = "bef_plot_many_sites")
})

test_that("highlighted sites, true effects and the reference line are drawn as asked", {
  fit <- plot_fit(12L, labels = TRUE)
  labels <- names(coef(fit))
  colors <- getFromNamespace(".bef_plot_colors", "bayesEfron")()

  marked <- plot_sites(fit, highlight = labels[c(2L, 5L)], sort_by = "none")
  points <- plot_layer_data(marked, 5L)
  expect_identical(points$colour[12:1],
                   ifelse(seq_len(12L) %in% c(2L, 5L), colors[["highlight"]], colors[["fit"]]))
  expect_match(plot_caption(marked), "Orange: the highlighted sites")
  expect_identical(plot_sites(fit, highlight = c(2, 5), sort_by = "none")$data$highlight,
                   marked$data$highlight)

  truth <- seq(-1, 1, length.out = 12L)
  with_truth <- plot_sites(fit, truth = truth, sort_by = "none", scale = 10)
  expect_identical(unname(plot_geoms(with_truth)),
                   c("GeomVline", "GeomSegment", "GeomSegment", "GeomPoint", "GeomPoint",
                     "GeomPoint"))
  expect_equal(plot_layer_data(with_truth, 5L)$x, 10 * rev(truth))
  expect_match(plot_caption(with_truth), "Crosses: the true effects")
  expect_error(plot_sites(fit, truth = truth[-1L]), "each of the 12 sites",
               class = "bef_invalid_args")

  zero <- plot_sites(fit, ref = "zero", show_observed = FALSE)
  expect_identical(unname(plot_geoms(zero)),
                   c("GeomVline", "GeomSegment", "GeomSegment", "GeomPoint"))
  expect_equal(plot_layer_data(zero, 1L)$xintercept, 0)
  expect_match(plot_caption(zero), "Dashed line: zero")
  expect_no_match(plot_caption(zero), "Open circles")
  none <- plot_sites(fit, ref = "none")
  expect_false("GeomVline" %in% plot_geoms(none))
  expect_no_match(plot_caption(none), "Dashed line")
})

test_that("wrong arguments are errors", {
  fit <- plot_fit(12L)
  expect_error(plot_sites(fit, level = 0.5, inner_level = 0.5), "smaller than `level`",
               class = "bef_invalid_args")
  expect_error(plot_sites(fit, level = 1.2), class = "bef_invalid_args")
  expect_error(plot_sites(fit, sort_by = "rank"), class = "bef_invalid_args")
  expect_error(plot_sites(fit, max_sites = 0), class = "bef_invalid_args")
  expect_error(plot_sites(fit, show_observed = NA), class = "bef_invalid_args")
  expect_error(plot_sites(fit, colour = "red"), "no argument", class = "bef_invalid_args")
  expect_error(plot_sites(list()), class = "bef_invalid_fit")
})

test_that("long labels are shortened on the axis", {
  fit <- plot_fit(12L, labels = TRUE)
  short <- getFromNamespace(".bef_short_labels", "bayesEfron")
  expect_identical(short(c("Short", strrep("a", 40L))), c("Short", paste0(strrep("a", 25L), "...")))
  expect_true(all(nchar(plot_axis_labels(plot_sites(fit))) <= 28L))
})
