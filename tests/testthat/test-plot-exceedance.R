# plot_exceedance(): the probabilities are those of exceedance_prob(), and
# the count at the cutoff is the count of those numbers.

# The count that the caption of a plot gives.
exceedance_count <- function(plot) {
  as.integer(sub(".* ([0-9]+) of the [0-9]+ sites are at or above it.*", "\\1",
                 plot_caption(plot)))
}

test_that("the points are the probabilities of exceedance_prob()", {
  fit <- plot_fit(12L, labels = TRUE)
  probability <- exceedance_prob(fit)
  plot <- plot_exceedance(fit)
  expect_s3_class(plot, "ggplot")
  expect_identical(unname(plot_geoms(plot)), c("GeomSegment", "GeomPoint"))
  expect_null(plot$labels$title)
  expect_identical(plot$labels$x, "Probability that the effect is above 0")

  by_probability <- order(probability, -seq_along(probability))
  points <- plot_layer_data(plot, 2L)
  expect_equal(points$x, unname(probability[by_probability]))
  expect_identical(plot_axis_labels(plot), names(probability)[by_probability])
  expect_match(plot_caption(plot), "ordered by the probability, the largest at the top")
  expect_no_match(plot_caption(plot), "Dashed line")

  # In the order of the data the first site is at the top.
  unsorted <- plot_exceedance(fit, sort = FALSE)
  expect_equal(plot_layer_data(unsorted, 2L)$x, unname(rev(probability)))
  expect_identical(plot_axis_labels(unsorted), rev(names(probability)))
  expect_no_match(plot_caption(unsorted), "ordered")

  # Another threshold and the other direction; `scale` only changes how the
  # threshold is written.
  below <- plot_exceedance(fit, threshold = 0.5, direction = "below", scale = 100, sort = FALSE)
  expect_equal(plot_layer_data(below, 2L)$x,
               unname(rev(exceedance_prob(fit, threshold = 0.5, direction = "below"))))
  expect_identical(below$labels$x, "Probability that the effect is below 50")
})

test_that("the cutoff is a line, and the caption counts the sites that reach it", {
  fit <- plot_fit(12L, labels = TRUE)
  probability <- exceedance_prob(fit)
  cutoff <- 0.8
  plot <- plot_exceedance(fit, cutoff = cutoff)
  expect_identical(unname(plot_geoms(plot)), c("GeomSegment", "GeomVline", "GeomPoint"))
  expect_equal(plot_layer_data(plot, 2L)$xintercept, cutoff)
  expect_identical(exceedance_count(plot), sum(probability >= cutoff))
  expect_match(plot_caption(plot), "does not control an error rate over many sites")
  # The sites that reach it are filled points.
  points <- plot_layer_data(plot, 3L)
  expect_identical(points$shape == 16, points$x >= cutoff)

  # A probability equal to the cutoff counts.
  exact <- plot_exceedance(fit, cutoff = max(probability[probability < 1 & probability > 0]))
  expect_gte(exceedance_count(exact), 1L)

  stored <- plot_exceedance(bayesEfron::krw_fit, cutoff = 0.95)
  expect_identical(exceedance_count(stored), sum(exceedance_prob(bayesEfron::krw_fit) >= 0.95))
  expect_error(plot_exceedance(fit, cutoff = 1), class = "bef_invalid_args")
  expect_error(plot_exceedance(fit, cutoff = c(0.5, 0.9)), class = "bef_invalid_args")
})

test_that("many sites are drawn as a curve, with a message", {
  fit <- plot_fit(150L)
  probability <- unname(exceedance_prob(fit))
  expect_message(plot <- plot_exceedance(fit, cutoff = 0.9), "Raise `max_sites`",
                 class = "bef_plot_many_sites")
  expect_identical(unname(plot_geoms(plot)), c("GeomHline", "GeomStep"))
  curve <- plot_layer_data(plot, 2L)
  expect_equal(curve$y, sort(probability))
  expect_equal(curve$x, seq_len(150L))
  expect_identical(plot$labels$y, "Probability that the effect is above 0")
  expect_identical(exceedance_count(plot), sum(probability >= 0.9))
  expect_match(plot_caption(plot), "probabilities of the 150 sites, in increasing order")

  expect_message(unsorted <- plot_exceedance(fit, sort = FALSE), class = "bef_plot_many_sites")
  expect_identical(unname(plot_geoms(unsorted)), "GeomPoint")
  expect_equal(plot_layer_data(unsorted, 1L)$y, probability)

  expect_no_message(labeled <- plot_exceedance(fit, max_sites = 150))
  expect_length(plot_axis_labels(labeled), 150L)
})

test_that("wrong arguments are errors", {
  fit <- plot_fit(12L)
  expect_error(plot_exceedance(fit, threshold = c(0, 1)), class = "bef_invalid_args")
  expect_error(plot_exceedance(fit, direction = "over"), class = "bef_invalid_args")
  expect_error(plot_exceedance(fit, sort = "yes"), class = "bef_invalid_args")
  expect_error(plot_exceedance(fit, max_sites = 2.5), class = "bef_invalid_args")
  expect_error(plot_exceedance(fit, sites = 1:3), "no argument", class = "bef_invalid_args")
  expect_error(plot_exceedance(NULL), class = "bef_invalid_fit")
})
