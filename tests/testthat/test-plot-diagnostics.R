# plot_diagnostics(): the values are those of diagnose(), in blocks, and no
# panel is empty.

# The names of the panels of a plot, and the number of points in each.
diagnostic_panels <- function(plot) {
  built <- ggplot2::ggplot_build(plot)
  layout <- built$layout$layout
  points <- which(plot_geoms(plot) == "GeomPoint")
  counts <- table(factor(built$data[[points]]$PANEL, levels = layout$PANEL))
  stats::setNames(as.integer(counts), as.character(layout$panel))
}

test_that("the panels hold the diagnostics of diagnose()", {
  fit <- plot_fit(12L)
  diagnostics <- diagnose(fit)
  plot <- plot_diagnostics(fit)
  expect_s3_class(plot, "ggplot")
  expect_null(plot$labels$title)
  panels <- diagnostic_panels(plot)
  expect_identical(names(panels),
                   c("R-hat", "Bulk ESS", "Tail ESS", "Mass at the ends of the grid", "SD of g"))
  # No panel is empty, with no divergent transition to show.
  expect_true(all(panels > 0L))
  expect_identical(sum(diagnostics$divergences), 0)

  data <- plot$data
  rhat <- data[data$panel == "R-hat", ]
  finite <- diagnostics$rhat[is.finite(diagnostics$rhat)]
  expect_identical(rhat$variable, names(finite))
  expect_equal(rhat$value, unname(finite))
  expect_identical(rhat$flagged, unname(finite > 1.01))
  bulk <- data[data$panel == "Bulk ESS", ]
  expect_equal(bulk$value, unname(diagnostics$ess_bulk[is.finite(diagnostics$ess_bulk)]))
  expect_identical(bulk$flagged, bulk$value < 400)
  tail <- data[data$panel == "Tail ESS", ]
  expect_equal(tail$value, unname(diagnostics$ess_tail[is.finite(diagnostics$ess_tail)]))

  # The worst quantity of each block of more than one is named.
  worst <- rhat[rhat$worst, ]
  by_block <- split(rhat, as.character(rhat$row))
  several <- by_block[vapply(by_block, nrow, integer(1L)) > 1L]
  expect_setequal(as.character(worst$row), names(several))
  expect_equal(worst$value[match(names(several), worst$row)],
               unname(vapply(several, function(block) max(block$value), numeric(1L))))
  low <- bulk[bulk$worst, ]
  expect_equal(low$value[low$row == "log_g"], min(bulk$value[bulk$row == "log_g"]))

  # A block of more than 30 quantities is a box.
  expect_true("GeomBoxplot" %in% plot_geoms(plot))
  expect_identical(unique(as.character(data$row[data$boxed])), "log_g")

  # The grid: the mass at each end with its site, and the two intervals.
  ends <- data[data$panel == "Mass at the ends of the grid", ]
  expect_equal(ends$value, diagnostics$grid$ends$mass)
  expect_identical(ends$variable, sprintf("site %s", diagnostics$grid$ends$label))
  interval <- plot_layer_data(plot, which(plot_geoms(plot) == "GeomSegment"))
  expect_equal(interval$x, unname(c(diagnostics$grid$sd_g_prior[["q5"]], diagnostics$grid$sd_g[["q5"]])))
  expect_equal(interval$xend,
               unname(c(diagnostics$grid$sd_g_prior[["q95"]], diagnostics$grid$sd_g[["q95"]])))

  expect_match(plot_caption(plot), "Dashed lines: 1.01 for R-hat, 400 for the effective sample")
  expect_match(plot_caption(plot), "Dotted line: 1.05")
  expect_match(plot_caption(plot), "Divergent transitions: 0. Iterations at the maximum tree depth: 0")
  expect_match(plot_caption(plot),
               sprintf("Smallest E-BFMI: %s", format(min(diagnostics$ebfmi), digits = 3)), fixed = TRUE)
})

test_that("quantities are put in blocks", {
  blocks <- getFromNamespace(".bef_diagnostic_blocks", "bayesEfron")
  expect_identical(
    blocks(c("lp__", "alpha[1]", "lambda", "mean_g", "var_g", "sd_g", "log_g[3]", "g[2]",
             "theta_rep[2]", "theta_mean[1]", "log_lik[4]", "effective_params",
             "log_marginal_likelihood")),
    c("other", "alpha", "lambda", "summaries of g", "summaries of g", "summaries of g", "log_g",
      "log_g", "sites", "sites", "sites", "other", "other")
  )

  fit <- plot_fit(12L)
  plot <- plot_diagnostics(fit, what = "rhat")
  # From the top: the coefficients first.
  expect_identical(rev(plot_axis_labels(plot)),
                   c("alpha", "lambda", "summaries of g", "log_g", "other"))
  one_row <- plot_diagnostics(fit, what = "rhat", by_block = FALSE)
  expect_identical(plot_axis_labels(one_row), "All quantities")
  expect_identical(sum(one_row$data$worst), 1L)
  expect_identical(one_row$data$variable[one_row$data$worst],
                   names(which.max(diagnose(fit)$rhat)))
})

test_that("`what` chooses the panels", {
  fit <- plot_fit(12L)
  expect_identical(names(diagnostic_panels(plot_diagnostics(fit, what = "rhat"))), "R-hat")
  expect_identical(names(diagnostic_panels(plot_diagnostics(fit, what = c("grid", "ess")))),
                   c("Bulk ESS", "Tail ESS", "Mass at the ends of the grid", "SD of g"))
  grid <- plot_diagnostics(fit, what = "grid")
  expect_identical(names(diagnostic_panels(grid)), c("Mass at the ends of the grid", "SD of g"))
  expect_no_match(plot_caption(grid), "monitored quantity")
  expect_no_match(plot_caption(grid), "Dotted line")
  expect_match(plot_caption(grid), "Dashed lines: 0.05 for the mass at an end of the grid")

  expect_error(plot_diagnostics(fit, what = "divergences"), class = "bef_invalid_args")
  expect_error(plot_diagnostics(fit, what = character()), class = "bef_invalid_args")
  expect_error(plot_diagnostics(fit, by_block = 1), class = "bef_invalid_args")
  expect_error(plot_diagnostics(fit, level = 0.9), "no argument", class = "bef_invalid_args")
  expect_error(plot_diagnostics(list()), class = "bef_invalid_fit")
})

test_that("values beyond a threshold are marked", {
  # A fit of a few short chains, whose diagnostics are poor.
  fit <- readRDS(test_path("_fixtures", "short_fit.rds"))
  diagnostics <- suppressWarnings(diagnose(fit))
  plot <- suppressWarnings(plot_diagnostics(fit))
  data <- plot$data
  expect_true(any(data$flagged))
  expect_identical(data$mark, ifelse(data$flagged, "beyond", "within"))
  expect_true("sites" %in% as.character(data$row))
  expect_match(plot_caption(plot), "Orange triangles: beyond the dashed line")
  expect_true(all(diagnostic_panels(plot) > 0L))
  expect_identical(sum(data$flagged[data$panel == "R-hat"]), sum(diagnostics$rhat > 1.01, na.rm = TRUE))
})

test_that("the three stored fits are drawn", {
  for (fit in list(bayesEfron::krw_fit, bayesEfron::raudenbush_fit, bayesEfron::twin_towers_fit)) {
    plot <- plot_diagnostics(fit)
    expect_true(all(diagnostic_panels(plot) > 0L))
    expect_length(diagnostic_panels(plot), 5L)
  }
})

test_that("a quantity whose R-hat is missing is not drawn, and the caption names it and the failed check", {
  fit <- plot_fit(12L, n_draws = 600L)
  stored <- attr(fit$metadata, "diagnostics")
  stored$rhat["lambda"] <- NA_real_
  attr(fit$metadata, "diagnostics") <- stored
  attr(fit$metadata, "sampler_diagnostics_failed") <- "rhat"
  diagnostics <- diagnose(fit)
  expect_identical(diagnostics$sampler_diagnostics_failed, "rhat")
  plot <- plot_diagnostics(fit)
  rhat <- plot$data[plot$data$panel == "R-hat", ]
  expect_false("lambda" %in% rhat$variable)
  finite <- diagnostics$rhat[is.finite(diagnostics$rhat)]
  expect_equal(rhat$value, unname(finite))
  expect_match(plot_caption(plot), "Sampler checks failed: R-hat")
  expect_match(plot_caption(plot), "R-hat unavailable for lambda")
  expect_no_match(plot_caption(plot_diagnostics(fit, what = "ess")), "R-hat unavailable")
})

test_that("in a run of fewer than 400 draws the marks follow the values, and the caption says the checks were skipped", {
  fit <- plot_fit(12L, n_draws = 100L)
  diagnostics <- diagnose(fit)
  expect_setequal(diagnostics$diagnostic_skipped,
                  c("rhat_check", "ess_bulk_check", "ess_tail_check"))
  expect_length(diagnostics$sampler_diagnostics_failed, 0L)
  plot <- plot_diagnostics(fit)
  bulk <- plot$data[plot$data$panel == "Bulk ESS", ]
  expect_equal(bulk$value, unname(diagnostics$ess_bulk))
  expect_identical(bulk$flagged, bulk$value < 400)
  expect_match(plot_caption(plot), "Sampler checks skipped: R-hat, bulk ESS, tail ESS")
  expect_match(plot_caption(plot), "threshold marks describe the values")
  expect_match(plot_caption(plot), "fewer than 400 sampling draws")
  expect_no_match(plot_caption(plot), "Sampler checks failed")
  grid <- plot_diagnostics(fit, what = "grid")
  expect_match(plot_caption(grid), "Sampler checks skipped")
  expect_no_match(plot_caption(grid), "threshold marks describe the values")
  expect_no_match(plot_caption(grid), "fewer than 400 sampling draws")
})

test_that("long lists of unavailable quantities are abbreviated", {
  fit <- plot_fit(12L, n_draws = 600L)
  stored <- attr(fit$metadata, "diagnostics")
  absent <- names(stored$rhat)[seq_len(8L)]
  stored$rhat[absent] <- NA_real_
  attr(fit$metadata, "diagnostics") <- stored
  plot <- plot_diagnostics(fit, what = "rhat")
  expect_match(plot_caption(plot), "R-hat unavailable for")
  expect_match(plot_caption(plot), "(2 more)", fixed = TRUE)
  expect_true(all(is.finite(plot$data$value)))
})
