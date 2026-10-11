#' Plot sampler diagnostics and grid diagnostics
#'
#' @description
#' `plot_diagnostics()` plots the diagnostics returned by [diagnose()]:
#' R-hat and the bulk and tail effective sample sizes (ESS) of the monitored
#' quantities against their thresholds, and two diagnostics of the grid.
#' The caption gives the number of divergent transitions, the number of
#' iterations that reached the maximum tree depth, and the smallest E-BFMI.
#'
#' @details
#' By default the monitored quantities are grouped into blocks: the spline
#' coefficients `alpha`, the precision `lambda`, the mean, variance and
#' standard deviation of \eqn{g}, `log_g` at the grid points, the
#' site-level quantities (if the fit stores them, see `site_draws` in
#' [bayes_efron_fit()]), and all other quantities. A block of up to 30
#' quantities is shown as points. A larger block is shown as a box plot,
#' with its worst value as a point. A name next to a block identifies its
#' worst quantity, i.e., the one with the largest R-hat or the smallest
#' effective sample size.
#'
#' Dashed lines mark the thresholds used by [diagnose()]: 1.01 for R-hat
#' and 400 for the effective sample sizes. A dotted line marks an R-hat of
#' 1.05, above which the draws should not be used.
#'
#' Sampler diagnostics are computed from all iterations of the run, even
#' if the fit was thinned (`thin`) and stores only some of them. The two
#' grid panels are computed from the stored draws.
#'
#' Sampler checks that failed or were skipped are listed in the caption.
#' For a short run (fewer than 400 draws), [diagnose()] skips the checks
#' of R-hat and of the effective sample sizes. Their values are still
#' plotted against the thresholds, but a value beyond a threshold then
#' does not mean that a check failed. Values that are not available are
#' not plotted, and the caption names the quantities concerned (a long
#' list is abbreviated). This includes a quantity with an infinite R-hat:
#' its value is stored as `NA`, and the caption also reports it as a
#' failed check.
#'
#' If `what` includes `"grid"`, two more panels show the grid diagnostics.
#' The first gives, for each end of the grid (its lower and upper limit),
#' the largest posterior mass that any site has in the outer 2% of the
#' range, and names that site. A value above the dashed line (0.05) means
#' that the site has substantial posterior mass at that limit. The second
#' panel compares the posterior interval of the standard deviation of
#' \eqn{g} with its prior interval on the same grid. See the Details of
#' [diagnose()] for how to read both panels.
#'
#' @inheritParams plot_g
#' @param what Character vector with the diagnostics to plot: one or more
#'   of `"rhat"`, `"ess"` and `"grid"`. By default all three are plotted.
#' @param by_block If `TRUE` (the default), the quantities are grouped into
#'   blocks as described in Details. If `FALSE`, they are shown in a single
#'   row.
#'
#' @return A `ggplot` object.
#'
#' @seealso [diagnose()], [plot_mcmc()], `vignette("diagnostics")`
#'
#' @examples
#' plot_diagnostics(krw_fit)
#'
#' # R-hat alone, without blocks.
#' plot_diagnostics(twin_towers_fit, what = "rhat", by_block = FALSE)
#'
#' @export
plot_diagnostics <- function(x, what = c("rhat", "ess", "grid"), by_block = TRUE, ...) {
  .bef_check_dots_empty(list(...), "plot_diagnostics")
  x <- .bef_plot_fit(x)
  choices <- c("rhat", "ess", "grid")
  if (!is.character(what) || length(what) == 0L || anyNA(what) || !all(what %in% choices)) {
    .bef_abort_invalid_args(
      "`what` must hold one or more of: \"rhat\", \"ess\", \"grid\".", arg = "what"
    )
  }
  what <- choices[choices %in% what]
  by_block <- .bef_plot_flag(by_block, "by_block")

  diagnostics <- diagnose(x)
  rule <- .bef_plot_diagnostics_rule()
  values <- list()
  lines <- list()
  if ("rhat" %in% what) {
    values$rhat <- .bef_plot_diagnostic_values(diagnostics$rhat, "R-hat", by_block,
                                              worst = "largest", threshold = rule$rhat[[1L]])
    lines$rhat <- data.frame(panel = "R-hat", value = rule$rhat, kind = c("first", "second"))
  }
  if ("ess" %in% what) {
    values$bulk <- .bef_plot_diagnostic_values(diagnostics$ess_bulk, "Bulk ESS", by_block,
                                               worst = "smallest", threshold = rule$ess)
    values$tail <- .bef_plot_diagnostic_values(diagnostics$ess_tail, "Tail ESS", by_block,
                                               worst = "smallest", threshold = rule$ess)
    lines$ess <- data.frame(panel = c("Bulk ESS", "Tail ESS"), value = rule$ess, kind = "first")
  }
  grid <- if ("grid" %in% what) diagnostics$grid
  intervals <- NULL
  if (!is.null(grid)) {
    ends_panel <- "Mass at the ends of the grid"
    values$ends <- data.frame(
      panel = ends_panel,
      row = ifelse(grid$ends$end == "lower", "Lower end", "Upper end"),
      variable = sprintf("site %s", grid$ends$label), value = grid$ends$mass,
      flagged = grid$ends$mass > grid$threshold, boxed = FALSE, worst = TRUE,
      stringsAsFactors = FALSE
    )
    lines$ends <- data.frame(panel = ends_panel, value = grid$threshold, kind = "first")
    intervals <- data.frame(
      panel = "SD of g", row = c("Prior", "Posterior"),
      value = c(grid$sd_g_prior[["q50"]], grid$sd_g[["mean"]]),
      lower = c(grid$sd_g_prior[["q5"]], grid$sd_g[["q5"]]),
      upper = c(grid$sd_g_prior[["q95"]], grid$sd_g[["q95"]]),
      stringsAsFactors = FALSE
    )
    values$sd <- data.frame(
      panel = "SD of g", row = intervals$row, variable = "", value = intervals$value,
      flagged = FALSE, boxed = FALSE, worst = FALSE, stringsAsFactors = FALSE
    )
  }
  values <- do.call(rbind, values)
  lines <- do.call(rbind, lines)
  if (is.null(values) || nrow(values) == 0L) {
    .bef_abort_invalid_args(
      "The fit holds none of the diagnostics that `what` asks for.", arg = "what"
    )
  }

  panels <- unique(values$panel)
  rows <- c("Posterior", "Prior", "Upper end", "Lower end",
            rev(c(.bef_diagnostic_block_names(), "All quantities")))
  set_levels <- function(data) {
    data$panel <- factor(data$panel, levels = panels)
    if (!is.null(data$row)) {
      data$row <- factor(data$row, levels = rows)
    }
    data
  }
  values <- set_levels(values)
  values$mark <- ifelse(values$flagged, "beyond", "within")
  colors <- .bef_plot_colors()

  # The blank layer comes first so that the rows of every panel are in the
  # order of the blocks whatever the layers after it hold.
  plot <- ggplot2::ggplot(values, ggplot2::aes(x = .data$value, y = .data$row)) +
    ggplot2::geom_blank()
  if (!is.null(lines)) {
    plot <- plot +
      ggplot2::geom_vline(data = set_levels(lines),
                          ggplot2::aes(xintercept = .data$value, linetype = .data$kind),
                          color = "grey35") +
      ggplot2::scale_linetype_manual(values = c(first = "dashed", second = "dotted"),
                                     guide = "none")
  }
  boxed <- values[values$boxed, , drop = FALSE]
  if (nrow(boxed) > 0L) {
    plot <- plot +
      ggplot2::geom_boxplot(data = boxed, orientation = "y", width = 0.5, outlier.size = 0.5,
                            color = colors[["fit"]], fill = NA, linewidth = 0.4)
  }
  if (!is.null(intervals)) {
    plot <- plot +
      ggplot2::geom_segment(data = set_levels(intervals),
                            ggplot2::aes(x = .data$lower, xend = .data$upper, yend = .data$row),
                            color = colors[["fit"]], linewidth = 0.6)
  }
  drawn <- values[!values$boxed | values$worst, , drop = FALSE]
  named <- values[values$worst & nzchar(values$variable), , drop = FALSE]
  plot +
    ggplot2::geom_point(data = drawn, ggplot2::aes(color = .data$mark, shape = .data$mark),
                        size = 1.8, alpha = 0.8) +
    ggplot2::geom_text(data = named, ggplot2::aes(label = .data$variable), vjust = -0.9,
                       hjust = "inward", size = 2.6, color = "grey20") +
    ggplot2::scale_color_manual(
      values = c(within = colors[["fit"]], beyond = colors[["highlight"]]), guide = "none"
    ) +
    ggplot2::scale_shape_manual(values = c(within = 16, beyond = 17), guide = "none") +
    # Few breaks, and no label that would run into its neighbor: the panels
    # are narrow, and the lines of the thresholds are named in the caption.
    ggplot2::scale_x_continuous(n.breaks = 4,
                                guide = ggplot2::guide_axis(check.overlap = TRUE)) +
    ggplot2::facet_wrap(ggplot2::vars(.data$panel), scales = "free",
                        ncol = if (length(panels) == 4L) 2L else min(3L, length(panels))) +
    ggplot2::labs(
      x = NULL, y = NULL,
      caption = .bef_plot_diagnostics_caption(diagnostics, what, grid, any(values$boxed),
                                              any(values$flagged))
    ) +
    theme_bef() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

# The lines of the plot: the thresholds of diagnose(), and the R-hat beyond
# which draws should not be used. A block of more than `box` quantities is
# drawn as a box.
.bef_plot_diagnostics_rule <- function() {
  list(rhat = c(1.01, 1.05), ess = 400, box = 30L)
}

.bef_diagnostic_block_names <- function() {
  c("alpha", "lambda", "summaries of g", "log_g", "sites", "other")
}

# The block of each monitored quantity.
.bef_diagnostic_blocks <- function(variables) {
  base <- sub("[[].*$", "", variables)
  block <- rep("other", length(variables))
  block[base == "alpha"] <- "alpha"
  block[base == "lambda"] <- "lambda"
  block[base %in% c("mean_g", "var_g", "sd_g")] <- "summaries of g"
  block[base %in% c("log_g", "g", "log_w")] <- "log_g"
  block[grepl("^theta", base) | base == "log_lik"] <- "sites"
  block
}

# One diagnostic of the monitored quantities as rows of the plot data: the
# value of each quantity, whether it is beyond the threshold, whether its
# block is drawn as a box, and whether it is the worst of its block.
.bef_plot_diagnostic_values <- function(values, panel, by_block, worst, threshold) {
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(NULL)
  }
  larger_is_worse <- identical(worst, "largest")
  out <- data.frame(
    panel = panel,
    row = if (by_block) .bef_diagnostic_blocks(names(values)) else "All quantities",
    variable = names(values), value = unname(values),
    flagged = if (larger_is_worse) values > threshold else values < threshold,
    stringsAsFactors = FALSE
  )
  size <- table(out$row)
  out$boxed <- as.integer(size[out$row]) > .bef_plot_diagnostics_rule()$box
  out$worst <- FALSE
  for (row in names(size)) {
    at <- which(out$row == row)
    if (length(at) > 1L) {
      out$worst[at[if (larger_is_worse) which.max(out$value[at]) else which.min(out$value[at])]] <-
        TRUE
    }
  }
  out
}

.bef_plot_diagnostics_caption <- function(diagnostics, what, grid, boxed, flagged) {
  rule <- .bef_plot_diagnostics_rule()
  thresholds <- c(
    if ("rhat" %in% what) sprintf("%s for R-hat", format(rule$rhat[[1L]])),
    if ("ess" %in% what) sprintf("%s for the effective sample sizes", format(rule$ess)),
    if (!is.null(grid)) sprintf("%s for the mass at an end of the grid", format(grid$threshold))
  )
  .bef_caption(
    if (any(c("rhat", "ess") %in% what)) {
      c("Each point is a monitored quantity.",
        if (boxed) {
          sprintf("A box stands for a block of more than %d, with its worst value as a point.",
                  rule$box)
        },
        "A name is that of the worst quantity of its block.")
    },
    sprintf("Dashed lines: %s.", paste(thresholds, collapse = ", ")),
    if ("rhat" %in% what) sprintf("Dotted line: %s.", format(rule$rhat[[2L]])),
    if (flagged) "Orange triangles: beyond the dashed line.",
    .bef_plot_diagnostics_status(diagnostics, what),
    if (!is.null(grid)) {
      c(sprintf("Grid: the largest posterior mass that a site has in the outer %s of the range",
                .bef_percent(grid$share)),
        "at each end, and 90% intervals of the standard deviation of g under the prior and the",
        "posterior.")
    },
    sprintf("Divergent transitions: %d.", as.integer(sum(diagnostics$divergences))),
    sprintf("Iterations at the maximum tree depth: %d.",
            as.integer(sum(diagnostics$max_treedepth))),
    sprintf("Smallest E-BFMI: %s.", format(min(diagnostics$ebfmi), digits = 3))
  )
}

# Preserve diagnose()'s status when a value cannot be plotted or a threshold
# check was skipped. Missing quantity lists are bounded so that a fit with
# many constant site summaries still has a readable caption.
.bef_plot_diagnostics_status <- function(diagnostics, what) {
  labels <- c(rhat = "R-hat", ess_bulk = "bulk ESS", ess_tail = "tail ESS",
              divergences = "divergent transitions", max_treedepth = "maximum tree depth",
              ebfmi = "E-BFMI", grid_ends = "mass at the ends of the grid",
              sampler_diagnostics = "divergences, maximum tree depth and E-BFMI")
  describe <- function(fields) {
    fields <- sub("_check$", "", fields)
    text <- unname(labels[fields])
    text[is.na(text)] <- fields[is.na(text)]
    paste(text, collapse = ", ")
  }
  failed <- diagnostics$sampler_diagnostics_failed
  skipped <- diagnostics$diagnostic_skipped
  fields <- c(if ("rhat" %in% what) "rhat",
              if ("ess" %in% what) c("ess_bulk", "ess_tail"))
  missing <- unlist(lapply(fields, function(field) {
    values <- diagnostics[[field]]
    absent <- names(values)[!is.finite(values)]
    if (length(absent) == 0L) return(NULL)
    shown <- utils::head(absent, 6L)
    rest <- length(absent) - length(shown)
    sprintf("%s unavailable for %s%s.", labels[[field]], paste(shown, collapse = ", "),
            if (rest > 0L) sprintf(" (%d more)", rest) else "")
  }), use.names = FALSE)
  c(
    if (length(failed) > 0L) sprintf("Sampler checks failed: %s.", describe(failed)),
    if (length(skipped) > 0L) sprintf("Sampler checks skipped: %s.", describe(skipped)),
    if (any(paste0(fields, "_check") %in% skipped))
      sprintf("Their threshold marks describe the values; the run has fewer than %d sampling draws.",
              .bef_diagnostic_thresholds()$min_draws_for_check),
    missing
  )
}
