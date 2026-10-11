#' Plot the exceedance probabilities of the site effects
#'
#' @description
#' For each site, `plot_exceedance()` shows the exceedance probability,
#' i.e., the posterior probability that the effect is above (or below) a
#' threshold, as computed by [exceedance_prob()]. If `cutoff` is specified,
#' a line marks that probability and the caption gives the number of sites
#' that reach it.
#'
#' @details
#' Up to `max_sites` sites are shown as labeled points. If there are more,
#' the probabilities are shown as a curve in increasing order or, when
#' `sort = FALSE`, as points in the order of the data.
#'
#' `threshold` is specified in the units of the estimates, as in
#' [exceedance_prob()]. `scale` only changes the units in which the
#' threshold appears in the axis label: with `threshold = 0.02` and
#' `scale = 100`, the label reads "above 2".
#'
#' The number of sites at or above `cutoff` is based on the posterior
#' probability of each site separately. It is not the result of a
#' procedure that controls an error rate over many sites, as the caption
#' also notes.
#'
#' @inheritParams plot_g
#' @inheritParams plot_sites
#' @inheritParams site_summary
#' @param threshold Threshold for the effects, a number in the units of the
#'   estimates.
#' @param direction One of `"above"` (the default) for the probability
#'   that the effect is greater than `threshold`, or `"below"` for the
#'   probability that it is not.
#' @param sort If `TRUE` (the default), the sites are ordered by their
#'   probabilities. If `FALSE`, they are in the order of the data.
#' @param cutoff An optional probability between 0 and 1, which is marked
#'   by a dashed line.
#' @param scale Positive number by which `threshold` is multiplied in the
#'   axis label. See Details.
#'
#' @return A `ggplot` object.
#'
#' @seealso [exceedance_prob()], [plot_sites()], `vignette("plots")`
#'
#' @examples
#' # The probability that the callback gap of a firm is positive, and the
#' # firms for which it is at least 0.95.
#' plot_exceedance(krw_fit, cutoff = 0.95)
#'
#' # The probability that the gap exceeds 2 percentage points.
#' plot_exceedance(krw_fit, threshold = 0.02, scale = 100)
#'
#' # 1,500 sites: a curve.
#' plot_exceedance(twin_towers_fit, cutoff = 0.9)
#'
#' @export
plot_exceedance <- function(x,
                            threshold = 0,
                            direction = c("above", "below"),
                            sort = TRUE,
                            cutoff = NULL,
                            max_sites = 100,
                            scale = 1,
                            ...) {
  .bef_check_dots_empty(list(...), "plot_exceedance")
  x <- .bef_plot_fit(x)
  direction <- .bef_validate_method_choice(direction, c("above", "below"), arg = "direction")
  sort <- .bef_plot_flag(sort, "sort")
  cutoff <- .bef_plot_level(cutoff, "cutoff", null_ok = TRUE)
  max_sites <- .bef_plot_count(max_sites, "max_sites", lower = 1L)
  scale <- .bef_validate_scale(scale)

  probability <- exceedance_prob(x, threshold = threshold, direction = direction)
  K <- length(probability)
  table <- data.frame(site = seq_len(K), label = names(probability),
                      probability = unname(probability), stringsAsFactors = FALSE)
  table$reached <- if (is.null(cutoff)) TRUE else table$probability >= cutoff
  labeled <- K <= max_sites
  if (!labeled) {
    .bef_inform_many_sites(K, max_sites, "Raise `max_sites` to label them.")
  }
  event <- sprintf("Probability that the effect is %s %s", direction,
                   format(threshold * scale))

  plot <- if (labeled) {
    .bef_plot_exceedance_points(table, sort, cutoff, event)
  } else {
    .bef_plot_exceedance_curve(table, sort, cutoff, event)
  }
  plot + ggplot2::labs(caption = .bef_plot_exceedance_caption(table, labeled, sort, cutoff))
}

# Labeled points, the first site or the largest probability at the top.
.bef_plot_exceedance_points <- function(table, sort, cutoff, event) {
  colors <- .bef_plot_colors()
  key <- if (sort) table$probability else -table$site
  table <- table[order(key, -table$site), , drop = FALSE]
  table$position <- seq_len(nrow(table))
  table$mark <- ifelse(table$reached, "at or above", "below")
  plot <- ggplot2::ggplot(table, ggplot2::aes(x = .data$probability, y = .data$position)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$probability, yend = .data$position),
                          color = "grey85", linewidth = 0.3)
  if (!is.null(cutoff)) {
    plot <- plot +
      ggplot2::geom_vline(xintercept = cutoff, linetype = "dashed", color = "grey35")
  }
  plot +
    ggplot2::geom_point(ggplot2::aes(shape = .data$mark, color = .data$mark), size = 1.8) +
    ggplot2::scale_shape_manual(values = c("at or above" = 16, "below" = 1), guide = "none") +
    ggplot2::scale_color_manual(
      values = c("at or above" = colors[["fit"]], "below" = colors[["data"]]), guide = "none"
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25),
                                expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::scale_y_continuous(breaks = table$position,
                                labels = .bef_short_labels(table$label),
                                expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::labs(x = event, y = NULL) +
    theme_bef() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = .bef_plot_label_size(nrow(table)))
    )
}

# Many sites: the probabilities in increasing order as a curve, or in the
# order of the data as points.
.bef_plot_exceedance_curve <- function(table, sort, cutoff, event) {
  colors <- .bef_plot_colors()
  if (sort) {
    table <- table[order(table$probability, table$site), , drop = FALSE]
  }
  table$position <- seq_len(nrow(table))
  plot <- ggplot2::ggplot(table, ggplot2::aes(x = .data$position, y = .data$probability))
  if (!is.null(cutoff)) {
    plot <- plot +
      ggplot2::geom_hline(yintercept = cutoff, linetype = "dashed", color = "grey35")
  }
  plot <- if (sort) {
    plot + ggplot2::geom_step(color = colors[["fit"]], linewidth = 0.6, direction = "mid")
  } else {
    plot + ggplot2::geom_point(color = colors[["fit"]], size = 0.6,
                               alpha = .bef_plot_alpha(nrow(table)) + 0.25)
  }
  plot +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25),
                                expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::labs(
      x = if (sort) "Sites, in the order of the probability" else "Sites, in the order of the data",
      y = event
    ) +
    theme_bef()
}

.bef_plot_exceedance_caption <- function(table, labeled, sort, cutoff) {
  K <- nrow(table)
  .bef_caption(
    if (labeled) {
      "Each point is the posterior probability for one site."
    } else if (sort) {
      sprintf("The posterior probabilities of the %d sites, in increasing order.", K)
    } else {
      "Each point is the posterior probability for one site."
    },
    if (labeled && sort) "Sites are ordered by the probability, the largest at the top.",
    if (!is.null(cutoff)) {
      sprintf(
        paste("Dashed line: %s. %d of the %d sites are at or above it%s; that count does not",
              "control an error rate over many sites."),
        format(cutoff), sum(table$reached), K,
        if (labeled) " (filled points)" else ""
      )
    }
  )
}
