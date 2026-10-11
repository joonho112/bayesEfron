#' Plot the Pareto k values of leave-one-out cross-validation
#'
#' @description
#' `plot_loo()` plots the Pareto k value of each site from
#' [loo.bef_fit()], together with the threshold above which the
#' leave-one-out estimate for a site is not reliable. The threshold depends
#' on the number of posterior draws and is computed as in the loo package.
#' Sites above the threshold are labeled.
#'
#' @details
#' For \eqn{S} posterior draws, the threshold is
#' \eqn{\min(1 - 1 / \log_{10} S, 0.7)}{min(1 - 1 / log10(S), 0.7)}
#' (Vehtari et al., 2024). This is 0.7 for about 2,200 draws or more, and
#' 0.656 for the 800 draws of the fits stored in the package.
#'
#' Removing a site that is above the threshold would change the posterior
#' distribution considerably, so the leave-one-out estimate for that site
#' should not be trusted. See `vignette("diagnostics")` for a discussion
#' of this check.
#'
#' The Pareto k values can take some time to compute when there are many
#' sites. To avoid computing them again, pass a result of [loo.bef_fit()]
#' as `loo`. This function requires the loo package.
#'
#' @inheritParams plot_g
#' @param loo An optional `loo` object, the result of [loo.bef_fit()] for
#'   `x`. If `NULL` (the default), it is computed.
#' @param label If `TRUE` (the default), the sites above the threshold are
#'   labeled.
#'
#' @return A `ggplot` object.
#'
#' @references
#' Vehtari, A., Simpson, D., Gelman, A., Yao, Y. and Gabry, J. (2024).
#' Pareto smoothed importance sampling. *Journal of Machine Learning
#' Research*, 25(72), 1-58.
#'
#' @seealso [loo.bef_fit()], [compare_efron_grids()], `vignette("plots")`
#'
#' @examplesIf requireNamespace("loo", quietly = TRUE)
#' plot_loo(krw_fit)
#'
#' # Reuse a result that was computed earlier.
#' result <- loo::loo(raudenbush_fit)
#' plot_loo(raudenbush_fit, loo = result)
#'
#' @export
plot_loo <- function(x, loo = NULL, label = TRUE, ...) {
  .bef_check_dots_empty(list(...), "plot_loo")
  x <- .bef_plot_fit(x)
  label <- .bef_plot_flag(label, "label")
  .bef_require_loo()
  table <- site_summary(x)
  K <- nrow(table)
  if (is.null(loo)) {
    loo <- loo::loo(x)
  }
  k <- if (inherits(loo, "loo")) loo$diagnostics$pareto_k
  if (!is.numeric(k) || length(k) != K) {
    .bef_abort_invalid_args(
      sprintf("`loo` must be the result of `loo()` for this fit, with a Pareto k for each of the %d sites.",
              K),
      arg = "loo"
    )
  }
  draws <- attr(loo, "dims", exact = TRUE)[[1L]]
  threshold <- .bef_pareto_k_threshold(draws)
  data <- data.frame(site = table$site, label = table$label, pareto_k = as.numeric(k),
                     stringsAsFactors = FALSE)
  data$mark <- ifelse(data$pareto_k > threshold, "above", "below")
  above <- data[data$mark == "above", , drop = FALSE]
  colors <- .bef_plot_colors()

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$site, y = .data$pareto_k)) +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", color = "grey35") +
    ggplot2::geom_point(ggplot2::aes(color = .data$mark, shape = .data$mark),
                        size = if (K > 100L) 0.9 else 1.8, alpha = if (K > 100L) 0.7 else 1) +
    ggplot2::scale_color_manual(
      values = c(below = colors[["fit"]], above = colors[["highlight"]]), guide = "none"
    ) +
    ggplot2::scale_shape_manual(values = c(below = 16, above = 17), guide = "none")
  if (label && nrow(above) > 0L) {
    plot <- plot +
      ggplot2::geom_text(data = above, ggplot2::aes(label = .bef_short_labels(.data$label)),
                         vjust = -0.8, size = 2.8, color = "grey20", check_overlap = TRUE)
  }
  plot +
    ggplot2::labs(
      x = "Site, in the order of the data", y = "Pareto k",
      caption = .bef_caption(
        "Each point is the Pareto k of one site in leave-one-out cross-validation.",
        sprintf("Dashed line: %s, above which the estimate for a site is not reliable with %s draws.",
                format(round(threshold, 3L)), format(draws, big.mark = ",")),
        if (nrow(above) == 0L) {
          "No site is above it."
        } else {
          sprintf("%d of the %d sites %s above it (orange triangles).", nrow(above), K,
                  if (nrow(above) == 1L) "is" else "are")
        }
      )
    ) +
    theme_bef()
}

# The Pareto k above which a leave-one-out estimate from `draws` draws is not
# reliable, as the loo package computes it.
.bef_pareto_k_threshold <- function(draws) {
  min(1 - 1 / log10(draws), 0.7)
}
