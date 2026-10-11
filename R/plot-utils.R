# What the plot functions share.
#
# Conventions of every plot_*() function and of the plot() methods:
#
#   1. The result is a ggplot object, so that layers, scales and a theme can
#      be added to it. Nothing is drawn until it is printed.
#   2. There is no title. A caption says what the lines, bands and points
#      are, with the level of every interval.
#   3. The numbers that are drawn are those of the accessors: g_summary(),
#      site_summary(), exceedance_prob(), diagnose() and loo(). A plot
#      computes no summary of its own.
#   4. Arguments come in this order: the fit `x`; what is drawn (`type` or
#      `style`, and the sites); the levels (`level`, `inner_level`); what is
#      laid over it (`observed`, `normal`, `truth`, `highlight`, ...);
#      `max_sites`; `scale`; and `...`, which must be empty.
#   5. `scale` multiplies every value on the axis of the effects and divides
#      every density, so that `scale = 100` turns proportions into
#      percentage points.
#   6. Colors have fixed roles (Okabe and Ito): this fit is blue and solid,
#      a reference model is vermilion and dashed, the data and the true
#      effects are gray, other fits are bluish green, reddish purple, orange
#      and sky blue, and what is highlighted is orange. A difference in color
#      is repeated in the line type or the point shape. A continuous color
#      (the standard error) is viridis.
#   7. With more than `max_sites` sites and no choice of sites, a plot of
#      sites draws them without labels and says how to choose some.

.bef_plot_colors <- function() {
  c(
    fit = "#0072B2", reference = "#D55E00", data = "#7F7F7F", highlight = "#E69F00",
    other1 = "#009E73", other2 = "#CC79A7", other3 = "#E69F00", other4 = "#56B4E9"
  )
}

# The theme of the plots. It is added when a plot is built, so a theme that
# the caller adds afterwards replaces it.
theme_bef <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(hjust = 0, color = "grey30"),
      plot.caption.position = "plot",
      strip.text = ggplot2::element_text(hjust = 0)
    )
}

# The fit of a plot function, with the checks that the methods make.
.bef_plot_fit <- function(x, what) {
  .bef_require_inherits(x, "bef_fit", "x", "bef_invalid_fit")
  .bef_prepare_fit(x)
}

# A level, or `NULL` where `NULL` is allowed.
.bef_plot_level <- function(level, arg = "level", null_ok = FALSE) {
  if (is.null(level) && isTRUE(null_ok)) {
    return(NULL)
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    .bef_abort_invalid_args(
      sprintf("`%s` must be a single number between 0 and 1%s.", arg,
              if (isTRUE(null_ok)) ", or `NULL`" else ""),
      arg = arg
    )
  }
  level
}

.bef_plot_flag <- function(value, arg) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    .bef_abort_invalid_args(sprintf("`%s` must be TRUE or FALSE.", arg), arg = arg)
  }
  value
}

.bef_plot_count <- function(value, arg, lower = 0L) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value != floor(value) || value < lower) {
    .bef_abort_invalid_args(
      sprintf("`%s` must be a single whole number of at least %d.", arg, as.integer(lower)),
      arg = arg
    )
  }
  as.integer(value)
}

# A level as text: 0.9 is "90%".
.bef_percent <- function(level) {
  paste0(format(100 * level, trim = TRUE), "%")
}

# The sentences of a caption, wrapped to the width of a plot.
.bef_caption <- function(..., width = 100L) {
  text <- paste(c(...), collapse = " ")
  paste(strwrap(text, width = width), collapse = "\n")
}

# Labels short enough for an axis.
.bef_short_labels <- function(labels, width = 28L) {
  labels <- as.character(labels)
  long <- nchar(labels) > width
  labels[long] <- paste0(substr(labels[long], 1L, width - 3L), "...")
  labels
}

# The positions of the sites that a plot is asked for: by position or by
# label, in the order given.
.bef_plot_sites <- function(sites, labels, arg = "sites") {
  K <- length(labels)
  if (is.null(sites)) {
    return(NULL)
  }
  if (is.character(sites)) {
    at <- match(sites, labels)
    if (anyNA(at)) {
      .bef_abort_invalid_args(
        sprintf("`%s` names no site: %s.", arg,
                paste(sprintf("\"%s\"", sites[is.na(at)]), collapse = ", ")),
        arg = arg
      )
    }
    return(unique(at))
  }
  if (!is.numeric(sites) || length(sites) == 0L || anyNA(sites) ||
      any(sites != floor(sites)) || any(sites < 1L) || any(sites > K)) {
    .bef_abort_invalid_args(
      sprintf("`%s` must hold positions between 1 and %d, or labels of sites.", arg, K),
      arg = arg
    )
  }
  unique(as.integer(sites))
}

# The true effects of the sites of simulated data: one number for each site.
.bef_plot_truth <- function(truth, K, arg = "truth") {
  if (is.null(truth)) {
    return(NULL)
  }
  if (!is.numeric(truth) || length(truth) != K || any(!is.finite(truth))) {
    .bef_abort_invalid_args(
      sprintf("`%s` must hold one finite number for each of the %d sites.", arg, K),
      arg = arg
    )
  }
  as.numeric(truth)
}

# Said when a plot leaves the labels of the sites out because there are too
# many of them.
.bef_inform_many_sites <- function(K, max_sites,
                                   how = "Choose fewer with `sites` or `n_extreme`, or raise `max_sites`.") {
  rlang::inform(
    sprintf(
      "There are %d sites to draw, more than `max_sites` (%d), so they are drawn without labels. %s",
      K, max_sites, how
    ),
    class = "bef_plot_many_sites"
  )
}

# The normal random-effects model fitted by restricted maximum likelihood:
# the mean and the standard deviation of the effects. The restricted log
# likelihood is maximized over the variance of the effects, which is not
# allowed below zero; the mean is the weighted mean at that variance.
.bef_reml_normal <- function(theta_hat, sigma) {
  v <- sigma^2
  objective <- function(tau2) {
    w <- 1 / (v + tau2)
    mu <- sum(w * theta_hat) / sum(w)
    0.5 * (sum(log(v + tau2)) + log(sum(w)) + sum(w * (theta_hat - mu)^2))
  }
  upper <- 10 * max(stats::var(theta_hat), max(v))
  found <- stats::optimize(objective, c(0, upper), tol = 1e-12)
  tau2 <- if (objective(0) <= found$objective) 0 else found$minimum
  w <- 1 / (v + tau2)
  list(mean = sum(w * theta_hat) / sum(w), sd = sqrt(tau2))
}

# The mean and the standard deviation of a normal reference: `TRUE` for the
# REML estimates, or the two numbers.
.bef_plot_normal <- function(normal, theta_hat, sigma) {
  if (isFALSE(normal) || is.null(normal)) {
    return(NULL)
  }
  if (isTRUE(normal)) {
    out <- .bef_reml_normal(theta_hat, sigma)
    out$source <- "REML"
    return(out)
  }
  ok <- is.numeric(normal) && length(normal) == 2L && all(is.finite(normal)) &&
    setequal(names(normal), c("mean", "sd")) && normal[["sd"]] >= 0
  if (!ok) {
    .bef_abort_invalid_args(
      "`normal` must be FALSE, TRUE, or a mean and a standard deviation given as `c(mean = , sd = )`.",
      arg = "normal"
    )
  }
  list(mean = normal[["mean"]], sd = normal[["sd"]], source = "given")
}
