# The table of site effects: one row for each site, with its label, the
# estimate and standard error that were given, and the mean, standard
# deviation, interval limits and mode of its posterior distribution.
#
# The posterior distribution of theta[i] averages over the draws of g (Lee
# and Sui 2025, Section 4.3 and Appendix A), and every column is computed
# from that distribution by the functions of R/posterior-engine.R.
#
# The columns `hpdi_lower`, `hpdi_upper` and `map` repeat `lower`, `upper`
# and `mode` under the names of version 0.3.0. They stay for version 0.4 and
# are removed in 0.5.

.bef_site_table_columns <- function() {
  c("site", "label", "estimate", "se", "mean", "sd", "lower", "upper", "mode",
    "hpdi_lower", "hpdi_upper", "map")
}

# The level of the interval kept in the table of a fit.
.bef_site_table_level <- function() {
  0.9
}

.bef_theta_summary <- function(engine, labels = NULL, level = .bef_site_table_level()) {
  K <- engine$K
  weights <- .bef_site_weights(engine)
  moments <- .bef_site_moments(engine)
  limits <- .bef_site_quantile(engine, .bef_interval_probs(level), weights = weights)
  mode <- .bef_site_mode(engine, weights = weights)
  out <- data.frame(
    site = seq_len(K),
    label = .bef_site_labels(labels, K),
    estimate = engine$theta_hat,
    se = engine$sigma,
    mean = moments$mean,
    sd = sqrt(moments$var),
    lower = limits[, 1L],
    upper = limits[, 2L],
    mode = mode,
    hpdi_lower = limits[, 1L],
    hpdi_upper = limits[, 2L],
    map = mode,
    row.names = seq_len(K),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  attr(out, "level") <- level
  out
}

# The labels of the sites: those given, or the site numbers as text.
.bef_site_labels <- function(labels, K) {
  if (is.null(labels)) {
    return(as.character(seq_len(K)))
  }
  as.character(labels)
}

# The mode of the averaged posterior of each site. With point masses it is
# the grid point with the largest probability, the smallest one when
# several are equal, as in the Stan program. With normal components it is
# found by .bef_normal_mixture_mode().
.bef_site_mode <- function(engine, weights = NULL, sites = NULL, n = 1001L) {
  sites <- .bef_engine_sites(engine, sites)
  if (is.null(weights)) {
    weights <- .bef_site_weights(engine, sites = sites)
  }
  if (engine$tau_c == 0) {
    return(engine$grid[max.col(weights, ties.method = "first")])
  }
  comp <- .bef_engine_components(engine, sites)
  vapply(seq_along(sites), function(k) {
    .bef_normal_mixture_mode(
      .bef_normal_mixture_parts(comp$mean[k, ], comp$sd[k], weights[k, ]), n = n
    )
  }, numeric(1L))
}

# The table of a fit with all the columns of this version. A table that was
# built with the columns of version 0.3.0 gets the new names, and one built
# with the new names gets the old ones; the label, the estimate and the
# standard error are filled in from the fit when they are missing.
.bef_site_table <- function(fit) {
  table <- fit$metadata$theta_summary
  K <- nrow(table)
  data <- fit$metadata$data_list
  pairs <- c(lower = "hpdi_lower", upper = "hpdi_upper", mode = "map")
  for (new in names(pairs)) {
    old <- pairs[[new]]
    if (is.null(table[[new]]) && !is.null(table[[old]])) table[[new]] <- table[[old]]
    if (is.null(table[[old]]) && !is.null(table[[new]])) table[[old]] <- table[[new]]
  }
  if (is.null(table$label)) table$label <- as.character(table$site)
  if (is.null(table$estimate)) table$estimate <- as.numeric(data$theta_hat)
  if (is.null(table$se)) table$se <- as.numeric(data$sigma)
  level <- attr(table, "level", exact = TRUE)
  table <- table[, .bef_site_table_columns(), drop = FALSE]
  attr(table, "level") <- if (is.null(level)) .bef_site_table_level() else level
  table
}

# The table with the interval at another level.
.bef_site_table_at_level <- function(fit, level) {
  table <- .bef_site_table(fit)
  if (isTRUE(all.equal(level, attr(table, "level", exact = TRUE)))) {
    return(table)
  }
  engine <- .bef_engine(fit)
  limits <- .bef_site_quantile(engine, .bef_interval_probs(level))
  table$lower <- table$hpdi_lower <- limits[, 1L]
  table$upper <- table$hpdi_upper <- limits[, 2L]
  attr(table, "level") <- level
  table
}

# Do the labels differ from the site numbers?
.bef_has_site_labels <- function(table) {
  !identical(table$label, as.character(table$site))
}

# The names of the columns changed in version 0.4; say so once in an
# interactive session. Scripts and documents are not interrupted.
.bef_note_site_table_names <- function(interactive_fun = interactive) {
  if (!isTRUE(interactive_fun()) || isTRUE(.bef_notes[["site_table_names"]])) {
    return(invisible(FALSE))
  }
  .bef_notes[["site_table_names"]] <- TRUE
  message(
    "The interval limits and the mode of the site table are now in the columns ",
    "`lower`, `upper` and `mode`. The columns `hpdi_lower`, `hpdi_upper` and ",
    "`map` hold the same values and will be removed in version 0.5."
  )
  invisible(TRUE)
}

.bef_notes <- new.env(parent = emptyenv())
