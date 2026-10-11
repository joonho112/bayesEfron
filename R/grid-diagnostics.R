# What a fit says about its grid: whether a site leans on an end of the
# grid, how the spread of g compares with what the prior alone gives, and
# how much of the grid g occupies.

# The outer share of the range that counts as an end, and the share of a
# site's posterior there above which the site is reported. For the stored
# fits of the callback data and of the teacher expectancy data, on ranges
# half the width of the estimates beyond them, the largest share is 0.001
# and 0.016; for the teacher expectancy data with a quarter of the width
# it is 0.056. On the four ranges that tests/testthat/test-grid-diagnostics.R
# takes as too narrow for those data, fits with the default settings give
# 0.10 to 0.18.
.bef_grid_ends_rule <- function() {
  list(share = 0.02, threshold = 0.05)
}

# Who set each end of a grid. An end that the caller gave or that follows
# from the estimates alone may cut the distribution of the effects. An end
# that the rule "pilot" moved in is where a fit on the wider range put the
# end of g, and a range around the true effects of a simulation holds them
# by construction.
.bef_grid_ends_set_by <- function(grid, grid_method) {
  if (inherits(grid, "bef_grid")) {
    grid_method <- grid$grid_method
  }
  if (identical(grid_method, "support")) {
    return(c("support", "support"))
  }
  if (grid_method %in% c("paper_simulation", "paper_sensitivity")) {
    return(c("true effects", "true effects"))
  }
  if (identical(grid_method, "pilot") && inherits(grid, "bef_grid") &&
      !is.null(grid$pilot)) {
    moved <- c(grid$support[1L] > grid$pilot$wide[1L], grid$support[2L] < grid$pilot$wide[2L])
    return(ifelse(moved, "pilot", "estimates"))
  }
  c("estimates", "estimates")
}

# The prior distribution of the standard deviation of g on a grid, by
# simulation: its 5%, 50% and 95% quantiles. The state of R's generator is
# put back.
.bef_prior_sd_g <- function(grid, B, tau_c = 0, n = 4000L, seed = 1L) {
  old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      suppressWarnings(rm(".Random.seed", envir = globalenv(), inherits = FALSE))
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  B <- matrix(as.numeric(B), nrow = length(grid))
  lambda <- abs(stats::rcauchy(n, 0, 5))
  alpha <- matrix(stats::rnorm(n * ncol(B)), nrow = n) / sqrt(lambda)
  eta <- alpha %*% t(B)
  eta <- eta - apply(eta, 1L, max)
  g <- exp(eta)
  g <- g / rowSums(g)
  center <- (grid[1L] + grid[length(grid)]) / 2
  d <- grid - center
  mean <- as.numeric(g %*% d)
  sd <- sqrt(pmax(as.numeric(g %*% d^2) - mean^2, 0) + tau_c^2)
  stats::quantile(sd, c(0.05, 0.5, 0.95), names = FALSE)
}

# The diagnostics of the grid of a fit. `engine` is that of the fit; `grid`
# is the grid that the fit recorded, or `NULL` for a fit that recorded none,
# in which case `grid_method` names its rule. `level` is the level of the
# site intervals whose limits are looked for at the ends.
.bef_grid_diagnostics <- function(engine, grid = NULL, grid_method = "paper_realdata",
                                  labels = NULL, B = NULL, level = 0.9) {
  rule <- .bef_grid_ends_rule()
  L <- engine$L
  K <- engine$K
  ends <- c(engine$grid[1L], engine$grid[L])
  width <- ends[2L] - ends[1L]
  weights <- .bef_site_weights(engine)
  inner <- c(ends[1L] + rule$share * width, ends[2L] - rule$share * width)
  if (engine$tau_c > 0) {
    cdf <- .bef_site_cdf(engine, inner, weights = weights)
    mass <- cbind(cdf[, 1L], 1 - cdf[, 2L])
    at_end <- c(0L, 0L)
  } else {
    # The grid points in the outer share of the range, with a tolerance for
    # points that fall on its boundary.
    tolerance <- 1e-9 * width
    lower_points <- engine$grid <= inner[1L] + tolerance
    upper_points <- engine$grid >= inner[2L] - tolerance
    mass <- cbind(rowSums(weights[, lower_points, drop = FALSE]),
                  rowSums(weights[, upper_points, drop = FALSE]))
    # A limit of the interval lies on the end point when that point alone
    # holds the tail probability.
    tail <- (1 - level) / 2
    at_end <- c(sum(weights[, 1L] >= tail), sum(weights[, L] > tail))
  }
  site <- c(which.max(mass[, 1L]), which.max(mass[, 2L]))
  largest <- c(mass[site[1L], 1L], mass[site[2L], 2L])
  beyond <- c(
    (ends[1L] - engine$theta_hat[site[1L]]) / engine$sigma[site[1L]],
    (engine$theta_hat[site[2L]] - ends[2L]) / engine$sigma[site[2L]]
  )
  set_by <- .bef_grid_ends_set_by(grid, grid_method)
  over <- largest > rule$threshold
  status <- ifelse(!over, "ok", ifelse(set_by %in% c("support", "estimates"), "warning", "note"))
  site_labels <- .bef_site_labels(labels, K)

  # The spread of g: its posterior, and what the prior alone gives on this
  # grid.
  g <- exp(engine$log_g)
  center <- (ends[1L] + ends[2L]) / 2
  d <- engine$grid - center
  mean_d <- as.numeric(g %*% d)
  sd_draws <- sqrt(pmax(as.numeric(g %*% d^2) - mean_d^2, 0) + engine$tau_c^2)
  sd_posterior <- c(mean = mean(sd_draws),
                    stats::quantile(sd_draws, c(0.05, 0.95), names = FALSE))
  names(sd_posterior) <- c("mean", "q5", "q95")
  sd_prior <- if (is.null(B)) {
    c(q5 = NA_real_, q50 = NA_real_, q95 = NA_real_)
  } else {
    stats::setNames(.bef_prior_sd_g(engine$grid, B, tau_c = engine$tau_c), c("q5", "q50", "q95"))
  }

  # The share of the range between the grid points that hold the central
  # 99.9% of the posterior mean of g.
  cumulative <- cumsum(colMeans(g))
  cumulative <- cumulative / cumulative[L]
  first <- which(cumulative >= 0.0005)[1L]
  last <- which(cumulative >= 0.9995)[1L]
  if (is.na(last)) {
    last <- L
  }

  list(
    ends = data.frame(
      end = c("lower", "upper"),
      value = ends,
      set_by = set_by,
      mass = largest,
      site = site,
      label = site_labels[site],
      estimate_beyond = beyond,
      sites_over = c(sum(mass[, 1L] > rule$threshold), sum(mass[, 2L] > rule$threshold)),
      limits_at_end = as.integer(at_end),
      status = status,
      stringsAsFactors = FALSE
    ),
    share = rule$share,
    threshold = rule$threshold,
    sd_g = sd_posterior,
    sd_g_prior = sd_prior,
    occupied = (engine$grid[last] - engine$grid[first]) / width,
    tau_c_over_sd_g = if (engine$tau_c > 0) engine$tau_c / sd_posterior[["mean"]] else NA_real_
  )
}

# The diagnostics of the grid of a fit object: those stored when the fit was
# made, or computed now for a fit that stored none.
.bef_fit_grid_diagnostics <- function(fit) {
  stored <- .bef_metadata_attr(fit$metadata, "grid_diagnostics")
  if (is.list(stored) && is.data.frame(stored$ends)) {
    return(stored)
  }
  tryCatch(
    .bef_grid_diagnostics(
      .bef_engine(fit),
      grid = .bef_metadata_attr(fit$metadata, "grid"),
      grid_method = fit$metadata$grid_method,
      labels = fit$metadata$theta_summary$label,
      B = fit$metadata$data_list$B
    ),
    error = function(err) NULL
  )
}

# One sentence for each end that a site leans on and that was not set by a
# fit on a wider range.
.bef_grid_ends_sentences <- function(diagnostics) {
  ends <- diagnostics$ends
  flagged <- ends[ends$status == "warning", , drop = FALSE]
  vapply(seq_len(nrow(flagged)), function(i) {
    sprintf(
      "Site %s has %s%% of its posterior distribution in the outer %s%% of the grid at the %s end (%s).",
      .bef_quote_label(flagged$label[[i]]),
      format(round(100 * flagged$mass[[i]], 1), trim = TRUE),
      format(100 * diagnostics$share, trim = TRUE), flagged$end[[i]],
      format(signif(flagged$value[[i]], 3), trim = TRUE)
    )
  }, character(1L))
}

.bef_quote_label <- function(label) {
  sprintf("\"%s\"", label)
}

.bef_warn_grid_ends <- function(diagnostics) {
  sentences <- .bef_grid_ends_sentences(diagnostics)
  if (length(sentences) == 0L) {
    return(invisible(FALSE))
  }
  .bef_warn(
    paste(
      c(sentences,
        "The result for such a site depends on where the grid ends. A wider range may change it; compare_efron_grids() fits several ranges, and diagnose() shows both ends."),
      collapse = " "
    ),
    "bef_warn_grid_ends",
    ends = diagnostics$ends
  )
  invisible(TRUE)
}
