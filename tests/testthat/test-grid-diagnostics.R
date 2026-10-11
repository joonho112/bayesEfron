# What a fit says about its grid: the share of a site's posterior at each
# end, who set that end, and the spread of g next to its prior.

gd_ns <- function(name) getFromNamespace(name, "bayesEfron")

# An engine for `K` sites on a grid of 101 points from -1 to 1. Every draw
# of g is the same normal shape; `bump` adds mass on the last grid point.
gd_engine <- function(theta_hat, sigma, center = 0, sd = 0.2, bump = 0, tau_c = 0,
                      n_draws = 20L) {
  grid <- seq(-1, 1, length.out = 101L)
  g <- stats::dnorm(grid, center, sd)
  g <- g / sum(g)
  g[101L] <- g[101L] + bump
  g <- g / sum(g)
  gd_ns(".bef_engine_from_parts")(
    log_g = matrix(log(g), nrow = n_draws, ncol = 101L, byrow = TRUE),
    grid = grid, theta_hat = theta_hat, sigma = sigma, tau_c = tau_c
  )
}

gd_grid <- function(support, wide = NULL, method = "paper_realdata") {
  grid <- make_efron_grid(c(-0.5, -0.2, 0, 0.2, 0.5), rep(0.2, 5), support = support, L = 101L)
  grid$grid_method <- method
  if (!is.null(wide)) {
    grid$pilot <- list(wide = wide, narrowed = TRUE)
  }
  grid
}

test_that("sites inside the bulk of g leave both ends alone", {
  engine <- gd_engine(c(-0.3, -0.1, 0, 0.1, 0.3), rep(0.1, 5))
  out <- gd_ns(".bef_grid_diagnostics")(engine, labels = letters[1:5])
  expect_named(out, c("ends", "share", "threshold", "sd_g", "sd_g_prior", "occupied",
                      "tau_c_over_sd_g"))
  expect_named(out$ends, c("end", "value", "set_by", "mass", "site", "label",
                           "estimate_beyond", "sites_over", "limits_at_end", "status"))
  expect_identical(out$ends$end, c("lower", "upper"))
  expect_identical(out$ends$value, c(-1, 1))
  expect_true(all(out$ends$mass < 1e-6))
  expect_identical(out$ends$status, c("ok", "ok"))
  expect_identical(out$ends$sites_over, c(0L, 0L))
  expect_identical(out$ends$limits_at_end, c(0L, 0L))
  expect_true(all(out$ends$estimate_beyond < 0))
  expect_identical(c(out$share, out$threshold), c(0.02, 0.05))
  expect_true(is.na(out$tau_c_over_sd_g))
  # g is a normal shape with standard deviation 0.2 in every draw.
  expect_equal(unname(out$sd_g), rep(0.2, 3L), tolerance = 1e-3)
  expect_true(all(is.na(out$sd_g_prior)))
  # The central 99.9% of that shape: about 3.3 standard deviations each way.
  expect_equal(out$occupied, 2 * 3.3 * 0.2 / 2, tolerance = 0.05)
})

test_that("a site that leans on an end is found, with its share and its label", {
  # The last site has an estimate beyond the upper end, and g has mass there.
  theta_hat <- c(-0.3, -0.1, 0, 0.1, 1.3)
  sigma <- c(0.1, 0.1, 0.1, 0.1, 0.15)
  engine <- gd_engine(theta_hat, sigma, bump = 0.01)
  out <- gd_ns(".bef_grid_diagnostics")(engine, labels = c("a", "b", "c", "d", "far"))
  upper <- out$ends[2L, ]
  expect_identical(upper$site, 5L)
  expect_identical(upper$label, "far")
  expect_gt(upper$mass, 0.9)
  expect_equal(upper$estimate_beyond, (1.3 - 1) / 0.15)
  expect_identical(upper$sites_over, 1L)
  expect_identical(upper$limits_at_end, 1L)
  expect_identical(upper$status, "warning")
  expect_identical(out$ends$status[[1L]], "ok")

  # The share is that of the grid points in the outer 2% of the range.
  weights <- gd_ns(".bef_site_weights")(engine)
  expect_equal(upper$mass, sum(weights[5L, engine$grid >= 0.96 - 1e-9]))
  expect_equal(out$ends$mass[[1L]], max(rowSums(weights[, engine$grid <= -0.96 + 1e-9])))
})

test_that("an end is a warning only where nothing but the caller or the estimates set it", {
  theta_hat <- c(-0.3, -0.1, 0, 0.1, 1.3)
  sigma <- c(0.1, 0.1, 0.1, 0.1, 0.15)
  engine <- gd_engine(theta_hat, sigma, bump = 0.01)
  status <- function(...) gd_ns(".bef_grid_diagnostics")(engine, ...)$ends[, c("set_by", "status")]

  out <- status(grid_method = "paper_realdata")
  expect_identical(out$set_by, c("estimates", "estimates"))
  expect_identical(out$status, c("ok", "warning"))
  out <- status(grid = gd_grid(c(-1, 1), method = "support"))
  expect_identical(out$set_by, c("support", "support"))
  expect_identical(out$status, c("ok", "warning"))

  for (method in c("paper_simulation", "paper_sensitivity")) {
    out <- status(grid_method = method)
    expect_identical(out$set_by, c("true effects", "true effects"))
    expect_identical(out$status, c("ok", "note"))
  }

  # The rule "pilot": an end that it moved in is a note, one that it left
  # where the estimates put it is a warning.
  out <- status(grid = gd_grid(c(-1, 1), wide = c(-3, 3), method = "pilot"))
  expect_identical(out$set_by, c("pilot", "pilot"))
  expect_identical(out$status, c("ok", "note"))
  out <- status(grid = gd_grid(c(-1, 1), wide = c(-3, 1), method = "pilot"))
  expect_identical(out$set_by, c("pilot", "estimates"))
  expect_identical(out$status, c("ok", "warning"))
  out <- status(grid = gd_grid(c(-1, 1), wide = c(-1, 1), method = "pilot"))
  expect_identical(out$set_by, c("estimates", "estimates"))
})

test_that("with normal components the share is read from the mixture and their width is reported", {
  theta_hat <- c(-0.3, -0.1, 0, 0.1, 1.3)
  sigma <- c(0.1, 0.1, 0.1, 0.1, 0.15)
  engine <- gd_engine(theta_hat, sigma, bump = 0.01, tau_c = 0.02)
  out <- gd_ns(".bef_grid_diagnostics")(engine)
  expect_equal(out$tau_c_over_sd_g, 0.02 / out$sd_g[["mean"]])
  expect_gt(out$ends$mass[[2L]], 0.5)
  expect_identical(out$ends$limits_at_end, c(0L, 0L))
  cdf <- gd_ns(".bef_site_cdf")(engine, 0.96)
  expect_equal(out$ends$mass[[2L]], max(1 - cdf))
})

test_that("the prior of the standard deviation of g is simulated without touching R's generator", {
  grid <- make_efron_grid(krw_discrimination, grid_method = "paper_realdata")
  withr::local_seed(4)
  before <- .Random.seed
  prior <- gd_ns(".bef_prior_sd_g")(grid$grid, grid$B)
  expect_identical(.Random.seed, before)
  expect_identical(prior, gd_ns(".bef_prior_sd_g")(grid$grid, grid$B))
  expect_true(all(diff(prior) > 0))
  # The prior draws g toward the uniform distribution on the grid, whose
  # standard deviation is the width over the square root of 12.
  uniform <- diff(grid$support) / sqrt(12)
  expect_gt(prior[[2L]], 0.8 * uniform)
  expect_lt(prior[[2L]], 1.2 * uniform)
  # Normal components add their variance.
  wider <- gd_ns(".bef_prior_sd_g")(grid$grid, grid$B, tau_c = uniform)
  expect_equal(wider, sqrt(prior^2 + uniform^2), tolerance = 1e-7)
})

test_that("diagnose() carries the grid diagnostics of the stored fits", {
  krw <- diagnose(krw_fit)
  expect_named(krw$grid, c("ends", "share", "threshold", "sd_g", "sd_g_prior", "occupied",
                           "tau_c_over_sd_g"))
  expect_identical(krw$grid$ends$status, c("ok", "ok"))
  expect_identical(krw$grid$ends$set_by, c("estimates", "estimates"))
  # On the range of version 0.3.0, which the rule "pilot" leaves for these
  # data, no site leans on an end: the largest shares are about 0.001 and
  # 0.016 (0.002 and 0.013 with the draws of version 0.3.0).
  expect_lt(max(krw$grid$ends$mass), 0.005)
  expect_equal(round(max(diagnose(raudenbush_fit)$grid$ends$mass), 2), 0.02)
  expect_false("grid_ends" %in% krw$sampler_diagnostics_warned)
  # The posterior of the standard deviation of g is far below its prior.
  expect_lt(krw$grid$sd_g[["q95"]], krw$grid$sd_g_prior[["q5"]])

  # The simulated sites on the range that the rule narrowed: sites lean on
  # the lower end, whose estimates lie beyond it, and that is a note.
  twin <- diagnose(twin_towers_fit)
  expect_identical(twin$grid$ends$set_by, c("pilot", "pilot"))
  expect_gt(twin$grid$ends$mass[[1L]], 0.05)
  expect_identical(twin$grid$ends$status[[1L]], "note")
  expect_gt(twin$grid$ends$estimate_beyond[[1L]], 0)
  expect_lt(max(twin$grid$ends$estimate_beyond), 4)
  expect_false("grid_ends" %in% twin$sampler_diagnostics_warned)

  table <- as.data.frame(krw)
  expect_identical(table$diagnostic[[7L]], "grid_ends")
  expect_identical(table$threshold[[7L]], 0.05)
  expect_identical(table$status[[7L]], "ok")
  expect_equal(table$value[[7L]], max(krw$grid$ends$mass))
  expect_identical(as.data.frame(twin)$status[[7L]], "note")

  shown <- format(krw)
  expect_true(any(shown == "Grid:"))
  expect_true(any(grepl("upper end 0.1585: at most 0.1% of a site's posterior in the outer 2%",
                        shown, fixed = TRUE)))
  expect_true(any(grepl("sd of g: posterior 0.017", shown, fixed = TRUE)))
  expect_true(any(grepl("that end was set by the pilot fit on a wider range", format(twin),
                        fixed = TRUE)))
  expect_identical(summary(krw)$grid, krw$grid)

  # A range around the true effects of a simulation is a note as well.
  simulated <- twin
  simulated$grid$ends$set_by <- c("true effects", "true effects")
  expect_true(any(grepl("that end was set from the true effects", format(simulated),
                        fixed = TRUE)))
})

test_that("a new fit stores the diagnostics of its grid and warns about an end that a site leans on", {
  check <- gd_ns(".bef_fit_check_grid")
  fit <- krw_fit
  # A fit whose grid ends where the largest estimate lies: its draws of g on
  # a narrower support.
  narrow <- fit
  data <- narrow$metadata$data_list
  y <- data$theta_hat
  pad <- 0.05 * diff(range(y))
  data$grid <- seq(min(y) - pad, max(y) + pad, length.out = data$L)
  narrow$metadata$data_list <- data
  narrow$metadata$theta_summary$label <- sprintf("firm %d", seq_along(y))

  caught <- NULL
  out <- withCallingHandlers(
    check(narrow),
    bef_warn_grid_ends = function(w) {
      caught <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_s3_class(caught, "bef_warn_grid_ends")
  expect_s3_class(caught, "bef_warning")
  largest <- which.max(y)
  expect_match(conditionMessage(caught), sprintf("Site \"firm %d\" has", largest), fixed = TRUE)
  expect_match(conditionMessage(caught), "at the upper end", fixed = TRUE)
  expect_match(conditionMessage(caught), "compare_efron_grids()", fixed = TRUE)
  expect_identical(caught$ends, attr(out$metadata, "grid_diagnostics", exact = TRUE)$ends)
  expect_true("grid_ends" %in% attr(out$metadata, "sampler_diagnostics_warned", exact = TRUE))
  expect_true("grid_ends" %in% diagnose(out)$sampler_diagnostics_warned)
  expect_identical(as.data.frame(diagnose(out))$status[[7L]], "warning")
  expect_true(any(grepl("Warnings: .*grid_ends", format(out))))
  expect_true(any(grepl("the result for that site depends on where the grid ends",
                        format(diagnose(out)), fixed = TRUE)))

  # The fit on its own range: nothing is said, and the diagnostics are
  # stored all the same.
  expect_no_warning(same <- check(fit))
  stored <- attr(same$metadata, "grid_diagnostics", exact = TRUE)
  expect_identical(stored$ends$status, c("ok", "ok"))
  expect_false("grid_ends" %in% attr(same$metadata, "sampler_diagnostics_warned", exact = TRUE))
  expect_identical(diagnose(same)$grid, stored)
})

test_that("ranges that are too narrow for the small data sets are reported, and their defaults are not", {
  skip_if_no_cmdstan()
  fit_short <- function(...) {
    caught <- list()
    fit <- withCallingHandlers(
      bayes_efron_fit(..., chains = 2L, iter_warmup = 300L, iter_sampling = 300L),
      warning = function(w) {
        caught[[length(caught) + 1L]] <<- w
        invokeRestart("muffleWarning")
      }
    )
    list(fit = fit,
         ends = any(vapply(caught, inherits, logical(1L), "bef_warn_grid_ends")),
         message = paste(vapply(caught, conditionMessage, ""), collapse = " "))
  }
  cases <- list(
    list(data = krw_discrimination, seed = 20250807L, narrow = c(-0.040, 0.084)),
    list(data = raudenbush1985, seed = 1985L, narrow = c(-0.33, 0.55))
  )
  for (case in cases) {
    largest <- which.max(case$data$yi)
    # The default, and the ranges of 0.3.0 with its width and twice it.
    for (args in list(list(), list(grid_method = "paper_realdata"),
                      list(grid_method = "paper_realdata", expansion = 1))) {
      out <- do.call(fit_short, c(list(case$data, seed = case$seed), args))
      expect_false(out$ends)
      expect_lt(max(diagnose(out$fit)$grid$ends$mass), 0.05)
    }
    # A range that stops just beyond the estimates, and one that the largest
    # estimate lies outside.
    for (args in list(list(grid_method = "paper_realdata", expansion = 0.05),
                      list(support = case$narrow))) {
      out <- do.call(fit_short, c(list(case$data, seed = case$seed), args))
      expect_true(out$ends)
      ends <- diagnose(out$fit)$grid$ends
      expect_identical(ends$status[[2L]], "warning")
      expect_gt(ends$mass[[2L]], 0.05)
      expect_match(out$message, sprintf("Site \"%d\" has", ends$site[[2L]]), fixed = TRUE)
      # The site that is named has one of the largest estimates.
      expect_gte(case$data$yi[[ends$site[[2L]]]], stats::quantile(case$data$yi, 0.85))
    }
  }
  # For the callback data on the range that stops just beyond the estimates,
  # the firm with the largest estimate is the one that is named.
  out <- fit_short(krw_discrimination, grid_method = "paper_realdata", expansion = 0.05,
                   seed = 20250807L)
  expect_identical(diagnose(out$fit)$grid$ends$site[[2L]], which.max(krw_discrimination$yi))
})

test_that("a range that the rule \"pilot\" narrowed raises no warning about its ends", {
  skip_if_no_cmdstan()
  # The pilot fit is replaced by an estimate of g that is uniform between
  # -2.46 and 3.13. From it the rule gives the range of `twin_towers_fit`,
  # -3.003 to 3.674, so that only the short fit on that range runs.
  testthat::local_mocked_bindings(
    .bef_pilot_fit = function(theta_hat, sigma, grid, seed, ...) {
      g <- as.numeric(grid >= -2.46 & grid <= 3.13)
      list(g = g / sum(g), diagnostics = list(max_rhat = NA_real_), seconds = 0)
    },
    .package = "bayesEfron"
  )
  caught <- list()
  fit <- withCallingHandlers(
    bayes_efron_fit(twin_towers$yi, twin_towers$sei, grid_method = "pilot", M = 8L,
                    chains = 2L, iter_warmup = 300L, iter_sampling = 300L, seed = 20250627L),
    warning = function(w) {
      caught[[length(caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(vapply(caught, inherits, logical(1L), "bef_warn_grid_ends")))
  ends <- diagnose(fit)$grid$ends
  expect_identical(ends$set_by, c("pilot", "pilot"))
  expect_true(all(ends$status %in% c("ok", "note")))
  expect_false("grid_ends" %in% diagnose(fit)$sampler_diagnostics_warned)
  # Sites do lean on the ends; no estimate is far beyond them.
  expect_gt(max(ends$mass), 0.05)
  expect_lt(max(ends$estimate_beyond), 4)
})
