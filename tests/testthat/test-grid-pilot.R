# The range rule "pilot": the range of "paper_realdata", narrowed to where a
# short fit puts g. The arithmetic of the rule is tested with a stand-in for
# the fit; the tests at the end run CmdStan.

pilot_ns <- function(name) getFromNamespace(name, "bayesEfron")

# A stand-in for the pilot fit that returns the probabilities `g_fun(grid)`.
local_pilot <- function(g_fun, record = new.env(parent = emptyenv()),
                        env = parent.frame()) {
  record$calls <- list()
  testthat::local_mocked_bindings(
    .bef_pilot_fit = function(theta_hat, sigma, grid, seed, parallel_chains = NULL, ...) {
      record$calls[[length(record$calls) + 1L]] <- list(
        theta_hat = theta_hat, sigma = sigma, grid = grid, seed = seed,
        parallel_chains = parallel_chains
      )
      g <- g_fun(grid)
      list(
        g = g / sum(g),
        diagnostics = list(max_rhat = 1, min_ess_bulk = 1000, divergences = 0L,
                           max_treedepth = 0L),
        seconds = 0
      )
    },
    .package = "bayesEfron",
    .env = env
  )
  record
}

pilot_flat <- function(grid) rep(1, length(grid))
pilot_bell <- function(center, sd) function(grid) stats::dnorm(grid, center, sd)

pilot_data <- function() {
  list(
    y = c(-1.9, -1.2, -0.6, -0.3, 0, 0.1, 0.4, 0.8, 1.3, 2.1),
    s = c(0.5, 0.4, 0.3, 0.2, 0.2, 0.2, 0.3, 0.3, 0.4, 0.5)
  )
}

test_that("a pilot estimate that fills the range leaves the grid of \"paper_realdata\"", {
  d <- pilot_data()
  record <- local_pilot(pilot_flat)
  grid <- make_efron_grid(d$y, d$s, grid_method = "pilot", seed = 7L)
  old <- make_efron_grid(d$y, d$s, grid_method = "paper_realdata")

  expect_s3_class(grid, "bef_grid")
  expect_identical(grid$grid_method, "pilot")
  expect_identical(grid$grid, old$grid)
  expect_identical(grid$B, old$B)
  expect_identical(grid$support, old$support)
  expect_identical(grid$expansion, 0.5)
  expect_false(grid$pilot$narrowed)
  expect_identical(grid$pilot$share, 1)
  expect_identical(grid$pilot$wide, old$support)
  expect_identical(grid$pilot$seed, 7L)
  expect_identical(grid$pilot$margin, 0.1)
  expect_identical(grid$pilot$protected, integer())
  expect_identical(grid$pilot[c("chains", "iter_warmup", "iter_sampling")],
                   list(chains = 4L, iter_warmup = 500L, iter_sampling = 500L))
  expect_named(grid$pilot$diagnostics,
               c("max_rhat", "min_ess_bulk", "divergences", "max_treedepth"))

  # The pilot fit is run on the wide range with 101 points, whatever `L` is.
  expect_length(record$calls, 1L)
  expect_equal(record$calls[[1L]]$grid, seq(old$support[1L], old$support[2L], length.out = 101L))
  expect_identical(record$calls[[1L]]$seed, 7L)
  make_efron_grid(d$y, d$s, grid_method = "pilot", seed = 7L, L = 201L, M = 8L)
  expect_length(record$calls[[2L]]$grid, 101L)

  # Other rules have no such element and run no fit.
  expect_null(old$pilot)
  expect_false("pilot" %in% names(old))
  expect_length(record$calls, 2L)
})

test_that("the rule is the default of make_efron_grid() and of bayes_efron_fit()", {
  expect_identical(formals(make_efron_grid)$grid_method, "pilot")
  expect_identical(formals(bayes_efron_fit)$grid_method, "pilot")
  d <- pilot_data()
  record <- local_pilot(pilot_flat)
  grid <- make_efron_grid(d$y, d$s, seed = 3L)
  expect_identical(grid$grid_method, "pilot")
  expect_length(record$calls, 1L)
  # A range that is given takes the place of the rule, and no fit is run.
  given <- make_efron_grid(d$y, d$s, support = c(-3, 3))
  expect_identical(given$grid_method, "support")
  expect_length(record$calls, 1L)
})

test_that("a concentrated pilot estimate narrows the range to its central 99.9% and a margin", {
  d <- pilot_data()
  local_pilot(pilot_bell(0.2, 0.5))
  wide <- make_efron_grid(d$y, d$s, grid_method = "paper_realdata")$support
  pilot_grid <- seq(wide[1L], wide[2L], length.out = 101L)
  g <- stats::dnorm(pilot_grid, 0.2, 0.5)
  cumulative <- cumsum(g / sum(g))
  core <- c(pilot_grid[which(cumulative >= 0.0005)[1L]],
            pilot_grid[which(cumulative >= 0.9995)[1L]])
  expect_lt(core[1L], 0.2 - 3 * 0.5)
  expect_gt(core[2L], 0.2 + 3 * 0.5)

  grid <- make_efron_grid(d$y, d$s, grid_method = "pilot", seed = 1L)
  expect_true(grid$pilot$narrowed)
  expect_equal(grid$pilot$core, core, tolerance = 1e-12)
  expect_equal(grid$support, core + c(-1, 1) * 0.1 * diff(core), tolerance = 1e-12)
  expect_equal(grid$pilot$share, diff(grid$support) / diff(wide), tolerance = 1e-12)
  expect_length(grid$grid, 101L)
  expect_equal(grid$spacing, diff(grid$support) / 100, tolerance = 1e-12)
  # The basis of rule "grid" does not depend on where the grid lies.
  expect_equal(unclass(grid$B), unclass(make_efron_grid(d$y, d$s, grid_method = "paper_realdata")$B),
               tolerance = 1e-10, ignore_attr = TRUE)

  # A larger margin gives a wider range; no margin still adds one spacing of
  # the pilot grid.
  wider <- make_efron_grid(d$y, d$s, grid_method = "pilot", seed = 1L, margin = 0.3)
  expect_equal(wider$support, core + c(-1, 1) * 0.3 * diff(core), tolerance = 1e-12)
  none <- make_efron_grid(d$y, d$s, grid_method = "pilot", seed = 1L, margin = 0)
  expect_equal(none$support, core + c(-1, 1) * diff(wide) / 100, tolerance = 1e-12)

  # `L`, `M` and the component apply to the grid that is returned.
  fine <- make_efron_grid(d$y, d$s, grid_method = "pilot", seed = 1L, L = 201L, M = 8L,
                          component = "normal")
  expect_equal(fine$support, grid$support, tolerance = 1e-12)
  expect_identical(c(fine$L, fine$M), c(201L, 8L))
  expect_equal(fine$tau_c, 2 / 3 * diff(grid$support) / 200, tolerance = 1e-12)
})

test_that("the rule never makes the range wider than that of \"paper_realdata\"", {
  d <- pilot_data()
  wide <- make_efron_grid(d$y, d$s, grid_method = "paper_realdata")$support
  # g against the upper end: that end stays, the lower one moves in.
  local_pilot(pilot_bell(wide[2L], 0.4))
  grid <- make_efron_grid(d$y, d$s, grid_method = "pilot", seed = 1L)
  expect_identical(grid$support[2L], wide[2L])
  expect_gt(grid$support[1L], wide[1L])
  expect_true(grid$pilot$narrowed)

  for (margin in c(0, 0.1, 0.5, 1)) {
    for (sd in c(0.01, 0.3, 2, 50)) {
      local_pilot(pilot_bell(-0.5, sd))
      one <- make_efron_grid(d$y, d$s, grid_method = "pilot", seed = 1L, margin = margin)
      expect_gte(one$support[1L], wide[1L])
      expect_lte(one$support[2L], wide[2L])
      expect_gt(diff(one$support), 0)
    }
  }
})

test_that("a pilot estimate on one grid point gives a range of two pilot spacings", {
  d <- pilot_data()
  wide <- make_efron_grid(d$y, d$s, grid_method = "paper_realdata")$support
  spacing <- diff(wide) / 100
  local_pilot(function(grid) as.numeric(seq_along(grid) == 41L))
  # The estimates are imprecise next to that range, so none moves an end;
  # the range is then narrower than a standard error, which is reported.
  expect_warning(
    grid <- make_efron_grid(d$y, 10 * d$s, grid_method = "pilot", seed = 1L),
    class = "bef_warn_grid_resolution"
  )
  point <- wide[1L] + 40 * spacing
  expect_equal(grid$pilot$core, c(point, point), tolerance = 1e-12)
  expect_equal(grid$support, point + c(-1, 1) * spacing, tolerance = 1e-12)
})

test_that("an estimate far outside the narrowed range in its own standard errors moves the end out", {
  support <- pilot_ns(".bef_pilot_support")
  grid <- seq(-10, 10, length.out = 101L)
  g <- stats::dnorm(grid, 0, 1)
  base <- support(grid, g, theta_hat = c(-1, 0, 1), sigma = c(1, 1, 1), margin = 0.1)
  expect_identical(base$protected, integer())
  upper <- base$support[2L]
  lower <- base$support[1L]

  # 3.9 standard errors outside: left alone. 4.1: the end goes to the
  # estimate plus two standard errors.
  near <- support(grid, g, theta_hat = c(-1, 0, upper + 3.9 * 0.2), sigma = c(1, 1, 0.2),
                  margin = 0.1)
  expect_identical(near$support, base$support)
  expect_identical(near$protected, integer())
  far <- support(grid, g, theta_hat = c(-1, 0, upper + 4.1 * 0.2), sigma = c(1, 1, 0.2),
                 margin = 0.1)
  expect_equal(far$support, c(lower, upper + 4.1 * 0.2 + 2 * 0.2), tolerance = 1e-12)
  expect_identical(far$protected, 3L)

  # Both ends, and the farthest site of each side decides.
  both <- support(grid, g, theta_hat = c(lower - 1, lower - 2, 0, upper + 1.5),
                  sigma = c(0.1, 0.2, 1, 0.1), margin = 0.1)
  expect_equal(both$support, c(lower - 2 - 0.4, upper + 1.5 + 0.2), tolerance = 1e-12)
  expect_identical(both$protected, c(1L, 2L, 4L))

  # The end still does not pass the wide range.
  capped <- support(grid, g, theta_hat = c(-1, 0, 9.9), sigma = c(1, 1, 0.5), margin = 0.1)
  expect_identical(capped$support[2L], 10)
  expect_identical(capped$protected, 3L)

  # An imprecise site far away in its own units of the data, but within four
  # of its standard errors, moves nothing.
  loose <- support(grid, g, theta_hat = c(-1, 0, 9), sigma = c(1, 1, 3), margin = 0.1)
  expect_identical(loose$support, base$support)
})

test_that("a precise outlier ends up inside the range with two standard errors to spare", {
  d <- pilot_data()
  y <- c(d$y, max(d$y) + 5 * stats::mad(d$y))
  s <- c(d$s, min(d$s) / 4)
  # The worst case: a pilot estimate with nothing near the outlier.
  local_pilot(pilot_bell(0, 0.6))
  grid <- make_efron_grid(y, s, grid_method = "pilot", seed = 1L)
  added <- length(y)
  expect_identical(grid$pilot$protected, added)
  expect_lte(grid$support[1L], y[added] - 2 * s[added])
  expect_gte(grid$support[2L], y[added] + 2 * s[added] - 1e-12)
  expect_true(grid$pilot$narrowed)
})

test_that("an imprecise outlier never widens the range beyond that of \"paper_realdata\"", {
  d <- pilot_data()
  y <- c(d$y, max(d$y) + 5 * stats::mad(d$y))
  s <- c(d$s, 4 * max(d$s))
  local_pilot(pilot_bell(0, 0.6))
  plain <- make_efron_grid(d$y, d$s, grid_method = "pilot", seed = 1L)
  with_outlier <- make_efron_grid(y, s, grid_method = "pilot", seed = 1L)
  old <- make_efron_grid(y, s, grid_method = "paper_realdata")
  expect_gte(with_outlier$support[1L], old$support[1L])
  expect_lte(with_outlier$support[2L], old$support[2L])
  # The outlier is within four of its standard errors of the range, so the
  # range is the one that g gives.
  expect_identical(with_outlier$pilot$protected, integer())
  expect_lt(diff(with_outlier$support), 0.75 * diff(old$support))
  expect_gt(diff(plain$support), 0)
})

test_that("shifting and rescaling the data shifts and rescales the grid", {
  d <- pilot_data()
  for (shape in list(pilot_flat, function(grid) {
    # The same shape whatever the units: a bell at 60% of the range.
    at <- (grid - grid[1L]) / (grid[length(grid)] - grid[1L])
    stats::dnorm(at, 0.6, 0.08)
  })) {
    local_pilot(shape)
    first <- make_efron_grid(d$y, d$s, grid_method = "pilot", seed = 1L)
    for (move in list(c(0, 1), c(1000, 1), c(0, 250), c(-3.5, 0.004))) {
      a <- move[[1L]]
      b <- move[[2L]]
      second <- make_efron_grid(a + b * d$y, b * d$s, grid_method = "pilot", seed = 1L)
      expect_equal(second$grid, a + b * first$grid, tolerance = 1e-9)
      expect_equal(second$pilot$share, first$pilot$share, tolerance = 1e-9)
      expect_identical(second$pilot$narrowed, first$pilot$narrowed)
    }
  }
})

test_that("the arguments of the rule are checked, and other rules ignore them", {
  d <- pilot_data()
  record <- local_pilot(pilot_flat)
  for (bad in list(-0.1, 1.5, NA_real_, "0.1", c(0.1, 0.2))) {
    err <- tryCatch(make_efron_grid(d$y, d$s, grid_method = "pilot", margin = bad),
                    error = identity)
    expect_s3_class(err, "bef_grid_error")
    expect_match(conditionMessage(err), "margin")
  }
  for (bad in list(-1, 1.5, NA_real_, "1", c(1, 2), 2^31)) {
    err <- tryCatch(make_efron_grid(d$y, d$s, grid_method = "pilot", seed = bad),
                    error = identity)
    expect_s3_class(err, "bef_grid_error")
    expect_match(conditionMessage(err), "seed")
  }
  expect_length(record$calls, 0L)

  # Without a seed one is taken from the clock and recorded.
  grid <- make_efron_grid(d$y, d$s, grid_method = "pilot")
  expect_type(grid$pilot$seed, "integer")
  expect_identical(record$calls[[1L]]$seed, grid$pilot$seed)

  plain <- make_efron_grid(d$y, d$s, margin = 0.3, seed = 5L, grid_method = "paper_realdata")
  expect_identical(plain$grid, make_efron_grid(d$y, d$s, grid_method = "paper_realdata")$grid)
  expect_length(record$calls, 1L)

  # `support` and a rule cannot both be given.
  expect_error(make_efron_grid(d$y, d$s, grid_method = "pilot", support = c(-3, 3)),
               "not both")
})

test_that("print() says what the rule did", {
  d <- pilot_data()
  local_pilot(pilot_flat)
  same <- format(make_efron_grid(d$y, d$s, grid_method = "pilot", seed = 1L))
  expect_match(same[[2L]], "\"pilot\", expansion 0.5, margin 0.1: not narrowed", fixed = TRUE)
  local_pilot(pilot_bell(0.2, 0.5))
  narrowed <- format(make_efron_grid(d$y, d$s, grid_method = "pilot", seed = 1L, margin = 0.2))
  expect_match(narrowed[[2L]], "margin 0.2: narrowed to [0-9.]+% of the range of \"paper_realdata\"")
})

test_that("equal estimates and a failed pilot fit give errors that say what to do", {
  local_pilot(pilot_flat)
  for (rule in c("pilot", "paper_realdata")) {
    err <- tryCatch(make_efron_grid(rep(0.3, 6), rep(0.1, 6), grid_method = rule),
                    error = identity)
    expect_s3_class(err, "bef_grid_error")
    expect_match(conditionMessage(err), "all equal")
    expect_match(conditionMessage(err), "`support`", fixed = TRUE)
  }

  d <- pilot_data()
  local_pilot(function(grid) rep(NA_real_, length(grid)))
  err <- tryCatch(make_efron_grid(d$y, d$s, grid_method = "pilot"), error = identity)
  expect_s3_class(err, "bef_pilot_failed")
  expect_s3_class(err, "bef_grid_error")
})

test_that("without CmdStan the rule stops and names the rules that need none", {
  d <- pilot_data()
  fit <- pilot_ns(".bef_pilot_fit")
  err <- tryCatch(
    fit(d$y, d$s, grid = seq(-4, 4, length.out = 101L), seed = 1L,
        check_installed_fun = function() stop("cmdstanr is not installed")),
    error = identity
  )
  expect_s3_class(err, "bef_pilot_failed")
  expect_match(conditionMessage(err), "needs CmdStan")
  expect_match(conditionMessage(err), "paper_realdata")
  expect_match(conditionMessage(err), "`support`", fixed = TRUE)

  # A sampler that fails is reported as a failure of the pilot fit.
  err <- tryCatch(
    fit(d$y, d$s, grid = seq(-4, 4, length.out = 101L), seed = 1L,
        check_installed_fun = function() invisible(TRUE),
        model_fun = function(...) list(sample = function(...) stop("no output"))),
    error = identity
  )
  expect_s3_class(err, "bef_pilot_failed")
  expect_match(conditionMessage(err), "pilot fit")
})

test_that("a fit with the rule runs the pilot fit with its own seed and chains", {
  d <- pilot_data()
  record <- local_pilot(pilot_bell(0.2, 0.5))
  prepare <- function(...) {
    pilot_ns(".bef_fit_prepare_context")(
      call = quote(bayes_efron_fit(theta_hat, sigma)), theta_hat = d$y, sigma = d$s,
      grid_method = "pilot", check_installed = FALSE, ...
    )
  }
  context <- prepare(seed = 31L, chains = 4L, parallel_chains = 2L)
  expect_length(record$calls, 1L)
  expect_identical(record$calls[[1L]]$seed, 31L)
  expect_identical(record$calls[[1L]]$parallel_chains, 2L)
  expect_identical(context$effective_seed, 31L)
  expect_identical(context$grid$grid_method, "pilot")
  expect_true(context$grid$pilot$narrowed)
  expect_identical(context$stan_data$grid, context$grid$grid)
  expect_identical(
    context$grid$grid,
    make_efron_grid(d$y, d$s, grid_method = "pilot", seed = 31L)$grid
  )

  # Without a seed, the fit and its pilot fit share the one taken from the
  # clock.
  context <- prepare(now = function() as.POSIXct("2026-10-09 12:00:00", tz = "UTC"))
  expect_identical(record$calls[[3L]]$seed, context$effective_seed)
  expect_identical(context$grid$pilot$seed, context$effective_seed)

  # Other rules run no pilot fit.
  pilot_ns(".bef_fit_prepare_context")(
    call = quote(bayes_efron_fit(theta_hat, sigma)), theta_hat = d$y, sigma = d$s,
    check_installed = FALSE
  )
  expect_length(record$calls, 3L)
})

test_that("the pilot fit runs with the settings of the rule", {
  d <- pilot_data()
  received <- NULL
  draws <- matrix(log(1 / 101), nrow = 2000L, ncol = 101L)
  fake <- list(sample = function(...) {
    received <<- list(...)
    list(
      draws = function(variables, format) draws,
      diagnostic_summary = function(quiet) list(num_divergent = c(0, 1, 0, 0),
                                                num_max_treedepth = c(2, 0, 0, 3)),
      summary = function(...) data.frame(rhat = c(1.01, 1.03), ess_bulk = c(900, 250)),
      output_files = function() character(),
      data_file = function() character()
    )
  })
  grid <- seq(-4, 4, length.out = 101L)
  out <- pilot_ns(".bef_pilot_fit")(
    d$y, d$s, grid = grid, seed = 11L, parallel_chains = 2L,
    check_installed_fun = function() invisible(TRUE),
    model_fun = function(program) {
      expect_identical(program, "efron_mix")
      fake
    }
  )
  expect_equal(out$g, rep(1 / 101, 101L), tolerance = 1e-12)
  expect_identical(out$diagnostics,
                   list(max_rhat = 1.03, min_ess_bulk = 250, divergences = 1, max_treedepth = 5))
  expect_identical(
    received[c("chains", "parallel_chains", "iter_warmup", "iter_sampling", "adapt_delta",
               "max_treedepth", "seed", "refresh", "init", "sig_figs")],
    list(chains = 4L, parallel_chains = 2L, iter_warmup = 500L, iter_sampling = 500L,
         adapt_delta = 0.9, max_treedepth = 10L, seed = 11L, refresh = 0L, init = 0.5,
         sig_figs = 12L)
  )
  expect_false(received$show_messages)
  expect_null(received$diagnostics)
  data <- received$data
  expect_named(data, c("K", "theta_hat", "sigma", "L", "grid", "M", "B", "tau_c",
                       "store_grid_quantities"))
  expect_identical(c(data$K, data$L, data$M), c(10L, 101L, 6L))
  expect_identical(data$tau_c, 0)
  expect_identical(data$grid, grid)
  expect_identical(dim(data$B), c(101L, 6L))
  expect_null(attr(data$B, "knots"))
  expect_equal(data$B, unclass(splines::ns(grid, df = 6)), ignore_attr = TRUE)
})

test_that("on the two small example data sets the rule leaves the range as it is", {
  skip_if_no_cmdstan()
  for (case in list(list(krw_discrimination, 20250807L), list(raudenbush1985, 1985L))) {
    grid <- make_efron_grid(case[[1L]], grid_method = "pilot", seed = case[[2L]])
    old <- make_efron_grid(case[[1L]], grid_method = "paper_realdata")
    expect_false(grid$pilot$narrowed)
    expect_identical(grid$grid, old$grid)
    expect_identical(grid$B, old$B)
    expect_true(is.finite(grid$pilot$diagnostics$max_rhat))
    expect_gt(grid$pilot$seconds, 0)
  }
  # The same seed gives the same pilot fit, and so the same grid.
  first <- make_efron_grid(raudenbush1985, grid_method = "pilot", seed = 1985L)
  expect_identical(first$pilot$core, grid$pilot$core)
  expect_identical(first$pilot$diagnostics, grid$pilot$diagnostics)
})

test_that("an outlier added to the callback data: an imprecise one widens the range no more than \"paper_realdata\" does, a precise one is covered", {
  skip_if_no_cmdstan()
  y <- krw_discrimination$yi
  s <- krw_discrimination$sei
  at <- max(y) + 5 * stats::mad(y)
  old <- make_efron_grid(y, s, grid_method = "paper_realdata")

  # An imprecise outlier, with ten times the median standard error: the
  # range grows by no more than that of 0.3.0.
  wide_se <- 10 * stats::median(s)
  loose <- make_efron_grid(c(y, at), c(s, wide_se), grid_method = "pilot", seed = 20250807L)
  loose_old <- make_efron_grid(c(y, at), c(s, wide_se), grid_method = "paper_realdata")
  expect_gte(loose$support[1L], loose_old$support[1L])
  expect_lte(loose$support[2L], loose_old$support[2L])
  expect_lte(diff(loose$support) / diff(old$support), diff(loose_old$support) / diff(old$support))

  # A precise outlier: its estimate and two standard errors on each side lie
  # inside.
  tight <- make_efron_grid(c(y, at), c(s, min(s)), grid_method = "pilot", seed = 20250807L)
  expect_lt(tight$support[1L], at - 2 * min(s))
  expect_gt(tight$support[2L], at + 2 * min(s))
})

test_that("a fit with the rule uses the grid that make_efron_grid() gives for its seed", {
  skip_if_no_cmdstan()
  fit <- suppressWarnings(bayes_efron_fit(
    raudenbush1985, grid_method = "pilot", chains = 2L, iter_warmup = 150L,
    iter_sampling = 50L, seed = 1985L
  ))
  grid <- make_efron_grid(raudenbush1985, grid_method = "pilot", seed = 1985L)
  expect_identical(fit$metadata$grid_method, "pilot")
  expect_identical(fit$metadata$data_list$grid, grid$grid)
  expect_identical(fit$metadata$seed, 1985L)
})

test_that("on the 1,500 simulated sites the rule narrows the range to about 0.3 of the wide one", {
  skip_if_no_cmdstan()
  grid <- make_efron_grid(twin_towers$yi, twin_towers$sei, grid_method = "pilot",
                          seed = 20250627L)
  old <- make_efron_grid(twin_towers$yi, twin_towers$sei, grid_method = "paper_realdata")
  expect_true(grid$pilot$narrowed)
  # With this seed the range is that of `twin_towers_fit`, -3.003 to 3.674,
  # 0.30 of the wide range. Another platform may put an end one or two
  # points of the pilot grid away.
  pilot_spacing <- diff(old$support) / 100
  expect_lt(abs(grid$support[1L] - (-3.003)), 2.2 * pilot_spacing)
  expect_lt(abs(grid$support[2L] - 3.674), 2.2 * pilot_spacing)
  expect_gt(grid$pilot$share, 0.25)
  expect_lt(grid$pilot$share, 0.36)
  expect_identical(grid$pilot$protected, integer())
  # No estimate is more than four of its standard errors outside.
  outside <- pmax((grid$support[1L] - twin_towers$yi) / twin_towers$sei,
                  (twin_towers$yi - grid$support[2L]) / twin_towers$sei)
  expect_lt(max(outside), 4)
})
