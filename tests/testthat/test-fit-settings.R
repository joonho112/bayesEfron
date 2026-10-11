# The arguments of bayes_efron_fit() that describe the grid, its defaults,
# and `settings = "0.3.0"`. Sampling is replaced by a stand-in that records
# what the fit was prepared with; the tests at the end run CmdStan.

settings_ns <- function(name) getFromNamespace(name, "bayesEfron")

settings_data <- function() {
  list(
    y = c(-1.9, -1.2, -0.6, -0.3, 0, 0.1, 0.4, 0.8, 1.3, 2.1),
    s = c(0.5, 0.4, 0.3, 0.2, 0.2, 0.2, 0.3, 0.3, 0.4, 0.5)
  )
}

# bayes_efron_fit() up to the point where sampling would start: the context
# that it prepared is returned in place of a fit. The pilot fit is replaced
# by an estimate of g that fills the range, or by `g_fun`.
settings_prepared <- function(..., g_fun = function(grid) rep(1, length(grid)),
                              record = new.env(parent = emptyenv())) {
  record$pilot <- 0L
  testthat::local_mocked_bindings(
    .bef_check_cmdstanr_installed = function() invisible(TRUE),
    .bef_pilot_fit = function(theta_hat, sigma, grid, seed, parallel_chains = NULL, ...) {
      record$pilot <- record$pilot + 1L
      record$pilot_seed <- seed
      g <- g_fun(grid)
      list(g = g / sum(g),
           diagnostics = list(max_rhat = 1, min_ess_bulk = 1000, divergences = 0L,
                              max_treedepth = 0L),
           seconds = 0)
    },
    .bef_fit_run_sampler = function(context, ...) context,
    .bef_fit_assemble = function(context, ...) context,
    .package = "bayesEfron"
  )
  suppressMessages(bayes_efron_fit(...), classes = "bef_default_changed")
}

settings_reset_messages <- function(env = parent.frame()) {
  emitted <- settings_ns(".bayesEfron_msgs_emitted")
  old <- emitted[["default_grid_0.4"]]
  emitted[["default_grid_0.4"]] <- NULL
  withr::defer(emitted[["default_grid_0.4"]] <- old, envir = env)
}

test_that("the default fit uses the rule \"pilot\", the program of this version and the draws of g", {
  d <- settings_data()
  record <- new.env(parent = emptyenv())
  context <- settings_prepared(d$y, d$s, seed = 5L, record = record)
  expect_identical(record$pilot, 1L)
  expect_identical(record$pilot_seed, 5L)
  expect_identical(context$grid$grid_method, "pilot")
  expect_identical(context$args$grid_method, "pilot")
  expect_identical(context$program, "efron_mix")
  expect_null(context$args$settings)
  expect_false(context$args$site_draws)
  expect_identical(context$grid$L, 101L)
  expect_identical(context$grid$M, 6L)
  expect_identical(context$grid$component, "point")
  expect_identical(context$stan_data$tau_c, 0)
  # A pilot estimate that fills the range leaves the grid of 0.3.0.
  expect_identical(context$stan_data$grid,
                   make_efron_grid(d$y, d$s, grid_method = "paper_realdata")$grid)

  formals <- formals(bayes_efron_fit)
  expect_identical(formals$grid_method, "pilot")
  expect_null(formals$L)
  expect_identical(formals$M, 6L)
  expect_identical(formals$max_treedepth, 10L)
  expect_false(formals$site_draws)
  expect_null(formals$settings)
  expect_null(formals$grid)
})

test_that("the arguments of the grid reach the grid of the fit", {
  d <- settings_data()
  record <- new.env(parent = emptyenv())
  context <- settings_prepared(d$y, d$s, support = c(-3, 3), L = 121L, M = 8L,
                               component = "normal", width = 0.5, record = record)
  expect_identical(record$pilot, 0L)
  expect_identical(context$grid$grid_method, "support")
  expect_identical(context$args$grid_method, "support")
  expect_identical(context$grid$support, c(-3, 3))
  expect_identical(c(context$grid$L, context$grid$M), c(121L, 8L))
  expect_identical(context$grid$component, "normal")
  expect_equal(context$stan_data$tau_c, 0.5 * 6 / 120)
  expect_identical(context$stan_data$L, 121L)
  expect_identical(dim(context$stan_data$B), c(121L, 8L))

  bell <- function(grid) stats::dnorm(grid, 0.2, 0.5)
  narrowed <- settings_prepared(d$y, d$s, margin = 0.3, seed = 2L, g_fun = bell)
  expect_true(narrowed$grid$pilot$narrowed)
  expect_identical(narrowed$grid$pilot$margin, 0.3)
  # The range that the rule gives for that estimate of g and that margin.
  wide <- make_efron_grid(d$y, d$s, grid_method = "paper_realdata")$support
  pilot_grid <- seq(wide[1L], wide[2L], length.out = 101L)
  found <- settings_ns(".bef_pilot_support")(pilot_grid, bell(pilot_grid), d$y, d$s,
                                             margin = 0.3, wide = wide)
  expect_identical(narrowed$grid$support, found$support)
  expect_identical(narrowed$stan_data$grid,
                   seq(found$support[1L], found$support[2L], length.out = 101L))

  expect_error(settings_prepared(d$y, d$s, support = c(-3, 3), grid_method = "pilot"),
               "not both", class = "bef_grid_error")
  expect_error(settings_prepared(d$y, d$s, knots = "estimates"), "`knots` must be",
               class = "bef_grid_error")
  expect_error(settings_prepared(d$y, d$s, component = "bump"), class = "bef_grid_error")
  expect_error(settings_prepared(d$y, d$s, margin = 2), "margin", class = "bef_grid_error")
})

test_that("a grid that was built beforehand is used as it is, and alone", {
  d <- settings_data()
  record <- new.env(parent = emptyenv())
  grid <- make_efron_grid(d$y, d$s, support = c(-4, 4), L = 151L, M = 7L)
  context <- settings_prepared(d$y, d$s, grid = grid, seed = 9L, record = record)
  expect_identical(record$pilot, 0L)
  expect_identical(context$stan_data$grid, grid$grid)
  expect_identical(context$stan_data$B, unname(matrix(as.numeric(grid$B), nrow = 151L)))
  expect_identical(c(context$args$L, context$args$M), c(151L, 7L))
  expect_identical(context$args$grid_method, "support")
  expect_identical(context$grid$knots, grid$knots)

  # A grid of the rule "pilot" keeps its record, and no second pilot fit is
  # run.
  made <- settings_ns(".bef_build_grid")
  testthat::local_mocked_bindings(
    .bef_pilot_fit = function(theta_hat, sigma, grid, seed, ...) {
      g <- stats::dnorm(grid, 0, 0.6)
      list(g = g / sum(g), diagnostics = list(max_rhat = 1), seconds = 0)
    },
    .package = "bayesEfron"
  )
  pilot_grid <- made(d$y, d$s, grid_method = "pilot", seed = 4L)
  context <- settings_prepared(d$y, d$s, grid = pilot_grid, record = record)
  expect_identical(record$pilot, 0L)
  expect_identical(context$args$grid_method, "pilot")
  expect_identical(context$grid$pilot, pilot_grid$pilot)

  for (also in list(list(L = 101L), list(M = 6L), list(grid_method = "pilot"),
                    list(support = c(-1, 1)), list(component = "normal"),
                    list(expansion = 1), list(margin = 0.2), list(knots = "grid"),
                    list(width = 0.5), list(theta_true = d$y), list(bound_expansion = 1))) {
    err <- tryCatch(
      do.call(settings_prepared, c(list(d$y, d$s, grid = grid), also)),
      error = identity
    )
    expect_s3_class(err, "bef_invalid_args")
    expect_identical(err$arg, "grid")
    expect_match(conditionMessage(err), sprintf("`%s` must be left out", names(also)), fixed = TRUE)
  }
  err <- tryCatch(settings_prepared(d$y, d$s, grid = grid, L = 101L, M = 6L), error = identity)
  expect_match(conditionMessage(err), "`L`, `M` must be left out", fixed = TRUE)

  for (bad in list(grid$grid, list(grid = grid$grid), "pilot", TRUE)) {
    err <- tryCatch(settings_prepared(d$y, d$s, grid = bad), error = identity)
    expect_s3_class(err, "bef_invalid_args")
    expect_match(conditionMessage(err), "made by make_efron_grid()", fixed = TRUE)
  }

  # The grid is compared with the data of the fit, which may be other data
  # than it was built from.
  expect_warning(
    settings_prepared(d$y, d$s / 20, grid = grid),
    class = "bef_warn_grid_resolution"
  )
})

test_that("settings = \"0.3.0\" replaces the defaults that were not given", {
  d <- settings_data()
  record <- new.env(parent = emptyenv())
  context <- settings_prepared(d$y, d$s, settings = "0.3.0", record = record)
  expect_identical(record$pilot, 0L)
  expect_identical(context$program, "efron_re")
  expect_identical(context$args$settings, "0.3.0")
  expect_identical(context$args$grid_method, "paper_realdata")
  expect_identical(context$args$L, 101L)
  expect_true(context$args$site_draws)
  expect_false("tau_c" %in% names(context$stan_data))
  expect_identical(context$stan_data$grid,
                   make_efron_grid(d$y, d$s, grid_method = "paper_realdata", L = 101L)$grid)

  # With standard errors that ask for a finer grid, version 0.3.0 still has
  # 101 points; it says that they are too few.
  expect_warning(
    coarse <- settings_prepared(d$y, d$s / 20, settings = "0.3.0"),
    class = "bef_warn_grid_resolution"
  )
  expect_identical(coarse$args$L, 101L)

  # An argument that is given keeps its value.
  context <- settings_prepared(d$y, d$s, settings = "0.3.0", site_draws = FALSE, L = 151L,
                               M = 8L, expansion = 1)
  expect_false(context$args$site_draws)
  expect_identical(c(context$args$L, context$args$M), c(151L, 8L))
  expect_identical(context$program, "efron_re")
  expect_identical(context$args$grid_method, "paper_realdata")
  context <- settings_prepared(d$y, d$s, settings = "0.3.0", grid_method = "pilot",
                               record = record)
  expect_identical(context$args$grid_method, "pilot")
  expect_identical(record$pilot, 1L)
  context <- settings_prepared(d$y, d$s, settings = "0.3.0", support = c(-3, 3))
  expect_identical(context$args$grid_method, "support")
  expect_identical(context$args$L, 101L)

  # The program of that version has point masses only.
  err <- tryCatch(settings_prepared(d$y, d$s, settings = "0.3.0", component = "normal"),
                  error = identity)
  expect_s3_class(err, "bef_error")
  expect_match(conditionMessage(err), "point masses only")

  for (bad in list("0.4.0", "old", 0.3, TRUE)) {
    err <- tryCatch(settings_prepared(d$y, d$s, settings = bad), error = identity)
    expect_s3_class(err, "bef_invalid_args")
    expect_identical(err$arg, "settings")
  }
})

test_that("the change of the default grid is announced once in a session", {
  d <- settings_data()
  settings_reset_messages()
  testthat::local_mocked_bindings(
    .bef_check_cmdstanr_installed = function() invisible(TRUE),
    .bef_pilot_fit = function(theta_hat, sigma, grid, seed, ...) {
      list(g = rep(1 / length(grid), length(grid)), diagnostics = list(), seconds = 0)
    },
    .bef_fit_run_sampler = function(context, ...) context,
    .bef_fit_assemble = function(context, ...) context,
    .package = "bayesEfron"
  )
  emitted <- settings_ns(".bayesEfron_msgs_emitted")

  # No message for a call that names its rule, its range, its grid or the
  # settings of 0.3.0.
  expect_no_message(bayes_efron_fit(d$y, d$s, grid_method = "pilot"))
  expect_no_message(bayes_efron_fit(d$y, d$s, grid_method = "paper_realdata"))
  expect_no_message(bayes_efron_fit(d$y, d$s, support = c(-3, 3)))
  expect_no_message(bayes_efron_fit(d$y, d$s, settings = "0.3.0"))
  expect_no_message(bayes_efron_fit(d$y, d$s, grid = make_efron_grid(d$y, d$s, grid_method = "paper_realdata")))
  expect_null(emitted[["default_grid_0.4"]])

  message <- tryCatch(bayes_efron_fit(d$y, d$s), message = identity)
  expect_s3_class(message, "bef_default_changed")
  expect_match(conditionMessage(message), "settings = \"0.3.0\"", fixed = TRUE)
  expect_match(conditionMessage(message), "once in a session")
  expect_no_message(bayes_efron_fit(d$y, d$s))
  expect_no_message(bayes_efron_fit(d$y, d$s, seed = 3L))
})

test_that("a fit records its grid, its program and its settings", {
  skip_if_no_cmdstan()
  fit <- suppressWarnings(suppressMessages(bayes_efron_fit(
    raudenbush1985, chains = 2L, iter_warmup = 150L, iter_sampling = 50L, seed = 1985L
  )))
  grid <- attr(fit$metadata, "grid", exact = TRUE)
  expect_s3_class(grid, "bef_grid")
  expect_identical(grid$grid_method, "pilot")
  expect_identical(fit$metadata$grid_method, "pilot")
  expect_false(grid$pilot$narrowed)
  expect_identical(grid$pilot$seed, 1985L)
  expect_identical(grid$grid, fit$metadata$data_list$grid)
  expect_identical(attr(fit$metadata, "program", exact = TRUE), "efron_mix")
  expect_null(attr(fit$metadata, "settings", exact = TRUE))
  expect_identical(
    fit$metadata$stan_file_sha256,
    digest::digest(settings_ns(".bef_stan_file")(program = "efron_mix"), algo = "sha256",
                   file = TRUE)
  )
  expect_identical(fit$metadata$cmdstan_version, as.character(cmdstanr::cmdstan_version()))
  expect_false(attr(fit$metadata, "storage", exact = TRUE)$site_draws)
  shown <- format(fit)
  expect_true(any(grepl("rule \"pilot\"", shown, fixed = TRUE)))

  # The grid of one fit can be given to another.
  again <- suppressWarnings(bayes_efron_fit(
    raudenbush1985, grid = grid, chains = 2L, iter_warmup = 150L, iter_sampling = 50L,
    seed = 1985L
  ))
  expect_identical(again$metadata$data_list$grid, fit$metadata$data_list$grid)
  expect_identical(unclass(again$draws), unclass(fit$draws))

  # Normal components and a range that was given.
  normal <- suppressWarnings(bayes_efron_fit(
    raudenbush1985, support = c(-1, 1.5), component = "normal", chains = 2L,
    iter_warmup = 150L, iter_sampling = 50L, seed = 1985L
  ))
  expect_identical(normal$metadata$grid_method, "support")
  expect_gt(normal$metadata$data_list$tau_c, 0)
  expect_equal(mean(normal$posterior$var_g) >= normal$metadata$data_list$tau_c^2, TRUE)
  expect_true(all(is.finite(normal$metadata$theta_summary$mean)))
  expect_true(any(grepl("ends given in `support`", format(normal), fixed = TRUE)))
})

test_that("settings = \"0.3.0\" gives the draws of version 0.3.0", {
  skip_if_no_cmdstan()
  ref <- readRDS(test_path("_fixtures", "ref-v030.rds"))
  # The same draws need the same CmdStan on the same kind of machine as the
  # reference. Elsewhere only the summaries can be compared.
  same_platform <- identical(as.character(cmdstanr::cmdstan_version()),
                             ref$conditions$cmdstan_version) &&
    identical(R.version$arch, ref$conditions$arch) &&
    identical(R.version$platform, ref$conditions$platform)
  draw_matrix <- function(fit) {
    m <- posterior::as_draws_matrix(fit$draws)
    matrix(as.numeric(m), nrow = nrow(m), dimnames = list(NULL, colnames(m)))
  }
  for (name in names(ref$cases)) {
    case <- ref$cases[[name]]
    # The call of the reference, with `settings` in place of the rule and
    # the number of points where they are the defaults of 0.3.0.
    args <- case$args
    if (identical(args$L, 101L)) {
      args$L <- NULL
    }
    args$grid_method <- NULL
    fit <- suppressWarnings(do.call(bayes_efron_fit, c(args, list(settings = "0.3.0"))))
    expect_identical(attr(fit$metadata, "program", exact = TRUE), "efron_re")
    expect_identical(attr(fit$metadata, "settings", exact = TRUE), "0.3.0")
    expect_identical(fit$metadata$stan_file_sha256, ref$conditions$stan_file_sha256)
    expect_identical(fit$metadata$grid_method, "paper_realdata")
    expect_identical(attr(fit$metadata, "sampler_settings", exact = TRUE)[names(case$sampler_settings)],
                     case$sampler_settings)
    draws <- draw_matrix(fit)
    expect_identical(colnames(draws), case$variables)
    expect_equal(range(fit$metadata$data_list$grid), case$grid_range)
    if (same_platform) {
      expect_identical(draws[, colnames(case$draws)], case$draws)
      blocks <- sub("\\[.*$", "", colnames(draws))
      hashes <- vapply(unique(blocks), function(b) {
        digest::digest(unname(draws[, blocks == b, drop = FALSE]), algo = "sha256")
      }, character(1L))
      expect_identical(hashes, case$block_sha256)
    } else {
      for (variable in c("mean_g", "sd_g")) {
        expect_equal(mean(draws[, variable]), mean(case$draws[, variable]),
                     tolerance = 0.2)
      }
    }
  }
  if (!same_platform) {
    skip(sprintf(
      "The draws were compared in summary only: the reference was made with CmdStan %s on %s.",
      ref$conditions$cmdstan_version, ref$conditions$platform
    ))
  }
})

test_that("thin keeps one iteration in every thin, with the diagnostics of the whole run", {
  skip_if_no_cmdstan()
  run <- function(...) {
    suppressWarnings(bayes_efron_fit(
      raudenbush1985, grid_method = "paper_realdata", chains = 2L, iter_warmup = 200L,
      iter_sampling = 200L, seed = 1985L, ...
    ))
  }
  full <- run()
  thinned <- run(thin = 4L)
  expect_identical(dim(thinned$draws)[1:2], c(50L, 2L))
  expect_identical(dimnames(thinned$draws)[[3L]], dimnames(full$draws)[[3L]])
  kept <- seq(1L, 200L, by = 4L)
  expect_equal(unclass(thinned$draws)[, , "lambda"], unclass(full$draws)[kept, , "lambda"],
               ignore_attr = TRUE)
  expect_equal(unclass(thinned$draws)[, , "log_g[40]"], unclass(full$draws)[kept, , "log_g[40]"],
               ignore_attr = TRUE)
  expect_equal(unclass(thinned$draws)[, , "effective_params"],
               unclass(full$draws)[kept, , "effective_params"], ignore_attr = TRUE)
  expect_length(thinned$posterior$sd_g, 100L)
  storage <- attr(thinned$metadata, "storage", exact = TRUE)
  expect_identical(c(storage$iterations, storage$chains, storage$thin), c(50L, 2L, 4L))
  expect_identical(attr(full$metadata, "storage", exact = TRUE)$thin, 1L)
  expect_true(any(grepl("Draws: 100 (2 chains, one iteration in 4 kept)", format(thinned),
                        fixed = TRUE)))
  expect_true(any(grepl("Draws: 400 (2 chains)", format(full), fixed = TRUE)))
  # The sampler quantities of each iteration go with the draws that are kept;
  # the counts are those of the whole run.
  expect_identical(dim(attr(thinned$metadata, "nuts", exact = TRUE))[1:2], c(50L, 2L))
  d_full <- attr(full$metadata, "diagnostics", exact = TRUE)
  d_thin <- attr(thinned$metadata, "diagnostics", exact = TRUE)
  expect_identical(d_thin$divergences, d_full$divergences)
  expect_identical(d_thin$max_treedepth, d_full$max_treedepth)
  expect_equal(d_thin$ebfmi, d_full$ebfmi)
  # So are R-hat and the effective sample sizes: thinning reduces what is
  # stored and nothing else.
  expect_identical(d_thin$rhat, d_full$rhat)
  expect_identical(d_thin$ess_bulk, d_full$ess_bulk)
  expect_identical(d_thin$ess_tail, d_full$ess_tail)
  expect_identical(attr(thinned$metadata, "diagnostic_skipped", exact = TRUE),
                   attr(full$metadata, "diagnostic_skipped", exact = TRUE))
  # The table of the sites is that of the kept draws.
  expect_equal(site_summary(thinned)$mean, thinned$metadata$theta_summary$mean)
  expect_true(all(is.finite(site_summary(thinned)$lower)))

  # With site-level draws the blocks are computed for the kept draws only.
  both <- run(thin = 4L, site_draws = TRUE)
  expect_identical(dim(both$draws)[1:2], c(50L, 2L))
  shared <- dimnames(thinned$draws)[[3L]]
  expect_identical(unclass(both$draws)[, , shared], unclass(thinned$draws)[, , shared])

  # The program of version 0.3.0 is thinned in the same way.
  old_full <- run(settings = "0.3.0")
  old_thin <- run(settings = "0.3.0", thin = 4L)
  expect_identical(dim(old_thin$draws)[1:2], c(50L, 2L))
  expect_equal(unclass(old_thin$draws)[, , "theta_rep[4]"],
               unclass(old_full$draws)[kept, , "theta_rep[4]"], ignore_attr = TRUE)
})
