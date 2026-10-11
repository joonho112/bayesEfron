# compare_efron_grids(): the specifications it fits, what it keeps of each
# fit, and what it does when a fit fails. The fits are replaced by a stored
# one except in the tests at the end, which run CmdStan.

comparison_ns <- function(name) getFromNamespace(name, "bayesEfron")

# bayes_efron_fit() replaced by a function that records its arguments and
# returns the stored fit of the 19 studies, or fails where `fail` says so.
local_comparison_fits <- function(record = new.env(parent = emptyenv()),
                                  fail = function(args) FALSE,
                                  env = parent.frame()) {
  record$calls <- list()
  testthat::local_mocked_bindings(
    bayes_efron_fit = function(...) {
      args <- list(...)
      record$calls[[length(record$calls) + 1L]] <- args
      if (isTRUE(fail(args))) {
        stop("the sampler stopped")
      }
      raudenbush_fit
    },
    .package = "bayesEfron",
    .env = env
  )
  record
}

# A list of labelled estimates, as as_bef_data() takes it.
metafor_free_escalc <- function(data) {
  list(theta_hat = stats::setNames(data$yi, rownames(data)), sigma = data$sei)
}

comparison_quiet <- function(expr) {
  suppressMessages(expr, classes = "bef_progress")
}

test_that("every combination of the values is fitted and reported, the best first", {
  skip_if_not_installed("loo")
  record <- local_comparison_fits()
  out <- comparison_quiet(compare_efron_grids(raudenbush1985, seed = 3L, chains = 2L))
  expect_s3_class(out, c("bef_grid_comparison", "data.frame"), exact = TRUE)
  expect_identical(nrow(out), 9L)
  expect_length(record$calls, 9L)
  expect_named(out, c("spec", "L", "M", "expansion", "grid_method", "component", "lower",
                      "upper", "elpd_loo", "se_elpd_loo", "p_loo", "elpd_diff", "se_diff",
                      "mean_g", "sd_g", "site_8", "site_4", "max_rhat", "min_ess_bulk",
                      "divergences", "max_pareto_k", "end_mass", "adequate",
                      "sampler_seconds"))
  expect_setequal(out$spec, 1:9)
  expect_setequal(paste(out$L, out$M), paste(rep(c(51, 101, 201), 3), rep(c(4, 6, 8), each = 3)))
  expect_identical(unique(out$grid_method), "paper_realdata")
  expect_identical(unique(out$component), "point")
  expect_identical(unique(out$expansion), 0.5)

  # What each fit was called with: the data as two vectors, the values of
  # the specification, and the other arguments as they were given.
  first <- record$calls[[1L]]
  expect_equal(unname(first$theta_hat), raudenbush1985$yi)
  expect_equal(first$sigma, raudenbush1985$sei)
  expect_identical(first[c("L", "M", "component", "grid_method", "expansion", "seed", "chains")],
                   list(L = 51L, M = 4L, component = "point", grid_method = "paper_realdata",
                        expansion = 0.5, seed = 3L, chains = 2L))

  # The same fit every time: the same numbers in every row.
  expect_equal(out$sd_g, rep(mean(raudenbush_fit$posterior$sd_g), 9L))
  expect_equal(out$site_4, rep(site_summary(raudenbush_fit, sites = 4L)$mean, 9L))
  expect_equal(out$elpd_diff, rep(0, 9L))
  expect_equal(out$end_mass, rep(max(diagnose(raudenbush_fit)$grid$ends$mass), 9L))
  expect_identical(out$lower, rep(min(raudenbush_fit$metadata$data_list$grid), 9L))

  # What is kept of each fit.
  g <- attr(out, "g", exact = TRUE)
  expect_named(g, as.character(1:9))
  expect_identical(g[["1"]], g_summary(raudenbush_fit, type = "density"))
  sites <- attr(out, "sites", exact = TRUE)
  expect_identical(nrow(sites), 18L)
  expect_identical(sites$spec, rep(1:9, each = 2L))
  expect_identical(sites$site, rep(c(8L, 4L), 9L))
  expect_null(attr(out, "fits", exact = TRUE))

  kept <- comparison_quiet(compare_efron_grids(raudenbush1985, L = 101L, M = 6L, keep_fits = TRUE))
  expect_named(attr(kept, "fits", exact = TRUE), "1")
  expect_s3_class(attr(kept, "fits", exact = TRUE)[["1"]], "bef_fit")
})

test_that("ranges, rules and components can be compared, and the sites can be chosen", {
  skip_if_not_installed("loo")
  record <- local_comparison_fits()
  out <- comparison_quiet(compare_efron_grids(
    raudenbush1985$yi, raudenbush1985$sei, L = 101L, M = 6L, expansion = c(0.25, 0.5, 1),
    sites = c(4L, 10L)
  ))
  expect_identical(nrow(out), 3L)
  expect_setequal(out$expansion, c(0.25, 0.5, 1))
  expect_true(all(c("site_4", "site_10") %in% names(out)))
  expect_identical(attr(out, "sites", exact = TRUE)$site, rep(c(4L, 10L), 3L))

  # Ranges given directly take the place of the rule and of `expansion`.
  out <- comparison_quiet(compare_efron_grids(
    raudenbush1985, L = 101L, M = c(6L, 8L), support = list(c(-1, 2), c(-0.5, 1.5))
  ))
  expect_identical(nrow(out), 4L)
  expect_identical(unique(out$grid_method), "support")
  expect_true(all(is.na(out$expansion)))
  calls <- record$calls[4:7]
  expect_identical(lapply(calls, `[[`, "support"),
                   list(c(-1, 2), c(-1, 2), c(-0.5, 1.5), c(-0.5, 1.5)))
  expect_false(any(vapply(calls, function(args) "grid_method" %in% names(args), logical(1))))
  one <- comparison_quiet(compare_efron_grids(raudenbush1985, L = 101L, M = 6L,
                                              support = c(-1, 2)))
  expect_identical(nrow(one), 1L)

  # Rules and components.
  out <- comparison_quiet(compare_efron_grids(
    raudenbush1985, L = 101L, M = 6L, grid_method = c("paper_realdata", "pilot"),
    component = c("point", "normal")
  ))
  expect_identical(nrow(out), 4L)
  expect_setequal(paste(out$grid_method, out$component),
                  c("paper_realdata point", "pilot point", "paper_realdata normal",
                    "pilot normal"))

  # Sites by label.
  labelled <- raudenbush1985
  rownames(labelled) <- make.unique(labelled$author)
  dat <- as_bef_data(metafor_free_escalc(labelled))
  out <- comparison_quiet(compare_efron_grids(dat, L = 101L, M = 6L, sites = dat$names[[4L]]))
  expect_true(paste0("site_", dat$names[[4L]]) %in% names(out))
})

test_that("what every fit would reject is rejected once, before any fit", {
  skip_if_not_installed("loo")
  record <- local_comparison_fits()
  # Three sites: a fit needs at least five.
  err <- tryCatch(compare_efron_grids(c(0.1, 0.2, 0.3), c(0.1, 0.1, 0.1)), error = identity)
  expect_s3_class(err, "bef_invalid_args")
  expect_match(conditionMessage(err), "at least 5 values")
  # A grid cannot be given, and its name is not taken for `grid_method`.
  grid <- make_efron_grid(raudenbush1985, grid_method = "paper_realdata")
  err <- tryCatch(compare_efron_grids(raudenbush1985, grid = grid), error = identity)
  expect_s3_class(err, "bef_invalid_args")
  expect_match(conditionMessage(err), "builds the grids that it compares")
  cases <- list(
    quote(compare_efron_grids(raudenbush1985, L = 20L)),
    quote(compare_efron_grids(raudenbush1985, M = 12L)),
    quote(compare_efron_grids(raudenbush1985, expansion = -1)),
    quote(compare_efron_grids(raudenbush1985, grid_method = "widest")),
    quote(compare_efron_grids(raudenbush1985, grid_method = "kl_target_experimental")),
    quote(compare_efron_grids(raudenbush1985, component = "bump")),
    quote(compare_efron_grids(raudenbush1985, support = c(2, 1))),
    quote(compare_efron_grids(raudenbush1985, sites = 25L)),
    quote(compare_efron_grids(raudenbush1985, sites = "no such study")),
    quote(compare_efron_grids(raudenbush1985, keep_fits = NA)),
    quote(compare_efron_grids(raudenbush1985, grid = make_efron_grid(raudenbush1985,
                                                                    support = c(-1, 2)))),
    quote(compare_efron_grids(raudenbush1985$yi, raudenbush1985$sei[-1]))
  )
  for (case in cases) {
    err <- tryCatch(suppressMessages(eval(case)), error = identity)
    expect_s3_class(err, "bef_error")
  }
  expect_length(record$calls, 0L)
})

test_that("a fit that fails gives a warning and a row of missing values, and the others are made", {
  skip_if_not_installed("loo")
  record <- local_comparison_fits(fail = function(args) args$L == 101L && args$M == 6L)
  caught <- list()
  out <- withCallingHandlers(
    comparison_quiet(compare_efron_grids(raudenbush1985)),
    bef_grid_comparison_fit_failed = function(w) {
      caught[[length(caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_length(caught, 1L)
  expect_s3_class(caught[[1L]], "bef_warning")
  expect_match(conditionMessage(caught[[1L]]),
               "The fit failed for L = 101, M = 6, rule \"paper_realdata\", expansion 0.5: the sampler stopped",
               fixed = TRUE)
  expect_length(record$calls, 9L)
  expect_identical(nrow(out), 9L)
  failed <- out[out$L == 101L & out$M == 6L, ]
  expect_identical(failed$spec, caught[[1L]]$spec)
  expect_true(all(is.na(failed[c("elpd_loo", "sd_g", "site_4", "max_rhat", "end_mass",
                                 "adequate", "lower", "elpd_diff")])))
  # The failed fit is the last row, and the others are complete.
  expect_identical(which(is.na(out$elpd_loo)), 9L)
  expect_false(anyNA(out$sd_g[1:8]))
  expect_null(attr(out, "g", exact = TRUE)[[as.character(failed$spec)]])
  expect_identical(nrow(attr(out, "sites", exact = TRUE)), 16L)

  # All fits fail: a table of missing values, and no error.
  local_comparison_fits(fail = function(args) TRUE)
  none <- suppressWarnings(comparison_quiet(compare_efron_grids(raudenbush1985, L = 101L)))
  expect_identical(nrow(none), 3L)
  expect_true(all(is.na(none$elpd_loo)))
  expect_output(print(none), "3 specifications")
})

test_that("a line says which fit is running, and print() reads the table", {
  skip_if_not_installed("loo")
  local_comparison_fits()
  lines <- character()
  out <- withCallingHandlers(
    compare_efron_grids(raudenbush1985, L = 101L, M = c(6L, 8L)),
    bef_progress = function(m) {
      lines <<- c(lines, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_identical(trimws(lines), c(
    "Fit 1 of 2: L = 101, M = 6, rule \"paper_realdata\", expansion 0.5",
    "Fit 2 of 2: L = 101, M = 8, rule \"paper_realdata\", expansion 0.5"
  ))
  shown <- utils::capture.output(print(out))
  expect_match(shown[[1L]], "<bef_grid_comparison> 2 specifications, best first", fixed = TRUE)
  expect_true(any(grepl("In every row: grid_method = paper_realdata, component = point, expansion = 0.5, L = 101",
                        shown, fixed = TRUE)))
  expect_true(any(grepl("elpd_loo does not separate these specifications", shown, fixed = TRUE)))
  expect_true(any(grepl("sd_g", shown, fixed = TRUE)))

  # A table made by version 0.3.0, which has fewer columns, still prints.
  old <- readRDS(test_path("_fixtures", "grid_comparison_v030.rds"))$raudenbush
  expect_output(print(old), "9 specifications, best first")
})

test_that("for the 19 studies the default call runs and the range moves the standard deviation of g", {
  skip_if_no_cmdstan()
  skip_if_not_installed("loo")
  short <- list(chains = 2L, iter_warmup = 300L, iter_sampling = 300L, seed = 1985L)
  out <- suppressWarnings(comparison_quiet(do.call(
    compare_efron_grids, c(list(raudenbush1985), short)
  )))
  expect_identical(nrow(out), 9L)
  expect_false(anyNA(out$elpd_loo))
  expect_true(all(is.finite(out$sd_g)))
  expect_identical(out$elpd_diff[[1L]], 0)
  # The coarsest grids are too coarse for the spread of the effects.
  expect_false(any(out$adequate[out$L == 51L]))
  expect_true(all(out$adequate[out$L == 201L]))

  ranges <- suppressWarnings(comparison_quiet(do.call(
    compare_efron_grids,
    c(list(raudenbush1985, L = 101L, M = 6L, expansion = c(0.25, 0.5, 1), sites = 4L), short)
  )))
  ranges <- ranges[order(ranges$expansion), ]
  expect_identical(ranges$expansion, c(0.25, 0.5, 1))
  expect_equal(ranges$upper - ranges$lower, 1.5 * (1 + 2 * c(0.25, 0.5, 1)))
  # The table shows what the range does to the study with the largest
  # estimate: at full length its posterior mean is 0.78, 0.66 and 0.44 on
  # the three ranges, and the standard deviation of g 0.32, 0.31 and 0.30.
  expect_gt(ranges$site_4[[1L]] - ranges$site_4[[3L]], 0.15)
  expect_true(all(ranges$sd_g > 0.2 & ranges$sd_g < 0.45))
  # The narrowest range is the one on which a study leans, and 101 points
  # are too few for the widest.
  expect_gt(ranges$end_mass[[1L]], ranges$end_mass[[3L]])
  expect_lt(ranges$end_mass[[3L]], 0.01)
  expect_identical(ranges$adequate, c(TRUE, TRUE, FALSE))
  g <- attr(ranges, "g", exact = TRUE)
  expect_length(g, 3L)
  expect_named(g[[1L]], c("value", "mean", "lower", "upper"))
})
