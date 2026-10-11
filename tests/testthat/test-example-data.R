test_that("the callback gaps and the simulated sites keep their sample, their units and the record of their source", {
  expect_equal(nrow(krw_discrimination), 97L)
  expect_equal(sum(krw_discrimination$n_obs), 78910)
  expect_equal(sum(krw_discrimination$n_jobs), 10453)
  expect_identical(krw_discrimination$site, seq_len(97L))
  expect_equal(length(unique(krw_discrimination$firm_id)), 97L)
  expect_true(all(diff(krw_discrimination$firm_id) > 0))
  expect_equal(nrow(twin_towers), 1500L)
  expect_identical(twin_towers$site, seq_len(1500L))
  expect_equal(max(twin_towers$sei) / min(twin_towers$sei), 9)
  expect_equal(1 / (1 + exp(mean(log(twin_towers$vi)))), .7)
  for (dat in list(krw_discrimination, twin_towers)) {
    expect_true(all(is.finite(dat$yi)))
    expect_true(all(is.finite(dat$sei) & dat$sei > 0))
    expect_equal(dat$vi, dat$sei^2)
    expect_true(is.list(attr(dat, "source")))
  }
})

test_that("the school estimates keep their sample, their units and the record of their source", {
  expect_identical(dim(star_schools), c(79L, 12L))
  expect_identical(star_schools$site, seq_len(79L))
  expect_true(all(diff(star_schools$school_id) > 0))
  expect_equal(sum(star_schools$n_students), 5786)
  expect_equal(sum(star_schools$n_classes), 325)
  expect_equal(sum(star_schools$n_classes_small), 127)
  expect_identical(range(star_schools$n_classes_small), c(1L, 4L))
  # Every school has classes and students of both kinds.
  expect_true(all(star_schools$n_classes_small >= 1L &
                    star_schools$n_classes > star_schools$n_classes_small))
  expect_true(all(star_schools$n_students_small >= 1L &
                    star_schools$n_students > star_schools$n_students_small))
  expect_setequal(star_schools$location, c("inner city", "suburban", "rural", "urban"))
  expect_true(all(is.finite(star_schools$yi)))
  expect_true(all(is.finite(star_schools$sei) & star_schools$sei > 0))
  expect_equal(star_schools$vi, star_schools$sei^2)
  # Taking the classes into account widens every standard error.
  expect_true(all(star_schools$sei > star_schools$sei_ols))
  expect_true(all(star_schools$sei_ols > 0))
  source <- attr(star_schools, "source")
  expect_true(is.list(source))
  expect_identical(source$license, "CC0 1.0")
  expect_match(source$data, "10.7910/DVN/SIWH9F", fixed = TRUE)
  expect_match(source$sample, "6,325 students in kindergarten; 5,786 with", fixed = TRUE)
  expect_true(all(c("lme4_version", "r_version") %in% names(source)))
  expect_false(any(grepl("/Users/|/home/", unlist(source))))

  # Both pairs of columns are read by name, and a grid follows from either.
  dat <- as_bef_data(star_schools)
  expect_identical(dat$theta_hat, star_schools$yi)
  expect_identical(dat$sigma, star_schools$sei)
  ols <- as_bef_data(star_schools, theta_hat = "yi_ols", sigma = "sei_ols")
  expect_identical(ols$theta_hat, star_schools$yi_ols)
  expect_identical(ols$sigma, star_schools$sei_ols)
  grid <- make_efron_grid(star_schools, grid_method = "paper_realdata")
  expect_s3_class(grid, "bef_grid")
  expect_true(grid$resolution$adequate)
})

test_that("the three stored fits hold their data in order, keep 800 of 12,000 draws and pass the checks of diagnose()", {
  fits <- list(krw_fit, raudenbush_fit, twin_towers_fit)
  inputs <- list(krw_discrimination, raudenbush1985, twin_towers)
  for (i in seq_along(fits)) {
    fit <- fits[[i]]; dat <- inputs[[i]]
    expect_silent(validate_bef_fit_re(fit))
    expect_equal(fit$metadata$data_list$theta_hat, dat$yi)
    expect_equal(fit$metadata$data_list$sigma, dat$sei)
    expect_equal(nobs(fit), nrow(dat))
    expect_null(fit$cmdstan_fit)
    # One draw in fifteen of four chains of 3,000, and the draws of g only.
    expect_equal(posterior::ndraws(fit$draws), 800L)
    storage <- attr(fit$metadata, "storage", exact = TRUE)
    expect_identical(c(storage$iterations, storage$chains, storage$thin), c(200L, 4L, 15L))
    expect_false(storage$site_draws)
    expect_false(any(grepl("^theta_", dimnames(fit$draws)[[3L]])))
    expect_identical(attr(fit$metadata, "sampler_settings", exact = TRUE)$iter_sampling, 3000L)
    # The defaults of this version: the rule "pilot" and its program.
    expect_identical(fit$metadata$grid_method, "pilot")
    expect_identical(attr(fit$metadata, "program", exact = TRUE), "efron_mix")
    expect_null(attr(fit$metadata, "settings", exact = TRUE))
    expect_s3_class(attr(fit$metadata, "grid", exact = TRUE), "bef_grid")
    expect_identical(attr(fit$metadata, "package_version", exact = TRUE),
                     as.character(utils::packageVersion("bayesEfron")))
    expect_equal(unname(coef(fit)), site_summary(fit)$mean)
    expect_equal(confint(fit)$site, seq_len(nrow(dat)))
    expect_true(all(is.finite(confint(fit)$lower)))
    # The diagnostics are those of the 12,000 draws that were run.
    diag <- diagnose(fit)
    expect_lt(max(diag$rhat, na.rm = TRUE), 1.01)
    expect_gte(min(diag$ess_bulk, na.rm = TRUE), 400)
    expect_gte(min(diag$ess_tail, na.rm = TRUE), 400)
    expect_gt(max(diag$ess_bulk, na.rm = TRUE), 800)
    expect_equal(diag$divergences, 0)
    expect_equal(diag$max_treedepth, 0)
    expect_identical(diag$sampler_diagnostics_failed, character())
    expect_identical(diag$sampler_diagnostics_warned, character())
    expect_false(any(diag$grid$ends$status == "warning"))
  }
  # The simulated sites are fitted without their true effects, on a range
  # that the rule narrowed, with eight spline functions; the two small data
  # sets keep the range of version 0.3.0 and six.
  twin <- attr(twin_towers_fit$metadata, "grid", exact = TRUE)
  expect_true(twin$pilot$narrowed)
  expect_lt(twin$pilot$share, 0.36)
  expect_identical(twin$M, 8L)
  for (fit in list(krw_fit, raudenbush_fit)) {
    grid <- attr(fit$metadata, "grid", exact = TRUE)
    expect_false(grid$pilot$narrowed)
    expect_identical(c(grid$L, grid$M), c(101L, 6L))
    expect_identical(grid$grid,
                     make_efron_grid(fit$metadata$data_list$theta_hat,
                                     fit$metadata$data_list$sigma,
                                     grid_method = "paper_realdata")$grid)
  }
  # Each stored fit is small enough to be included in the package.
  sizes <- vapply(list(krw_fit, raudenbush_fit, twin_towers_fit), function(fit) {
    as.numeric(utils::object.size(fit)) / 2^20
  }, numeric(1))
  expect_true(all(sizes < 1.5))
})

test_that("the stored probabilities follow from the stored fits, and the example is about as close to the truth as the reference", {
  evidence <- readRDS(system.file("examples", "lee_sui_examples.rds",
                                  package = "bayesEfron"))
  expect_equal(nrow(evidence$sensitivity), 6L)
  for (target in c("krw", "twin")) {
    fit <- if (target == "krw") krw_fit else twin_towers_fit
    dat <- if (target == "krw") krw_discrimination else twin_towers
    sd <- fit$metadata$data_list
    dm <- .bef_draws_matrix(fit$draws)
    p <- .bef_posterior_probabilities(.bef_draws_vector(dm, "log_g", sd$L),
                                     dat$yi, dat$sei, sd$grid)
    expect_equal(rowSums(p), rep(1, nrow(dat)), tolerance = 1e-10)
    expect_equal(rowSums(p[, sd$grid > 0, drop = FALSE]),
                 evidence$cases[[target]]$positive_probability, tolerance = 1e-10)
    expect_equal(unname(exceedance_prob(fit, 0)),
                 evidence$cases[[target]]$positive_probability_interpolated, tolerance = 1e-10)
  }
  expect_equal(sqrt(mean((coef(twin_towers_fit) - twin_towers$theta_true)^2)),
               evidence$cases$twin$compact_metrics$rmse, ignore_attr = TRUE)

  # The table that sets the fit of the simulated sites beside other ranges.
  comparison <- readRDS(system.file("examples", "twin_towers_comparison.rds",
                                    package = "bayesEfron"))
  table <- comparison$table
  expect_identical(nrow(table), 5L)
  expect_identical(table$uses_true_effects, c(FALSE, FALSE, TRUE, TRUE, FALSE))
  example <- table[1L, ]
  reference <- table[table$fit == comparison$reference, ]
  expect_identical(example$M, 8L)
  expect_identical(example$rule, "pilot")
  # The example recovers the true effects as well as the fit on the range
  # around them does, without using them.
  expect_lt(example$w1, 1.1 * reference$w1)
  expect_lt(abs(example$rmse / reference$rmse - 1), 0.01)
  expect_lt(example$treedepth_share, 0.01)
  expect_equal(example$divergences, 0)
  expect_lt(example$max_rhat_parameters, 1.01)
  # The range of version 0.3.0 without the true effects is the costly one.
  old <- table[table$rule == "paper_realdata", ]
  expect_gt(old$treedepth_share, 0.3)
  expect_gt(old$w1, reference$w1)
})
