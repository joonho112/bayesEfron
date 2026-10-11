# Divergent transitions, the maximum tree depth and E-BFMI are checked
# whatever the length of the run. R-hat and the effective sample sizes are
# compared with their thresholds only from 400 draws.
test_that("a run of 100 draws is checked for divergences, tree depth and E-BFMI, and not for Rhat or ESS", {
  check <- getFromNamespace(".bef_diagnostic_check", "bayesEfron")
  values <- list(
    rhat = 2.257, ess_bulk = 1.377, ess_tail = 1.894,
    divergences = 14, max_treedepth = 86, ebfmi = c(0.15, 0.4)
  )
  draws <- array(0, dim = c(100L, 1L, 1L))
  result <- check(values, draws)
  expect_identical(result$failed, "divergences")
  expect_setequal(result$warned, c("max_treedepth", "ebfmi"))
  expect_identical(
    result$skipped, c("rhat_check", "ess_bulk_check", "ess_tail_check")
  )
})

test_that("series thresholds apply once 400 draws are available", {
  check <- getFromNamespace(".bef_diagnostic_check", "bayesEfron")
  values <- list(
    rhat = c(1.005, 1.06), ess_bulk = c(500, 250), ess_tail = 300,
    divergences = 0, max_treedepth = 0, ebfmi = 0.3
  )
  result <- check(values, array(0, dim = c(200L, 2L, 1L)))
  expect_identical(result$failed, "rhat")
  expect_setequal(result$warned, c("ess_bulk", "ess_tail"))
  expect_identical(result$skipped, character())
})

test_that("a diagnostic whose check was skipped may keep its value, and one recorded as skipped may not", {
  validate <- getFromNamespace(".bef_validate_diagnostic_skip_consistency", "bayesEfron")
  values <- list(rhat = 2.2, ess_bulk = 2, ess_tail = 3)
  expect_error(
    validate(
      values, c("rhat_check", "ess_bulk_check", "ess_tail_check"),
      "bef_invalid_fit"
    ),
    NA
  )
  expect_error(validate(values, "rhat", "bef_invalid_fit"), class = "bef_invalid_fit")
})

test_that("postprocessing stores failures even for a short fit", {
  draws <- posterior::as_draws_array(array(
    seq_len(100L), dim = c(100L, 1L, 1L),
    dimnames = list(NULL, "chain1", "alpha[1]")
  ))
  sampler <- list(diagnostic_summary = function(...) list(
    num_divergent = 14, num_max_treedepth = 86, ebfmi = 0.15
  ))
  expect_warning(
    result <- getFromNamespace(".bef_postprocess_diagnostics", "bayesEfron")(
      sampler, draws
    ),
    class = "bef_sampler_diagnostics_failed"
  )
  expect_identical(result$sampler_diagnostics_failed, "divergences")
  expect_setequal(result$sampler_diagnostics_warned, c("max_treedepth", "ebfmi"))
  expect_true(all(c("rhat_check", "ess_bulk_check", "ess_tail_check") %in%
    result$diagnostic_skipped))
  expect_identical(result$values$divergences, 14)
  expect_identical(result$values$max_treedepth, 86)
})

test_that("an infinite Rhat is stored as missing with the name of its variable, and Rhat is not recorded as skipped", {
  local_mocked_bindings(
    summarize_draws = function(...) data.frame(
      variable = c("alpha[1]", "alpha[2]"), rhat = c(Inf, 1.01),
      ess_bulk = c(3, 4), ess_tail = c(3, 4)
    ),
    .package = "posterior"
  )
  compute <- getFromNamespace(".bef_per_parameter_diagnostics", "bayesEfron")
  # An infinite Rhat means that the chains did not mix: it is recorded as a
  # failed check and not as a check that could not be made.
  short <- compute(array(0, dim = c(4L, 1L, 2L)))
  expect_identical(unname(short$rhat), c(NA, 1.01))
  expect_identical(short$rhat_infinite, "alpha[1]")
  expect_identical(short$diagnostic_skipped, character())
  expect_identical(unname(short$ess_bulk), c(3, 4))
  expect_error(
    getFromNamespace(".bef_validate_diagnostic_skip_consistency", "bayesEfron")(
      short, short$diagnostic_skipped, "bef_invalid_fit"
    ),
    NA
  )
  # With enough draws for the checks the result is the same: no warning that
  # the diagnostic is unavailable, because it is available and it failed.
  expect_warning(long <- compute(array(0, dim = c(400L, 1L, 2L))), NA)
  expect_identical(unname(long$rhat), c(NA, 1.01))
  expect_identical(long$rhat_infinite, "alpha[1]")
})

test_that("E-BFMI below 0.3 is a warning, also in a stored fit that recorded none", {
  check <- getFromNamespace(".bef_diagnostic_check", "bayesEfron")
  values <- list(rhat = c(a = 1.001), ess_bulk = c(a = 900), ess_tail = c(a = 900),
                 divergences = 0, max_treedepth = 0, ebfmi = c(0.25, 0.6))
  result <- check(values, array(0, dim = c(500L, 2L, 1L)))
  expect_identical(result$failed, character())
  expect_identical(result$warned, "ebfmi")
  values$ebfmi <- c(0.31, 0.6)
  expect_identical(check(values, array(0, dim = c(500L, 2L, 1L)))$warned, character())

  # As in a fit stored by a version before 0.4.0, whose threshold was 0.2: a
  # smallest E-BFMI of 0.26 and no warning recorded. The present threshold
  # applies.
  stored <- bayesEfron::raudenbush_fit
  expect_false("ebfmi" %in% diagnose(stored)$sampler_diagnostics_warned)
  recorded <- attr(stored$metadata, "diagnostics", exact = TRUE)
  recorded$ebfmi[[1L]] <- 0.26
  attr(stored$metadata, "diagnostics") <- recorded
  diag <- diagnose(stored)
  expect_true("ebfmi" %in% diag$sampler_diagnostics_warned)
  expect_true(any(grepl("ebfmi", format(stored, use_cli = FALSE))))
  table <- as.data.frame(diag)
  expect_named(table, c("diagnostic", "value", "threshold", "rule", "status"))
  expect_identical(table$diagnostic,
                   c("rhat", "ess_bulk", "ess_tail", "divergences", "max_treedepth", "ebfmi",
                     "grid_ends"))
  expect_identical(table$status[table$diagnostic == "ebfmi"], "warning")
  expect_identical(table$threshold[table$diagnostic == "ebfmi"], 0.3)
  expect_identical(table$status[table$diagnostic == "divergences"], "ok")
})

test_that("a fit made before 0.4 has no sampler quantities for bayesplot, and says so", {
  skip_if_not_installed("bayesplot")
  fit <- readRDS(test_path("_fixtures", "short_fit.rds"))
  expect_error(bayesplot::nuts_params(fit), "Fit the model again", class = "bef_invalid_fit")
  expect_named(bayesplot::log_posterior(fit), c("Chain", "Iteration", "Value"))
  expect_identical(nrow(bayesplot::log_posterior(fit)), 100L)
  expect_identical(names(bayesplot::rhat(fit)), dimnames(fit$draws)[[3L]])
  expect_equal(bayesplot::neff_ratio(fit, pars = "lambda"),
               attr(fit$metadata, "diagnostics")$ess_bulk["lambda"] / 100)
})

test_that("the effective sample size of a thinned fit is divided by the length of the whole run", {
  skip_if_not_installed("bayesplot")
  # The stored fit keeps 800 of 12,000 draws; its diagnostics are those of
  # the whole run, and so is the ratio.
  fit <- bayesEfron::krw_fit
  expect_equal(prod(dim(fit$draws)[1:2]), 800)
  ratio <- bayesplot::neff_ratio(fit)
  expect_equal(ratio, diagnose(fit)$ess_bulk / 12000)
  expect_lt(max(ratio), 1.5)
  expect_identical(nrow(bayesplot::nuts_params(fit, pars = "energy__")), 800L)
})
