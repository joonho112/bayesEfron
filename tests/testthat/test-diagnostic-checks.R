# Sampler events remain meaningful even when a run is too short for
# threshold checks on R-hat and effective sample sizes.
test_that("short runs flag sampler events and identify unchecked series", {
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

test_that("short-run unchecked thresholds do not hide numerical diagnostics", {
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

test_that("nonfinite convergence diagnostics obey the unavailable-field contract", {
  local_mocked_bindings(
    summarize_draws = function(...) data.frame(
      variable = c("alpha[1]", "alpha[2]"), rhat = c(Inf, 1.01),
      ess_bulk = c(3, 4), ess_tail = c(3, 4)
    ),
    .package = "posterior"
  )
  compute <- getFromNamespace(".bef_per_parameter_diagnostics", "bayesEfron")
  short <- compute(array(0, dim = c(4L, 1L, 2L)))
  expect_true(all(is.na(short$rhat)))
  expect_identical(short$diagnostic_skipped, "rhat")
  expect_identical(unname(short$ess_bulk), c(3, 4))
  expect_error(
    getFromNamespace(".bef_validate_diagnostic_skip_consistency", "bayesEfron")(
      short, short$diagnostic_skipped, "bef_invalid_fit"
    ),
    NA
  )
  expect_warning(
    long <- compute(array(0, dim = c(400L, 1L, 2L))),
    class = "bef_diagnostic_skipped"
  )
  expect_true(all(is.na(long$rhat)))
})
