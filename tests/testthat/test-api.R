# User-facing behavior that is not tied to one file in R/: data frame input,
# the grid comparison table, refusal of fits saved by version 0.1, and the
# one-number views of per-parameter diagnostics.

api_ns <- function(name) getFromNamespace(name, "bayesEfron")

example_fit <- function() {
  readRDS(test_path("_fixtures", "short_fit.rds"))
}

test_that("as_bef_data.data.frame infers conventional column names", {
  df <- data.frame(theta_hat = seq(-0.4, 0.4, length.out = 8),
                   sigma = rep(0.3, 8))
  out <- as_bef_data(df)
  expect_s3_class(out, "bef_data")
  expect_length(out$theta_hat, 8L)
  expect_identical(out$source, "data.frame")
})

test_that("as_bef_data.data.frame accepts metafor column conventions", {
  df <- data.frame(yi = seq(-0.4, 0.4, length.out = 8), sei = rep(0.3, 8))
  expect_length(as_bef_data(df)$theta_hat, 8L)

  vi <- data.frame(yi = seq(-0.4, 0.4, length.out = 8), vi = rep(0.09, 8))
  out <- as_bef_data(vi, sigma = "vi", variance = TRUE)
  expect_equal(out$sigma, rep(0.3, 8), tolerance = 1e-12)
})

test_that("as_bef_data.data.frame carries site labels", {
  df <- data.frame(theta_hat = seq(-0.4, 0.4, length.out = 8),
                   sigma = rep(0.3, 8),
                   study = paste0("s", 1:8))
  expect_identical(as_bef_data(df, site = "study")$names, paste0("s", 1:8))
})

test_that("as_bef_data.data.frame refuses missing values rather than dropping them", {
  df <- data.frame(theta_hat = c(seq(-0.4, 0.4, length.out = 7), NA),
                   sigma = rep(0.3, 8))
  err <- tryCatch(as_bef_data(df), error = identity)
  expect_s3_class(err, "bef_invalid_args")
  expect_match(conditionMessage(err), "missing values")
})

test_that("as_bef_data.data.frame refuses uninferable columns", {
  df <- data.frame(a = 1:8, b = 1:8)
  err <- tryCatch(as_bef_data(df), error = identity)
  expect_s3_class(err, "bef_invalid_args")
})

test_that("the grid comparison table is ordered and flags unconverged fits", {
  skip_if_not_installed("loo")
  # A synthetic table: only the ordering and the printed note are checked here.
  cmp <- structure(
    data.frame(
      L = c(101L, 51L), M = c(6L, 4L), expansion = c(0.5, 0.5),
      elpd_loo = c(-334.1, -340.2), se_elpd_loo = c(9.3, 9.5),
      p_loo = c(2.4, 1.8), elpd_diff = c(0, -6.1), se_diff = c(0, 2.1),
      max_rhat = c(1.01, 1.15), min_ess_bulk = c(500, 12),
      divergences = c(0, 0), max_pareto_k = c(0.25, 0.4),
      sampler_seconds = c(0.9, 0.8)
    ),
    class = c("bef_grid_comparison", "data.frame")
  )
  expect_s3_class(cmp, "bef_grid_comparison")
  expect_true(all(diff(cmp$elpd_loo) <= 0))
  out <- capture.output(print(cmp))
  # A specification that did not converge must be called out next to its score.
  expect_true(any(grepl("convergence or PSIS problems", out)))
})

test_that("differences in elpd_loo follow the order of the fits, whatever loo_compare() calls them", {
  skip_if_not_installed("loo")
  set.seed(11)
  fake_loo <- function(shift) {
    log_lik <- matrix(rnorm(400 * 19, mean = -1 + shift, sd = 0.3), 400, 19)
    suppressWarnings(loo::loo(log_lik))
  }
  # The best fit is the second one, the worst the first.
  loos <- list(fake_loo(-0.5), fake_loo(0.3), fake_loo(0))
  elpd <- vapply(loos, function(x) x$estimates["elpd_loo", "Estimate"], numeric(1))
  expect_identical(order(-elpd), c(2L, 3L, 1L))

  out <- api_ns(".bef_elpd_differences")(loos)
  expect_equal(out$elpd_diff, elpd - max(elpd), tolerance = 1e-10)
  expect_identical(out$elpd_diff[[2L]], 0)
  expect_identical(out$se_diff[[2L]], 0)
  expect_true(all(out$se_diff[c(1L, 3L)] > 0))

  # The same numbers as loo::loo_compare(), which reports them best first.
  compared <- as.data.frame(loo::loo_compare(loos))
  expect_equal(sort(out$elpd_diff, decreasing = TRUE), compared$elpd_diff,
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(out$se_diff[order(-out$elpd_diff)], compared$se_diff,
               tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("loo() refuses a fit saved by version 0.1 and says to refit", {
  skip_if_not_installed("loo")
  fake <- structure(
    list(
      draws = array(0, dim = c(2L, 1L, 1L),
                    dimnames = list(NULL, NULL, "mean_g")),
      metadata = stats::setNames(
        vector("list", length(api_ns(".bef_v01_metadata_fields")())),
        api_ns(".bef_v01_metadata_fields")()
      ),
      posterior = list()
    ),
    class = c("bef_fit_re", "bef_fit")
  )
  expect_true(api_ns(".bef_is_v01_fit")(fake))
  err <- tryCatch(loo::loo(fake), error = identity)
  expect_s3_class(err, "bef_invalid_fit")
  expect_match(conditionMessage(err), "bayesEfron 0\\.1")
  expect_match(conditionMessage(err), "refit", ignore.case = TRUE)
})

test_that("chain ids match the draws matrix layout", {
  draws <- array(0, dim = c(3L, 2L, 1L),
                 dimnames = list(NULL, NULL, "x"))
  expect_identical(
    api_ns(".bef_draw_chain_id")(draws),
    c(1L, 1L, 1L, 2L, 2L, 2L)
  )
})

# `rhat`, `ess_bulk` and `ess_tail` are per-parameter vectors, and
# `divergences` and `max_treedepth` have one entry per chain, but print(),
# summary() and the diagnostic plot each show one number per diagnostic. The
# formatter returns "NA" for any input whose length is not 1, so a vector that
# is not reduced first prints NA without an error.

test_that(".bef_diagnostic_worst reduces per-parameter vectors", {
  worst <- api_ns(".bef_diagnostic_worst")
  expect_identical(worst(c(a = 1.00, b = 1.15, c = 1.02), max), 1.15)
  expect_identical(worst(c(a = 900, b = 42, c = 500), min), 42)
  expect_identical(worst(c(3L, 0L, 1L), sum), 4)
  # Non-finite entries are dropped, not propagated: rhat is legitimately NA
  # for constant variables.
  expect_identical(worst(c(1.01, NA, Inf, 1.20), max), 1.20)
  expect_identical(worst(c(NA_real_, NA_real_), max), NA_real_)
  expect_identical(worst(numeric(0), min), NA_real_)
})

test_that("text views report the worst per-parameter diagnostic, not NA", {
  fit <- example_fit()
  diag <- diagnose(fit)

  # The example fit has to have more than one monitored parameter for this
  # test to mean anything.
  expect_gt(length(diag$rhat), 1L)
  expect_true(any(is.finite(diag$rhat)))

  expected_rhat <- max(diag$rhat[is.finite(diag$rhat)])
  expected_ess <- min(diag$ess_bulk[is.finite(diag$ess_bulk)])

  summary_txt <- paste(
    utils::capture.output(print(summary(fit))),
    collapse = "\n"
  )
  print_txt <- paste(utils::capture.output(print(fit)), collapse = "\n")

  for (txt in list(summary_txt, print_txt)) {
    expect_false(grepl("Rhat[^\n]*NA", txt))
    expect_false(grepl("ESS bulk[^\n]*NA", txt))
    expect_false(grepl("ESS tail[^\n]*NA", txt))
  }
  expect_match(summary_txt, format(signif(expected_rhat, 4L), trim = TRUE),
               fixed = TRUE)
  expect_match(summary_txt, format(signif(expected_ess, 4L), trim = TRUE),
               fixed = TRUE)

  # The plot must use the same reduction as the two text views.
  plot_data <- api_ns(".bef_diagnostic_plot_data")(diag)
  expect_equal(plot_data$value[plot_data$metric == "rhat"], expected_rhat)
  expect_equal(plot_data$value[plot_data$metric == "ess_bulk"], expected_ess)
})

test_that("sampler counts from several chains are added, not dropped", {
  diag <- structure(
    list(
      rhat = c(`alpha[1]` = 1.001, `alpha[2]` = 1.004),
      ess_bulk = c(`alpha[1]` = 1200, `alpha[2]` = 950),
      ess_tail = c(`alpha[1]` = 1100, `alpha[2]` = 880),
      divergences = c(0L, 2L, 0L, 1L),
      max_treedepth = c(0L, 0L, 0L, 0L)
    ),
    class = "bef_diagnostic"
  )
  plot_data <- api_ns(".bef_diagnostic_plot_data")(diag)
  expect_equal(plot_data$value[plot_data$metric == "divergences"], 3)
  expect_equal(plot_data$value[plot_data$metric == "rhat"], 1.004)
  expect_equal(plot_data$value[plot_data$metric == "ess_tail"], 880)
})

test_that("effective_params is large when there is little shrinkage", {
  # The quantity is sum_i Var(theta_i | data) / sigma_i^2. It is zero under
  # complete shrinkage and grows as the site posteriors approach their
  # likelihoods. It is bounded below by zero but not above by K: under a
  # multimodal g, a site whose likelihood falls between two modes can have a
  # posterior variance larger than sigma_i^2.
  fit <- example_fit()
  K <- as.integer(fit$metadata$data_list$K)
  ep <- fit$posterior$effective_params
  expect_true(all(is.finite(ep)))
  expect_true(all(ep >= 0))
  # Five nearly independent sites under weak shrinkage sit near K.
  expect_gt(mean(ep), K / 2)
})
