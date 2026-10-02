# Synthetic log likelihoods let the tests change the density's units while
# keeping posterior draws and importance weights fixed. No sampler is run.
loo_stability_case <- function() {
  withr::local_seed(420)
  log_lik <- matrix(stats::rnorm(400L * 5L, -3, 0.25), 400L, 5L)
  colnames(log_lik) <- sprintf("log_lik[%d]", seq_len(ncol(log_lik)))
  fit <- readRDS(test_path("_fixtures", "short_fit.rds"))
  fit$draws <- posterior::as_draws_array(array(
    as.numeric(log_lik), dim = c(200L, 2L, 5L),
    dimnames = list(NULL, NULL, colnames(log_lik))
  ))
  list(fit = fit, log_lik = log_lik, chain_id = rep(1:2, each = 200L))
}

test_that("LOO centering preserves ordinary results and forwarded methods", {
  skip_if_not_installed("loo")
  case <- loo_stability_case()
  r_eff <- loo::relative_eff(exp(case$log_lik), chain_id = case$chain_id)
  for (method in c("psis", "tis", "sis")) {
    for (save in c(FALSE, TRUE)) {
      direct <- loo::loo(case$log_lik, r_eff = r_eff, cores = 1,
                         save_psis = save, is_method = method)
      result <- loo::loo(case$fit, cores = 1, save_psis = save,
                         is_method = method)
      expect_equal(result, direct, tolerance = 1e-9)
    }
  }
})

test_that("LOO restores shifted scores, standard errors and saved weights", {
  skip_if_not_installed("loo")
  case <- loo_stability_case()
  offsets <- c(-1000, -1100, -1200, -1300, -1400)
  shifted <- case$fit
  shifted$draws <- posterior::as_draws_array(array(
    as.numeric(sweep(case$log_lik, 2L, offsets, "+")),
    dim = dim(case$fit$draws), dimnames = dimnames(case$fit$draws)
  ))
  expect_true(all(exp(as.numeric(shifted$draws)) == 0))
  r_eff <- loo::relative_eff(exp(case$log_lik), chain_id = case$chain_id)

  for (method in c("psis", "tis", "sis")) {
    for (save in c(FALSE, TRUE)) {
      direct <- loo::loo(case$log_lik, r_eff = r_eff, cores = 1,
                         save_psis = save, is_method = method)
      result <- loo::loo(shifted, cores = 1, save_psis = save,
                         is_method = method)
      expected <- direct$pointwise
      expected[, "elpd_loo"] <- expected[, "elpd_loo"] + offsets
      expected[, "looic"] <- -2 * expected[, "elpd_loo"]

      expect_equal(result$pointwise, expected, tolerance = 1e-9)
      expect_true(all(is.finite(result$pointwise[, "mcse_elpd_loo"])))
      expect_true(all(is.finite(result$diagnostics$n_eff)))
      expect_equal(result$diagnostics, direct$diagnostics, tolerance = 1e-9)
      expect_identical(class(result), class(direct))
      expect_identical(dim(result), dim(direct))
      for (metric in c("elpd_loo", "p_loo", "looic")) {
        values <- expected[, metric]
        expect_equal(result$estimates[metric, "Estimate"], sum(values),
                     tolerance = 1e-9)
        expect_equal(result$estimates[metric, "SE"],
                     sqrt(nrow(expected) * stats::var(values)),
                     tolerance = 1e-9)
      }
      # Check deprecated aliases without invoking loo's deprecation warning.
      aliases <- c("elpd_loo", "p_loo", "looic", "se_elpd_loo",
                   "se_p_loo", "se_looic")
      expect_equal(unlist(unclass(result)[aliases], use.names = FALSE),
                   as.vector(result$estimates))

      slot <- paste0(method, "_object")
      if (save) {
        object <- result[[slot]]
        original <- direct[[slot]]
        for (log in c(FALSE, TRUE)) {
          expect_equal(stats::weights(object, log = log, normalize = TRUE),
                       stats::weights(original, log = log, normalize = TRUE),
                       tolerance = 1e-9)
        }
        expect_equal(
          stats::weights(object, log = TRUE, normalize = FALSE),
          sweep(stats::weights(original, log = TRUE, normalize = FALSE),
                2L, offsets, "-"),
          tolerance = 1e-9
        )
        expect_equal(attr(object, "norm_const_log"),
                     attr(original, "norm_const_log") - offsets,
                     tolerance = 1e-9)
        expect_identical(class(object), class(original))
        expect_identical(dim(object), dim(original))
        expect_equal(attr(object, "r_eff"), attr(original, "r_eff"),
                     tolerance = 1e-9)
      } else {
        expect_null(result[[slot]])
      }
    }
  }
})

test_that("LOO still rejects missing or nonfinite log likelihoods", {
  skip_if_not_installed("loo")
  case <- loo_stability_case()
  missing <- case$fit
  dimnames(missing$draws)[[3L]][1L] <- "not_log_lik"
  expect_error(loo::loo(missing), class = "bef_invalid_fit")

  for (value in c(NA_real_, Inf, -Inf)) {
    invalid <- case$fit
    invalid$draws[1L, 1L, 1L] <- value
    expect_error(loo::loo(invalid), "is not finite",
                 class = "bef_invalid_fit")
  }
})
