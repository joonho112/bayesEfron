# predict(), fitted() and waic() on a stored fit. None of them runs Stan.

stored_fit <- function() {
  readRDS(test_path("_fixtures", "short_fit.rds"))
}

test_that("predict() draws effects of a new site from the grid", {
  fit <- stored_fit()
  grid <- fit$metadata$data_list$grid
  n_draws <- posterior::ndraws(fit$draws)

  theta_new <- predict(fit, seed = 3)
  expect_length(theta_new, n_draws)
  expect_true(all(theta_new %in% grid))
  expect_identical(predict(fit, seed = 3), theta_new)
  expect_false(identical(predict(fit, seed = 4), theta_new))

  expect_length(predict(fit, n = 7, seed = 3), 7L)
  expect_length(predict(fit, n = 3L * n_draws, seed = 3), 3L * n_draws)

  spacing <- diff(grid)[1]
  jittered <- predict(fit, jitter = TRUE, seed = 3)
  nearest <- vapply(jittered, function(x) min(abs(x - grid)), numeric(1))
  expect_true(all(nearest <= spacing / 2 + 1e-12))
  expect_false(all(jittered %in% grid))
})

test_that("predict(type = \"theta_hat\") adds the sampling error of the new site", {
  fit <- stored_fit()
  n_draws <- posterior::ndraws(fit$draws)

  one <- predict(fit, newdata = 0.2, type = "theta_hat", seed = 3)
  expect_length(one, n_draws)
  expect_false(all(one %in% fit$metadata$data_list$grid))

  several <- predict(fit, newdata = c(small = 0.1, large = 0.5), type = "theta_hat",
                     n = 2000, seed = 3)
  expect_identical(dim(several), c(2000L, 2L))
  expect_identical(colnames(several), c("small", "large"))
  expect_lt(stats::sd(several[, "small"]), stats::sd(several[, "large"]))

  expect_error(predict(fit, type = "theta_hat"), class = "bef_invalid_args")
  expect_error(predict(fit, newdata = -1, type = "theta_hat"), class = "bef_invalid_args")
  expect_error(predict(fit, n = 0), class = "bef_invalid_args")
  expect_error(predict(fit, seed = 1.5), class = "bef_invalid_args")
  expect_error(predict(fit, 0.2, "theta_hat"), class = "bef_error")
})

test_that("predict(seed = ) leaves the random number generator as it found it", {
  fit <- stored_fit()
  withr::local_preserve_seed()

  set.seed(11)
  before <- get(".Random.seed", envir = globalenv())
  predict(fit, seed = 5)
  expect_identical(get(".Random.seed", envir = globalenv()), before)

  rm(".Random.seed", envir = globalenv())
  predict(fit, seed = 5)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("one predicted estimate keeps the documented vector or matrix shape", {
  fit <- stored_fit()
  withr::local_preserve_seed()
  set.seed(11)
  before <- get(".Random.seed", envir = globalenv())

  named <- predict(fit, newdata = c(a = 0.1, b = 0.3), type = "theta_hat",
                   n = 1, seed = 9)
  unnamed <- predict(fit, newdata = c(0.1, 0.3), type = "theta_hat",
                     n = 1, seed = 9)
  one <- predict(fit, newdata = 0.1, type = "theta_hat", n = 1, seed = 9)

  expect_identical(dim(named), c(1L, 2L))
  expect_identical(colnames(named), c("a", "b"))
  expect_identical(dim(unnamed), c(1L, 2L))
  expect_null(colnames(unnamed))
  expect_equal(as.numeric(named), as.numeric(unnamed))
  expect_true(all(is.finite(named)))
  expect_identical(
    predict(fit, newdata = c(a = 0.1, b = 0.3), type = "theta_hat",
            n = 1, seed = 9),
    named
  )
  expect_type(one, "double")
  expect_length(one, 1L)
  expect_null(dim(one))
  expect_identical(get(".Random.seed", envir = globalenv()), before)
})

test_that("fitted() returns the posterior means of the site effects", {
  fit <- stored_fit()
  expect_identical(fitted(fit), coef(fit))
  expect_identical(fitted(fit, "ignored"), coef(fit))
  expect_length(fitted(fit), nobs(fit))
})

test_that("loo() and waic() use the log density of each site", {
  skip_if_not_installed("loo")
  fit <- stored_fit()
  K <- nobs(fit)

  fit_loo <- suppressWarnings(loo::loo(fit))
  expect_s3_class(fit_loo, "psis_loo")
  expect_identical(nrow(fit_loo$pointwise), K)

  fit_waic <- suppressWarnings(loo::waic(fit))
  expect_s3_class(fit_waic, "waic")
  expect_identical(nrow(fit_waic$pointwise), K)

  log_lik <- posterior::as_draws_matrix(fit$draws)[, sprintf("log_lik[%d]", seq_len(K))]
  direct <- suppressWarnings(loo::waic(matrix(as.numeric(log_lik), nrow = nrow(log_lik))))
  expect_equal(fit_waic$estimates["elpd_waic", "Estimate"],
               direct$estimates["elpd_waic", "Estimate"])
})
