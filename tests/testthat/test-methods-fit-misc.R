# Methods added in version 0.4: prediction for new sites, the information
# criteria that are not defined, the generics of posterior, and arguments
# that a method does not have.

added_fit <- function() {
  readRDS(test_path("_fixtures", "short_fit.rds"))
}

test_that("predict() with a data frame conditions each new site separately", {
  fit <- added_fit()
  data <- fit$metadata$data_list
  # Each row is a new site: other rows do not update its posterior.
  same <- predict(fit, newdata = data.frame(estimate = data$theta_hat, se = data$sigma))
  expect_named(same, c("estimate", "se", "mean", "sd", "lower", "upper", "mode"))
  for (i in seq_len(nrow(same))) {
    one <- predict(fit, newdata = data.frame(estimate = data$theta_hat[i], se = data$sigma[i]))
    expect_equal(same[i, c("mean", "sd", "lower", "upper", "mode")],
                 one[, c("mean", "sd", "lower", "upper", "mode")], ignore_attr = TRUE)
  }
  expect_identical(attr(same, "level"), 0.9)
  # A more precise estimate is shrunk less and has a narrower posterior.
  two <- predict(fit, newdata = data.frame(estimate = c(0.5, 0.5), se = c(0.05, 0.5)),
                 level = 0.8)
  expect_lt(abs(two$mean[1L] - 0.5), abs(two$mean[2L] - 0.5))
  expect_lt(two$sd[1L], two$sd[2L])
  expect_identical(attr(two, "level"), 0.8)
  expect_error(predict(fit, newdata = data.frame(yi = 1, sei = 1)), "`estimate` and `se`",
               class = "bef_invalid_args")
  expect_error(predict(fit, newdata = data.frame(estimate = 1, se = -1)),
               class = "bef_invalid_args")
  # Draws of the effect of a new site, and of the estimate of a new site
  # with a standard error of 0.2.
  expect_length(predict(fit, n = 7L, seed = 2L), 7L)
  expect_length(predict(fit, newdata = 0.2, type = "theta_hat", n = 7L, seed = 2L), 7L)
})

test_that("AIC() and BIC() point to loo() and give no number", {
  fit <- added_fit()
  expect_error(AIC(fit), "loo::loo", class = "bef_invalid_args")
  expect_error(BIC(fit), "loo::loo", class = "bef_invalid_args")
  expect_s3_class(logLik(fit), "logLik")
})

test_that("ndraws() and variables() describe the draws of the fit", {
  fit <- added_fit()
  expect_identical(posterior::ndraws(fit), 100L)
  expect_identical(posterior::variables(fit), dimnames(fit$draws)[[3L]])
})

test_that("an argument that summary() or plot() does not have is an error", {
  fit <- added_fit()
  expect_error(summary(fit, levle = 0.5), "`levle`", class = "bef_invalid_args")
  expect_error(summary(fit, 0.9, 1), "without a name", class = "bef_invalid_args")
  expect_s3_class(summary(fit, level = 0.5), "summary.bef_fit_re")
  expect_error(plot(fit, sortby = "none"), "`sortby`", class = "bef_invalid_args")
})

test_that("print() says how precise the estimates are and how fine the grid is", {
  lines <- format(added_fit(), use_cli = FALSE)
  expect_true(any(grepl("^Standard errors: 0\\.15 to 0\\.24 \\(largest / smallest 1\\.6\\)$", lines)))
  expect_true(any(grepl("^Grid: 51 points from -0\\.62 to 1\\.02; smallest standard error 4\\.573 grid spacings$", lines)))
  expect_true(any(grepl("^Next: summary\\(\\) and site_summary\\(\\)", lines)))
  expect_true(any(grepl("^Draws: 100 \\(1 chain\\); site-level draws kept$", lines)))
})
