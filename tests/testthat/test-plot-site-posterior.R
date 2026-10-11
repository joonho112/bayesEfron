# plot_site_posterior(): the posterior distribution of a site is that of the
# engine, and the likelihood and g are what they are said to be.

site_engine <- function(fit) getFromNamespace(".bef_engine", "bayesEfron")(fit)
site_weights <- function(fit, sites) {
  getFromNamespace(".bef_site_weights", "bayesEfron")(site_engine(fit), sites = sites)
}

test_that("the three curves of a site are drawn", {
  fit <- plot_fit(12L, labels = TRUE)
  engine <- site_engine(fit)
  plot <- plot_site_posterior(fit, sites = c(3, 7))
  expect_s3_class(plot, "ggplot")
  expect_identical(unname(plot_geoms(plot)), c("GeomArea", "GeomLine"))
  expect_null(plot$labels$title)
  expect_identical(plot$labels$y, "Density")
  expect_identical(levels(plot$data$curve),
                   c("Posterior of the effect", "Likelihood of the estimate", "Estimate of g"))
  expect_identical(nlevels(plot$data$panel), 2L)
  expect_match(levels(plot$data$panel)[[1L]], "^Site 003: estimate .*, standard error ")

  data <- plot$data
  first <- data[as.integer(data$panel) == 1L, ]
  posterior <- first[first$curve == "Posterior of the effect", ]
  weights <- site_weights(fit, 3L)
  at <- match(posterior$value, engine$grid)
  expect_false(anyNA(at))
  expect_equal(posterior$density, weights[1L, at] / engine$spacing)
  # The posterior distribution integrates to one over the part that is drawn.
  expect_equal(sum(posterior$density) * engine$spacing, 1, tolerance = 1e-3)
  g <- first[first$curve == "Estimate of g", ]
  expect_equal(g$density, g_summary(fit)$mean[at])
  likelihood <- first[first$curve == "Likelihood of the estimate", ]
  expect_equal(likelihood$density,
               stats::dnorm(likelihood$value, engine$theta_hat[[3L]], engine$sigma[[3L]]))
  # The axis holds three standard errors on each side of both estimates.
  expect_lte(min(data$value), min(engine$theta_hat[c(3L, 7L)] - 3 * engine$sigma[c(3L, 7L)]))
  expect_gte(max(data$value), max(engine$theta_hat[c(3L, 7L)] + 3 * engine$sigma[c(3L, 7L)]))

  scaled <- plot_site_posterior(fit, sites = 3, scale = 100)$data
  one <- plot_site_posterior(fit, sites = 3)$data
  expect_equal(scaled$value, 100 * one$value)
  expect_equal(scaled$density, one$density / 100)
})

test_that("sites are chosen by label, and curves can be left out", {
  fit <- plot_fit(12L, labels = TRUE)
  by_label <- plot_site_posterior(fit, sites = c("Site 007", "Site 003"))
  expect_match(levels(by_label$data$panel), "^Site 00[73]:")
  expect_match(levels(by_label$data$panel)[[1L]], "^Site 007")

  alone <- plot_site_posterior(fit, sites = 2, likelihood = FALSE, prior = FALSE)
  expect_identical(as.character(unique(alone$data$curve)), "Posterior of the effect")
  expect_no_match(plot_caption(alone), "Dashed line")
  expect_no_match(plot_caption(alone), "Long dashes")
  no_prior <- plot_site_posterior(fit, sites = 2, prior = FALSE)
  expect_identical(as.character(unique(no_prior$data$curve)),
                   c("Posterior of the effect", "Likelihood of the estimate"))

  # A fit whose sites have no labels names them by their numbers.
  plain <- plot_site_posterior(plot_fit(12L), sites = 5)
  expect_match(levels(plain$data$panel), "^Site 5: estimate ")
})

test_that("with normal components the posterior distribution is a mixture of normals", {
  fit <- plot_fit(12L, component = "normal")
  engine <- site_engine(fit)
  expect_gt(engine$tau_c, 0)
  data <- plot_site_posterior(fit, sites = 4)$data
  posterior <- data[data$curve == "Posterior of the effect", ]
  # It integrates to one, and its mean is the posterior mean of the site.
  step <- diff(posterior$value[1:2])
  expect_equal(sum(posterior$density) * step, 1, tolerance = 1e-3)
  expect_equal(sum(posterior$value * posterior$density) * step, unname(coef(fit)[[4L]]),
               tolerance = 1e-3)
  g <- data[data$curve == "Estimate of g", ]
  expect_equal(sum(g$density) * step, 1, tolerance = 5e-3)
})

test_that("a posterior distribution with two modes is seen in the stored fit", {
  data <- plot_site_posterior(bayesEfron::twin_towers_fit, sites = 268)$data
  posterior <- data[data$curve == "Posterior of the effect", ]
  # Probability on both sides of the gap between the two groups, and
  # little in it.
  expect_gt(max(posterior$density[posterior$value < -0.5]), 0.3)
  expect_gt(max(posterior$density[posterior$value > 0.5]), 0.2)
  expect_lt(min(posterior$density[abs(posterior$value) < 0.5]), 0.05)
})

test_that("wrong arguments are errors", {
  fit <- plot_fit(150L)
  expect_error(plot_site_posterior(fit), "must name the sites", class = "bef_invalid_args")
  expect_error(plot_site_posterior(fit, sites = 1:13), "at most 12", class = "bef_invalid_args")
  expect_error(plot_site_posterior(fit, sites = 151), class = "bef_invalid_args")
  expect_error(plot_site_posterior(fit, sites = "none"), class = "bef_invalid_args")
  expect_error(plot_site_posterior(fit, sites = 1, likelihood = NA), class = "bef_invalid_args")
  expect_error(plot_site_posterior(fit, sites = 1, level = 0.9), "no argument",
               class = "bef_invalid_args")
  expect_s3_class(plot_site_posterior(fit, sites = 1:12), "ggplot")
})
