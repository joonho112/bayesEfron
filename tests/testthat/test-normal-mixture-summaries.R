# Numeric references here use the mixture formula directly, independently
# of the posterior engine's density, CDF and quantile helpers.
normal_summary_case <- function(mean, sd, weight) {
  engine <- .bef_engine_from_parts(matrix(log(weight), nrow = 1L),
                                   grid = 2 * mean, theta_hat = 0,
                                   sigma = sqrt(2) * sd, tau_c = sqrt(2) * sd)
  list(engine = engine, weights = matrix(weight / sum(weight), nrow = 1L),
       mean = mean, sd = sd, weight = weight / sum(weight))
}

normal_summary_cdf <- function(x, case) {
  vapply(x, function(value) sum(case$weight * pnorm((value - case$mean) / case$sd)),
         numeric(1L))
}

normal_summary_mode_reference <- function(case) {
  # Search separately around every component so even a remote, narrow peak
  # gets the same resolution as the dominant one.
  origin <- case$mean[which.max(case$weight)]
  mu <- (case$mean - origin) / case$sd
  density <- function(z) sum(case$weight * dnorm(z - mu))
  x <- sort(unique(as.vector(outer(mu, seq(-2, 2, length.out = 4001L), `+`))))
  d <- vapply(x, density, numeric(1L))
  peaks <- which(d > 0 & d >= c(-Inf, head(d, -1L)) & d >= c(tail(d, -1L), -Inf))
  candidates <- vapply(peaks, function(i) {
    a <- x[max(1L, i - 1L)]
    b <- x[min(length(x), i + 1L)]
    if (a == b) return(a)
    a + optimize(function(delta) density(a + delta), c(0, b - a),
                 maximum = TRUE, tol = 1e-12)$maximum
  }, numeric(1L))
  origin + case$sd * candidates[which.max(vapply(candidates, density, numeric(1L)))]
}

normal_summary_hpd_reference <- function(case, level) {
  origin <- case$mean[which.max(case$weight)]
  mu <- (case$mean - origin) / case$sd
  cdf <- function(z) sum(case$weight * pnorm(z - mu))
  quantile <- function(p) {
    if (p <= 0) return(-Inf)
    if (p >= 1) return(Inf)
    uniroot(function(z) cdf(z) - p,
            c(min(mu) + qnorm(p), max(mu) + qnorm(p)), tol = 1e-10)$root
  }
  width <- function(u) quantile(u + level) - quantile(u)
  starts <- seq(0, 1 - level, length.out = 2001L)
  widths <- vapply(starts, width, numeric(1L))
  minima <- which(widths <= c(Inf, head(widths, -1L)) &
                    widths <= c(tail(widths, -1L), Inf) & is.finite(widths))
  best <- vapply(minima, function(i) {
    optimize(width, starts[c(i - 1L, i + 1L)], tol = 1e-12)$minimum
  }, numeric(1L))
  u <- best[which.min(vapply(best, width, numeric(1L)))]
  origin + case$sd * c(quantile(u), quantile(u + level))
}

test_that("a narrow single normal has its exact mode and shortest interval", {
  for (sd in c(.1, .01, .001, 1e-6)) {
    case <- normal_summary_case(c(-10, .6, 10), sd, c(0, 1, 0))
    expect_equal(.bef_site_mode(case$engine, case$weights), .6, tolerance = 1e-12)
    limits <- .bef_site_hpd(case$engine, .9, case$weights)[1L, ]
    expected <- .6 + sd * qnorm(c(.05, .95))
    expect_equal(unname(limits), expected, tolerance = 1e-12)
    expect_equal(unname(diff(normal_summary_cdf(limits, case))), .9, tolerance = 1e-10)
  }
})

test_that("the mode of a mixture of normals is the highest of its peaks, however far apart or alike they are", {
  cases <- list(
    normal_summary_case(c(-8, 0, 7), .01, c(.2, .3, .5)),
    normal_summary_case(c(-1, 0, 1), .45, c(.3, .15, .55)),
    normal_summary_case(c(-1.001, 1.001), 1, c(.5, .5)),
    normal_summary_case(c(0, 1000), .001, c(1 - 1e-9, 1e-9))
  )
  for (case in cases) {
    actual <- .bef_site_mode(case$engine, case$weights)
    expected <- normal_summary_mode_reference(case)
    actual_density <- sum(case$weight * dnorm(actual, case$mean, case$sd))
    expected_density <- sum(case$weight * dnorm(expected, case$mean, case$sd))
    expect_equal(actual_density, expected_density, tolerance = 1e-10)
    # Symmetric mixtures can have two equally high modes.
    expect_lt(abs(abs(actual) - abs(expected)), 1e-6 * case$sd)
  }
})

test_that("the search for the shortest interval compares several minima", {
  cases <- list(
    normal_summary_case(c(-5, 3), .3, c(.35, .65)),
    normal_summary_case(c(-2, -.5, 2), .4, c(.2, .3, .5)),
    normal_summary_case(c(-3, 3), .2, c(.5, .5))
  )
  levels <- c(.5, .6, .4)
  for (i in seq_along(cases)) {
    case <- cases[[i]]
    mixture <- .bef_normal_mixture_parts(case$mean, case$sd, case$weight)
    limits <- .bef_normal_mixture_hpd_exact(mixture, levels[i])
    reference <- normal_summary_hpd_reference(case, levels[i])
    expect_equal(diff(limits), diff(reference), tolerance = 1e-7 * case$sd)
    expect_equal(diff(normal_summary_cdf(limits, case)), levels[i], tolerance = 1e-9)
    endpoint_density <- vapply(limits, function(x) {
      sum(case$weight * dnorm(x, case$mean, case$sd))
    }, numeric(1L))
    expect_equal(endpoint_density[1L], endpoint_density[2L], tolerance = 1e-6)
    # What the accessors use: the level exactly, and a length within one
    # hundredth of the shortest.
    quick <- .bef_site_hpd(case$engine, levels[i], case$weights)[1L, ]
    expect_equal(unname(diff(normal_summary_cdf(quick, case))), levels[i], tolerance = 1e-9)
    expect_gte(unname(diff(quick)), diff(reference) * (1 - 1e-7))
    expect_lte(unname(diff(quick)), diff(reference) * 1.0102)
  }
})

test_that("the search for a quantile does not stop at one of its own values", {
  # In these mixtures the distribution function at one of the values that
  # bracket the quantiles is a probability the search is asked for, and the
  # end of the bracket is then the quantile itself.
  mixtures <- list(
    list(mean = c(-0.21, 0.08), weight = c(0.24, 0.95)),
    list(mean = c(-0.6, -0.31), weight = c(0.51, 0.53)),
    list(mean = c(-0.04, 0), weight = c(0.56, 0.75)),
    list(mean = c(-0.41, -0.37, -0.17), weight = c(0.14, 0.69, 0.42))
  )
  for (m in mixtures) {
    mixture <- .bef_normal_mixture_parts(m$mean, 1, m$weight)
    cdf <- function(q) sum(mixture$weight * pnorm(q - mixture$mean))
    inverse <- .bef_normal_mixture_quantiles(mixture)
    probabilities <- inverse$cdf[inverse$cdf > 1e-12 & inverse$cdf < 1 - 1e-12]
    quantiles <- vapply(probabilities, inverse$quantile, numeric(1L))
    expect_lt(max(abs(vapply(quantiles, cdf, numeric(1L)) - probabilities)), 1e-9)
    for (level in c(.4, .5)) {
      limits <- (.bef_normal_mixture_hpd_exact(mixture, level) - mixture$origin) / mixture$sd
      expect_equal(cdf(limits[2L]) - cdf(limits[1L]), level, tolerance = 1e-9)
    }
  }
})

test_that("the limits of the quick shortest interval are where the search puts them", {
  # Many intervals between values of a mesh have the same length. For a
  # symmetric mixture the shortest one is symmetric.
  cases <- list(
    list(mean = c(-1, 1), weight = c(.5, .5)),
    list(mean = c(-.5, 0, .5), weight = c(.25, .5, .25)),
    list(mean = seq(-3, 3, by = .5), weight = dnorm(seq(-3, 3, by = .5)))
  )
  for (case in cases) {
    mixture <- .bef_normal_mixture_parts(case$mean, 1, case$weight)
    for (level in c(.5, .9, .95)) {
      quick <- .bef_normal_mixture_hpd(mixture, level)
      expect_equal(quick, .bef_normal_mixture_hpd_exact(mixture, level), tolerance = 1e-6)
      expect_lt(abs(sum(quick)), 1e-6)
    }
  }
  # A shortest interval is not longer than the interval with equal tails.
  fit <- plot_fit(12L, component = "normal")
  for (level in c(.5, .9)) {
    shortest <- confint(fit, level = level, interval = "hpd")
    equal <- confint(fit, level = level)
    expect_true(all(shortest$upper - shortest$lower <= (equal$upper - equal$lower) * (1 + 1e-7)))
  }
})

test_that("the quick shortest interval holds the level and is within a hundredth of the shortest", {
  withr::local_seed(4127)
  ratios <- numeric()
  for (i in seq_len(24L)) {
    components <- sample(c(3L, 12L, 60L), 1L)
    # Units of the component SD: overlapping, apart, and far apart.
    spacing <- sample(c(0.05, 0.7, 4, 300), 1L)
    weight <- rexp(components)^sample(c(1, 4), 1L)
    mixture <- .bef_normal_mixture_parts(spacing * seq_len(components), 1, weight)
    cdf <- function(q) sum(mixture$weight * pnorm(q - mixture$mean))
    for (level in c(.2, .5, .9)) {
      quick <- (.bef_normal_mixture_hpd(mixture, level) - mixture$origin) / mixture$sd
      exact <- (.bef_normal_mixture_hpd_exact(mixture, level) - mixture$origin) / mixture$sd
      expect_equal(cdf(quick[2L]) - cdf(quick[1L]), level, tolerance = 1e-8)
      ratios <- c(ratios, diff(quick) / diff(exact))
    }
  }
  expect_gte(min(ratios), 1 - 1e-8)
  expect_lte(max(ratios), 1.0102)
})

test_that("a remote component of negligible weight does not make the shortest interval less exact", {
  case <- normal_summary_case(c(0, 1000), .001, c(1 - 1e-9, 1e-9))
  limits <- .bef_site_hpd(case$engine, .9, case$weights)[1L, ]
  reference <- .001 * qnorm(c((1 - .9 / case$weight[1L]) / 2,
                              (1 + .9 / case$weight[1L]) / 2))
  expect_lt(max(abs(unname(limits) - reference)), 1e-10)
  expect_equal(unname(diff(normal_summary_cdf(limits, case))), .9, tolerance = 1e-9)
})

test_that("the mode and the shortest interval of a mixture of normals follow a change of location and scale", {
  case <- normal_summary_case(c(-2, -.5, 2), .4, c(.2, .3, .5))
  mode <- .bef_site_mode(case$engine, case$weights)
  limits <- .bef_site_hpd(case$engine, .6, case$weights)[1L, ]
  for (scale in c(1e-3, 1e3)) {
    shift <- 10000
    transformed <- normal_summary_case(shift + scale * case$mean,
                                       scale * case$sd, case$weight)
    expect_lt(abs(.bef_site_mode(transformed$engine, transformed$weights) -
                    (shift + scale * mode)), 1e-8 * scale)
    # The shortest interval is found to within a hundredth of its length.
    moved <- .bef_site_hpd(transformed$engine, .6, transformed$weights)[1L, ]
    expect_lt(max(abs(moved - (shift + scale * limits))), 0.0102 * scale * diff(limits))
    expect_equal(unname(diff(normal_summary_cdf(moved, transformed))), .6, tolerance = 1e-9)
    exact <- .bef_normal_mixture_hpd_exact(
      .bef_normal_mixture_parts(transformed$mean, transformed$sd, transformed$weight), .6
    )
    base <- .bef_normal_mixture_hpd_exact(
      .bef_normal_mixture_parts(case$mean, case$sd, case$weight), .6
    )
    expect_lt(max(abs(exact - (shift + scale * base))), 1e-7 * scale)
  }
})

test_that("for a posterior that is one narrow normal component, coef() and confint() give its mode and shortest interval", {
  context <- .bef_fit_prepare_context(
    call = quote(bayes_efron_fit(y, se)), theta_hat = rep(-9.8, 5L), sigma = rep(.2, 5L),
    support = c(-10, 10), component = "normal", width = .005, L = 101L,
    chains = 2L, iter_sampling = 200L, seed = 130L, check_installed = FALSE
  )
  data <- context$stan_data
  alpha <- coef(lm.fit(cbind(1, data$B), dnorm(data$grid, -10, .01, log = TRUE)))[-1L]
  context$cmdstan_fit <- plot_fit_cmdstan(matrix(rep(alpha, 400L), nrow = 400L, byrow = TRUE), data)
  context$runtime_seconds <- 0
  context$sampler_settings <- .bef_fit_sampler_settings(context, interactive_fun = function() FALSE)
  fit <- suppressWarnings(.bef_fit_assemble(context, cmdstan_version_fun = function() "2.38.0"))
  # Only the leftmost component has appreciable probability in this fit.
  component <- .bef_engine_components(.bef_engine(fit), seq_len(5L))
  expected_mode <- component$mean[, 1L]
  expected_interval <- outer(expected_mode, c(1, 1)) +
    outer(component$sd, qnorm(c(.05, .95)))
  expect_lt(max(abs(coef(fit, type = "mode") - expected_mode)), 1e-10)
  limits <- as.matrix(confint(fit, interval = "hpd")[, c("lower", "upper")])
  expect_lt(max(abs(limits - expected_interval)), 1e-10)
  expect_equal(site_summary(fit)$mode, unname(coef(fit, type = "mode")))
})

test_that("with normal components the shortest intervals are a matrix with a row per site, in the order asked for", {
  engine <- .bef_engine_from_parts(matrix(log(c(.2, .3, .5)), nrow = 1L),
                                   grid = c(-2, 0, 2), theta_hat = c(-1, 0, 1),
                                   sigma = c(.2, .3, .4), tau_c = .05)
  all <- .bef_site_hpd(engine, .8)
  selected <- .bef_site_hpd(engine, .8, sites = c(3, 1))
  expect_identical(dim(all), c(3L, 2L))
  expect_identical(colnames(all), c("lower", "upper"))
  expect_equal(selected, all[c(3, 1), , drop = FALSE])
})

test_that("point modes and shortest intervals follow the two readings", {
  engine <- .bef_engine_from_parts(matrix(log(c(.1, .4, .4, .1)), nrow = 1L),
                                   grid = seq(-1.5, 1.5, by = 1),
                                   theta_hat = 0, sigma = 1)
  weights <- matrix(c(.1, .4, .4, .1), nrow = 1L)
  expect_equal(.bef_site_mode(engine, weights), -.5)
  # On the points: the two middle points hold 0.8.
  expect_equal(unname(.bef_site_hpd(engine, .6, weights, reading = "step")[1L, ]), c(-.5, .5))
  # Over the cells: any interval of length 1.5 inside the two middle cells
  # holds 0.6, and the middle one is given.
  expect_equal(unname(.bef_site_hpd(engine, .6, weights, reading = "interpolated")[1L, ]),
               c(-.75, .75))
})
