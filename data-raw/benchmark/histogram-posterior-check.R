# Checks of histogram-posterior.R, from the package root:
#
#   Rscript data-raw/benchmark/histogram-posterior-check.R
#
# The formulas are set beside numerical integration, for cells of unequal
# width, for estimates far in a tail and for cells without probability.
# Then the 200 known-g rows of the third distribution of few-studies.R are
# computed again from their seeds and compared with the stored rows and
# with their row of inst/verification/few-studies.csv. The script stops at
# the first check that fails.

suppressMessages(devtools::load_all(".", quiet = TRUE))
source(file.path("data-raw", "benchmark", "histogram-posterior.R"))

# The same posterior by numerical integration. The likelihood is scaled
# before it is integrated, so that an estimate fifty standard errors away
# is not lost.
quadrature <- function(y, sigma, lower, upper, mass) {
  closest <- pmin(pmax(y, lower), upper)
  top <- max(stats::dnorm(y, closest, sigma, log = TRUE))
  piece <- function(j, end = upper[j], moment = FALSE) {
    if (end <= lower[j]) return(0)
    stats::integrate(function(x) {
      (if (moment) x else 1) * mass[j] / (upper[j] - lower[j]) *
        exp(stats::dnorm(y, x, sigma, log = TRUE) - top)
    }, lower[j], end, rel.tol = 1e-12, abs.tol = 1e-12)$value
  }
  cells <- vapply(seq_along(mass), piece, numeric(1L))
  total <- sum(cells)
  mean <- sum(vapply(seq_along(mass), function(j) piece(j, moment = TRUE), numeric(1L))) / total
  cdf <- function(q) {
    vapply(q, function(value) {
      sum(vapply(seq_along(mass), function(j) piece(j, min(value, upper[j])), numeric(1L))) / total
    }, numeric(1L))
  }
  quantile <- function(p) {
    vapply(p, function(probability) {
      stats::uniroot(function(q) cdf(q) - probability, range(c(lower, upper)), tol = 1e-12)$root
    }, numeric(1L))
  }
  list(mean = mean, cdf = cdf, quantile = quantile, weight = cells / total)
}

# Cells of unequal width: the prior density of a cell is its probability
# divided by its width.
lower <- c(-2, -.5, 1)
upper <- c(-.5, 1, 3)
mass <- c(.6, .1, .3)
probs <- c(.05, .5, .95)
for (y in c(-1.3, .4, 2.7)) {
  exact <- histogram_site(y, .6, lower, upper, mass)
  numeric <- quadrature(y, .6, lower, upper, mass)
  at <- c(-2, -1, 0, 2, 3)
  stopifnot(
    abs(exact$mean - numeric$mean) < 1e-10,
    max(abs(exact$weight - numeric$weight)) < 1e-10,
    max(abs(exact$cdf(at) - numeric$cdf(at))) < 1e-10,
    max(abs(exact$quantile(probs) - numeric$quantile(probs))) < 1e-9
  )
}

# Estimates far in a tail.
for (y in c(20, 50, -20, -50)) {
  exact <- histogram_site(y, 1, c(-1, 0), c(0, 1), c(.2, .8))
  numeric <- quadrature(y, 1, c(-1, 0), c(0, 1), c(.2, .8))
  stopifnot(
    abs(exact$mean - numeric$mean) < 1e-9,
    max(abs(exact$quantile(probs) - numeric$quantile(probs))) < 1e-9,
    max(abs(exact$cdf(exact$quantile(probs)) - probs)) < 1e-10,
    abs(sum(exact$weight) - 1) < 1e-12
  )
}

# A cell without probability, and the probabilities 0 and 1.
exact <- histogram_site(.2, .3, c(-2, -1, 0), c(-1, 0, 1), c(0, .25, .75))
stopifnot(
  identical(exact$quantile(c(0, 1)), c(-1, 1)),
  isTRUE(all.equal(exact$cdf(c(-2, -1, 1, 2)), c(0, 0, 1, 1))),
  max(abs(exact$cdf(exact$quantile(c(.05, .95))) - c(.05, .95))) < 1e-10
)

# The known-g rows of the third distribution, from their seeds as in
# few-studies.R.
g <- g_summary(raudenbush_fit, type = "mass")
grid <- g$value
mass <- g$mean / sum(g$mean)
spacing <- grid[2L] - grid[1L]
sigma <- raudenbush1985$sei
K <- length(sigma)
replications <- 200L
computed <- do.call(rbind, lapply(seq_len(replications), function(r) {
  set.seed(20261200L + 3000L + r)
  theta <- sample(grid, K, replace = TRUE, prob = mass) + stats::runif(K, -spacing / 2, spacing / 2)
  y <- stats::rnorm(K, theta, sigma)
  sites <- histogram_sites(y, sigma, grid, mass)
  inside <- theta >= sites$lower & theta <= sites$upper
  largest <- which.max(y)
  data.frame(replication = r, inclusion = mean(inside), length = mean(sites$upper - sites$lower),
             rmse = sqrt(mean((sites$mean - theta)^2)), largest_inside = inside[largest],
             largest_error = sites$mean[largest] - theta[largest])
}))
rows <- utils::read.csv(file.path("data-raw", "benchmark", "results-v040", "few-studies-rows.csv"),
                        stringsAsFactors = FALSE)
rows <- rows[rows$scenario == "fitted g" & rows$method == "known g", ]
rows <- rows[order(rows$replication), ]
measures <- c("inclusion", "length", "rmse", "largest_error")
stopifnot(
  identical(rows$replication, computed$replication),
  max(abs(as.matrix(rows[, measures]) - as.matrix(computed[, measures]))) < 1e-12,
  identical(as.logical(rows$largest_inside), computed$largest_inside)
)
se <- function(x) stats::sd(x) / sqrt(length(x))
expected <- c(
  inclusion = mean(computed$inclusion), inclusion_se = se(computed$inclusion),
  length = mean(computed$length), length_se = se(computed$length),
  rmse = mean(computed$rmse), rmse_se = se(computed$rmse),
  largest_inclusion = mean(computed$largest_inside),
  largest_error = mean(computed$largest_error), largest_error_se = se(computed$largest_error)
)
table <- utils::read.csv(file.path("inst", "verification", "few-studies.csv"), stringsAsFactors = FALSE)
row <- table[table$scenario == "fitted g" & table$method == "known g", ]
stopifnot(nrow(row) == 1L, max(abs(unlist(row[, names(expected)]) - expected)) < 1e-12)

cat("histogram-posterior.R: the formulas agree with numerical integration, and the",
    replications, "known-g rows and their row of the table follow from their seeds.\n")
