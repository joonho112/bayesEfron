# The posterior distribution of an effect whose prior distribution is
# uniform within each cell of a grid, given an estimate with a normal error:
# the weight of a cell is its prior probability times the likelihood
# averaged over the cell, and within a cell the posterior distribution is a
# normal distribution cut at the ends of the cell. few-studies.R sources this
# file for the rows called "known g" of its third distribution, whose effects
# are drawn in that way. histogram-posterior-check.R sets the formulas beside
# numerical integration.
normal_interval_logprob <- function(a, b) {
  right <- a >= 0
  top <- stats::pnorm(b, log.p = TRUE)
  bottom <- stats::pnorm(a, log.p = TRUE)
  top[right] <- stats::pnorm(a[right], lower.tail = FALSE, log.p = TRUE)
  bottom[right] <- stats::pnorm(b[right], lower.tail = FALSE, log.p = TRUE)
  top + log(-expm1(pmin(bottom - top, 0)))
}

histogram_site <- function(y, sigma, lower, upper, mass) {
  active <- mass > 0
  lower <- lower[active]
  upper <- upper[active]
  mass <- mass[active]
  a <- (lower - y) / sigma
  b <- (upper - y) / sigma
  log_probability <- normal_interval_logprob(a, b)
  log_weight <- log(mass) - log(upper - lower) + log_probability
  weight <- exp(log_weight - max(log_weight))
  weight <- weight / sum(weight)
  # First moments of the truncated normals, computed from log densities
  # so a cell in a far tail need not have a representable likelihood.
  log_a <- stats::dnorm(a, log = TRUE)
  log_b <- stats::dnorm(b, log = TRUE)
  top <- pmax(log_a, log_b)
  bottom <- pmin(log_a, log_b)
  ratio <- sign(log_a - log_b) * exp(top + log(-expm1(bottom - top)) - log_probability)
  ratio[log_a == log_b] <- 0
  cell_mean <- pmin(pmax(y + sigma * ratio, lower), upper)
  mean <- sum(weight[weight > 0] * cell_mean[weight > 0])
  cumulative <- cumsum(weight)
  cumulative[length(cumulative)] <- 1
  cdf <- function(q) {
    vapply(q, function(value) {
      z <- pmax(a, pmin((value - y) / sigma, b))
      within <- exp(normal_interval_logprob(a, z) - log_probability)
      within[value <= lower] <- 0
      within[value >= upper] <- 1
      sum(weight[weight > 0] * within[weight > 0])
    }, numeric(1L))
  }
  quantile <- function(p) {
    vapply(p, function(probability) {
      if (probability <= 0) return(lower[1L])
      if (probability >= 1) return(upper[length(upper)])
      j <- which(cumulative >= probability)[1L]
      before <- if (j == 1L) 0 else cumulative[j - 1L]
      fraction <- pmin(pmax((probability - before) / weight[j], 0), 1)
      if (fraction <= 0) return(lower[j])
      if (fraction >= 1) return(upper[j])
      if (a[j] >= 0) {
        log_tail <- stats::pnorm(a[j], lower.tail = FALSE, log.p = TRUE)
        log_tail <- log_tail + log1p(-fraction * exp(log_probability[j] - log_tail))
        z <- stats::qnorm(log_tail, lower.tail = FALSE, log.p = TRUE)
      } else {
        log_left <- stats::pnorm(a[j], log.p = TRUE)
        log_piece <- log(fraction) + log_probability[j]
        top <- max(log_left, log_piece)
        log_cdf <- top + log(exp(log_left - top) + exp(log_piece - top))
        z <- stats::qnorm(log_cdf, log.p = TRUE)
      }
      pmin(pmax(y + sigma * z, lower[j]), upper[j])
    }, numeric(1L))
  }
  list(mean = mean, cdf = cdf, quantile = quantile, weight = weight)
}

histogram_sites <- function(y, sigma, grid, mass, level = .9) {
  spacing <- grid[2L] - grid[1L]
  lower <- grid - spacing / 2
  upper <- grid + spacing / 2
  probs <- c((1 - level) / 2, (1 + level) / 2)
  rows <- lapply(seq_along(y), function(i) {
    posterior <- histogram_site(y[i], sigma[i], lower, upper, mass)
    limits <- posterior$quantile(probs)
    data.frame(mean = posterior$mean, lower = limits[1L], upper = limits[2L])
  })
  do.call(rbind, rows)
}
