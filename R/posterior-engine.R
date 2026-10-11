# The posterior distribution of the site effects, computed from the draws of g.
#
# g has one component at each grid point: a normal distribution with standard
# deviation tau_c, or a point mass when tau_c is 0. Its probabilities are
# exp(log_g), one row of `log_g` for each posterior draw. Given g, site i
# belongs to component j with probability proportional to
#
#   g[j] * normal(theta_hat[i] | grid[j], sqrt(sigma[i]^2 + tau_c^2)),
#
# and within component j the effect theta[i] is normal with mean
# grid[j] + b[i] * (theta_hat[i] - grid[j]) and variance b[i] * sigma[i]^2,
# where b[i] = tau_c^2 / (tau_c^2 + sigma[i]^2). With tau_c = 0 the effect is
# the grid point itself. The posterior distribution of theta[i] averages
# this over the draws of g (Lee and Sui 2025, Section 4.3), so it is a
# mixture over the grid points with the averaged probabilities as weights.
#
# Everything a method reports about a site (mean, standard deviation,
# quantiles, probabilities, log likelihood) is computed here from log_g and
# the data, and not from site-level draws.

# What the functions below need from a fit, in one list. Fits saved by
# versions 0.1 to 0.3 are accepted: they hold log_g (0.1 holds g).
.bef_engine <- function(fit) {
  data <- fit$metadata$data_list
  L <- as.integer(data$L)
  K <- as.integer(data$K)
  grid <- as.numeric(data$grid)
  tau_c <- if (is.null(data$tau_c)) 0 else as.numeric(data$tau_c)
  list(
    grid = grid,
    L = L,
    spacing = if (L > 1L) (grid[L] - grid[1L]) / (L - 1L) else 0,
    tau_c = tau_c,
    log_g = .bef_engine_log_g(fit$draws, L),
    theta_hat = as.numeric(data$theta_hat),
    sigma = as.numeric(data$sigma),
    K = K
  )
}

# log_g as a plain matrix with one row per draw and one column per grid
# point.
.bef_engine_log_g <- function(draws, L) {
  dm <- .bef_draws_matrix(draws)
  log_names <- paste0("log_g[", seq_len(L), "]")
  if (all(log_names %in% colnames(dm))) {
    return(.bef_plain_draw_matrix(dm[, log_names, drop = FALSE]))
  }
  g_names <- paste0("g[", seq_len(L), "]")
  if (!all(g_names %in% colnames(dm))) {
    .bef_abort_invalid_fit(
      "This fit holds neither `log_g` nor `g`, so its site summaries cannot be computed."
    )
  }
  log(.bef_plain_draw_matrix(dm[, g_names, drop = FALSE]))
}

# An engine from its parts, for computations that have no fit object.
.bef_engine_from_parts <- function(log_g, grid, theta_hat, sigma, tau_c = 0) {
  L <- length(grid)
  list(
    grid = as.numeric(grid),
    L = L,
    spacing = if (L > 1L) (grid[L] - grid[1L]) / (L - 1L) else 0,
    tau_c = tau_c,
    log_g = log_g,
    theta_hat = as.numeric(theta_hat),
    sigma = as.numeric(sigma),
    K = length(theta_hat)
  )
}

# The likelihood of each estimate at each component, with the largest value
# of each site divided out: `lik` has one row per site and one column per
# grid point, and `row_max` holds the logarithm of what was divided out.
.bef_engine_likelihood <- function(engine, sites) {
  K <- length(sites)
  L <- engine$L
  s_tot <- sqrt(engine$sigma[sites]^2 + engine$tau_c^2)
  log_lik <- matrix(
    stats::dnorm(
      rep(engine$theta_hat[sites], times = L),
      mean = rep(engine$grid, each = K),
      sd = rep(s_tot, times = L),
      log = TRUE
    ),
    nrow = K,
    ncol = L
  )
  row_max <- apply(log_lik, 1L, max)
  list(lik = exp(log_lik - row_max), row_max = row_max, log_lik = log_lik)
}

# The sites in blocks, so that no intermediate matrix has more than
# `max_cells` elements.
.bef_engine_blocks <- function(sites, n_draws, max_cells = 2e6) {
  per_block <- max(1L, as.integer(max_cells %/% n_draws))
  split(sites, ceiling(seq_along(sites) / per_block))
}

.bef_engine_sites <- function(engine, sites) {
  if (is.null(sites)) seq_len(engine$K) else as.integer(sites)
}

# P(site i belongs to component j | data), averaged over the draws: one row
# per site and one column per grid point. The rows sum to 1.
.bef_site_weights <- function(engine, sites = NULL, max_cells = 2e6) {
  sites <- .bef_engine_sites(engine, sites)
  g <- exp(engine$log_g)
  n_draws <- nrow(g)
  out <- matrix(NA_real_, nrow = length(sites), ncol = engine$L)
  at <- 0L
  for (block in .bef_engine_blocks(sites, n_draws, max_cells)) {
    rows <- at + seq_along(block)
    at <- at + length(block)
    like <- .bef_engine_likelihood(engine, block)
    # marginal[s, i] = sum over j of g[s, j] * lik[i, j]
    marginal <- g %*% t(like$lik)
    # mean over s of g[s, j] * lik[i, j] / marginal[s, i]
    out[rows, ] <- like$lik * t(crossprod(g, 1 / marginal)) / n_draws
    # If a product underflowed for some draw, redo that site on the log scale.
    for (k in which(rowSums(!is.finite(out[rows, , drop = FALSE])) > 0L)) {
      weights <- .bef_engine_site_weights_log(engine$log_g, like$log_lik[k, ])
      out[rows[k], ] <- colMeans(weights)
    }
  }
  out
}

# The weights of one site for every draw, computed on the log scale: one row
# per draw and one column per grid point.
.bef_engine_site_weights_log <- function(log_g, log_lik) {
  log_post <- sweep(log_g, 2L, log_lik, "+")
  weights <- exp(log_post - apply(log_post, 1L, max))
  weights / rowSums(weights)
}

# The mean and the variance of theta[i] given g, for every draw: two
# matrices with one row per draw and one column per site.
.bef_site_conditional <- function(engine, sites = NULL, max_cells = 2e6) {
  sites <- .bef_engine_sites(engine, sites)
  g <- exp(engine$log_g)
  n_draws <- nrow(g)
  grid <- engine$grid
  mean <- matrix(NA_real_, nrow = n_draws, ncol = length(sites))
  var <- mean
  # Moments are taken about the middle of the grid and not about zero. A
  # variance formed as a difference of raw moments loses its digits when the
  # effects are far from zero; about the middle of the grid the loss is at
  # most the squared ratio of the half-width of the grid to the standard
  # deviation, times the machine precision.
  center <- (grid[1L] + grid[length(grid)]) / 2
  d <- grid - center
  at <- 0L
  for (block in .bef_engine_blocks(sites, n_draws, max_cells)) {
    cols <- at + seq_along(block)
    at <- at + length(block)
    like <- .bef_engine_likelihood(engine, block)
    marginal <- g %*% t(like$lik)
    m1 <- (g %*% t(sweep(like$lik, 2L, d, "*"))) / marginal
    m2 <- (g %*% t(sweep(like$lik, 2L, d^2, "*"))) / marginal
    mean[, cols] <- center + m1
    var[, cols] <- pmax(m2 - m1^2, 0)
    # If a product underflowed for some draw, redo that site on the log scale.
    bad <- which(colSums(!is.finite(m1) | !is.finite(m2)) > 0L)
    for (k in bad) {
      weights <- .bef_engine_site_weights_log(engine$log_g, like$log_lik[k, ])
      a1 <- as.numeric(weights %*% d)
      mean[, cols[k]] <- center + a1
      var[, cols[k]] <- pmax(as.numeric(weights %*% d^2) - a1^2, 0)
    }
  }
  # With normal components the effect is shrunk toward the estimate within
  # each component.
  if (engine$tau_c > 0) {
    b <- engine$tau_c^2 / (engine$tau_c^2 + engine$sigma[sites]^2)
    mean <- sweep(sweep(mean, 2L, 1 - b, "*"), 2L, b * engine$theta_hat[sites], "+")
    var <- sweep(sweep(var, 2L, (1 - b)^2, "*"), 2L, b * engine$sigma[sites]^2, "+")
  }
  list(mean = mean, var = var)
}

# The effective number of parameters for every draw: the sum over sites of
# the variance of theta[i] given g, divided by sigma[i]^2. The program
# efron_re.stan computes it in Stan; for efron_mix.stan it is computed here.
.bef_effective_params <- function(engine, max_cells = 2e6) {
  n_draws <- nrow(engine$log_g)
  out <- numeric(n_draws)
  for (block in .bef_engine_blocks(seq_len(engine$K), n_draws, max_cells)) {
    conditional <- .bef_site_conditional(engine, sites = block, max_cells = max_cells)
    out <- out + as.numeric(conditional$var %*% (1 / engine$sigma[block]^2))
  }
  out
}

# The site-level draws that version 0.3.0 stored, computed from the draws of
# g: for every draw and site the mean, the standard deviation and the mode of
# theta[i] given g, one draw of theta[i], and the log likelihood. Each
# element has one row per draw and one column per site. `seed` makes the
# draws of theta[i] reproducible; the state of R's generator is restored.
.bef_site_blocks <- function(engine, seed = NULL) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      get(".Random.seed", envir = globalenv(), inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        suppressWarnings(rm(".Random.seed", envir = globalenv(), inherits = FALSE))
      } else {
        assign(".Random.seed", old_seed, envir = globalenv())
      }
    }, add = TRUE)
    set.seed(seed)
  }
  K <- engine$K
  L <- engine$L
  n_draws <- nrow(engine$log_g)
  conditional <- .bef_site_conditional(engine)
  like <- .bef_engine_likelihood(engine, seq_len(K))
  shrink <- if (engine$tau_c > 0) {
    engine$tau_c^2 / (engine$tau_c^2 + engine$sigma^2)
  } else {
    rep(0, K)
  }
  theta_map <- matrix(NA_real_, nrow = n_draws, ncol = K)
  theta_rep <- theta_map
  for (i in seq_len(K)) {
    log_post <- sweep(engine$log_g, 2L, like$log_lik[i, ], "+")
    # The most probable component, the first one when several are equal, as
    # in efron_re.stan.
    top <- max.col(log_post, ties.method = "first")
    weights <- exp(log_post - log_post[cbind(seq_len(n_draws), top)])
    cumulative <- weights
    for (j in seq_len(L)[-1L]) {
      cumulative[, j] <- cumulative[, j - 1L] + weights[, j]
    }
    u <- stats::runif(n_draws) * cumulative[, L]
    index <- pmin(rowSums(cumulative < u) + 1L, L)
    center_map <- engine$grid[top]
    center_rep <- engine$grid[index]
    if (shrink[i] > 0) {
      theta_map[, i] <- center_map + shrink[i] * (engine$theta_hat[i] - center_map)
      theta_rep[, i] <- stats::rnorm(
        n_draws, center_rep + shrink[i] * (engine$theta_hat[i] - center_rep),
        sqrt(shrink[i]) * engine$sigma[i]
      )
    } else {
      theta_map[, i] <- center_map
      theta_rep[, i] <- center_rep
    }
  }
  list(
    theta_map = theta_map,
    theta_mean = conditional$mean,
    theta_sd = sqrt(conditional$var),
    theta_rep = theta_rep,
    log_lik = .bef_log_lik(engine)
  )
}

# The posterior mean and variance of each site effect. The variance has two
# parts, returned as well: `var_within`, the variance of theta[i] given g
# averaged over the draws, and `var_between`, the variance over the draws of
# the mean of theta[i] given g.
.bef_site_moments <- function(engine, sites = NULL, max_cells = 2e6) {
  conditional <- .bef_site_conditional(engine, sites = sites, max_cells = max_cells)
  mean <- colMeans(conditional$mean)
  within <- colMeans(conditional$var)
  between <- colMeans(sweep(conditional$mean, 2L, mean)^2)
  data.frame(
    mean = mean,
    var = within + between,
    var_within = within,
    var_between = between
  )
}

# log p(theta_hat[i] | g) for every draw: one row per draw and one column
# per site, the shape that the loo package expects.
.bef_log_lik <- function(engine, sites = NULL, max_cells = 2e6) {
  sites <- .bef_engine_sites(engine, sites)
  g <- exp(engine$log_g)
  n_draws <- nrow(g)
  out <- matrix(NA_real_, nrow = n_draws, ncol = length(sites))
  at <- 0L
  for (block in .bef_engine_blocks(sites, n_draws, max_cells)) {
    cols <- at + seq_along(block)
    at <- at + length(block)
    like <- .bef_engine_likelihood(engine, block)
    out[, cols] <- sweep(log(g %*% t(like$lik)), 2L, like$row_max, "+")
    for (k in which(colSums(!is.finite(out[, cols, drop = FALSE])) > 0L)) {
      log_post <- sweep(engine$log_g, 2L, like$log_lik[k, ], "+")
      top <- apply(log_post, 1L, max)
      out[, cols[k]] <- top + log(rowSums(exp(log_post - top)))
    }
  }
  out
}

# The mean and the standard deviation of theta[i] within component j, for
# normal components: `mean` has one row per site and one column per grid
# point, `sd` one value per site.
.bef_engine_components <- function(engine, sites) {
  b <- engine$tau_c^2 / (engine$tau_c^2 + engine$sigma[sites]^2)
  list(
    mean = outer(1 - b, engine$grid) + b * engine$theta_hat[sites],
    sd = sqrt(b) * engine$sigma[sites]
  )
}

# Work in units of the common component SD. Centering on an active mean
# makes accuracy depend on the components, rather than the support range.
.bef_normal_mixture_parts <- function(mean, sd, weight) {
  active <- weight > 0
  mean <- mean[active]
  weight <- weight[active]
  origin <- mean[which.max(weight)]
  list(mean = (mean - origin) / sd, sd = sd, origin = origin,
       weight = weight / sum(weight))
}

# Component neighborhoods bracket quantiles even when a remote component
# makes the full range enormous relative to the width of the posterior.
.bef_normal_mixture_quantiles <- function(mixture) {
  mu <- mixture$mean
  w <- mixture$weight
  x <- sort(unique(as.vector(outer(mu, c(-8, -5, -3, -2, -1, 0, 1, 2, 3, 5, 8), `+`))))
  cdf <- vapply(x, function(value) sum(w * stats::pnorm(value - mu)), numeric(1L))
  survival <- vapply(x, function(value) {
    sum(w * stats::pnorm(value - mu, lower.tail = FALSE))
  }, numeric(1L))
  quantile <- function(p) {
    if (p <= 0) return(-Inf)
    if (p >= 1) return(Inf)
    if (length(mu) == 1L) return(mu + stats::qnorm(p))
    if (p <= 0.5) {
      below <- cdf < p
      above <- cdf >= p
      fun <- function(value) sum(w * stats::pnorm(value - mu)) - p
    } else {
      # Use the survival function for upper quantiles to retain small tails.
      below <- survival > 1 - p
      above <- survival <= 1 - p
      fun <- function(value) 1 - p - sum(w * stats::pnorm(value - mu, lower.tail = FALSE))
    }
    a <- if (any(below)) max(x[below]) else min(mu) + stats::qnorm(p)
    b <- if (any(above)) min(x[above]) else max(mu) + stats::qnorm(p)
    if (a == b) return(a)
    # An end of the bracket can be the quantile itself, because the starts of
    # the search include the distribution function at these values. The
    # function is then zero there, and its value at a + (b - a), which need
    # not be b, can have the other sign. So the ends are evaluated here and
    # handed to uniroot().
    f_a <- fun(a)
    f_b <- fun(b)
    if (f_a >= 0) return(a)
    if (f_b <= 0) return(b)
    # Solve in a local coordinate so uniroot's relative tolerance does not
    # grow with the distance between clusters of components.
    a + stats::uniroot(function(delta) fun(a + delta), c(0, b - a),
                       f.lower = f_a, f.upper = f_b,
                       tol = 1e-10, check.conv = TRUE)$root
  }
  list(quantile = quantile, cdf = cdf)
}

# The shortest interval of a normal mixture with probability `level`. The
# distribution function is first evaluated on equally spaced values, which
# is quick. A shortest interval between two of those values is at most two
# steps longer than the shortest interval, so the mesh is taken where that
# interval spans at least `steps` steps: the length is then within
# 2 / (steps - 2) of the shortest, one hundredth for 202 steps. Where it
# spans fewer, a finer mesh is tried, and where that would be too large, the
# search of .bef_normal_mixture_hpd_exact().
.bef_normal_mixture_hpd <- function(mixture, level, steps = 202, n = 2001L,
                                    max_cells = 4e6) {
  if (length(mixture$mean) == 1L) {
    return(mixture$origin + mixture$sd * stats::qnorm(c((1 - level) / 2,
                                                       (1 + level) / 2)))
  }
  # The mesh leaves out a probability of about 1e-10 at each side.
  if (level <= 1 - 1e-8) {
    mesh <- .bef_normal_mixture_hpd_mesh(mixture, level, n)
    if (mesh$steps < steps) {
      # The shortest interval is at most two steps shorter than the one
      # found between values of the mesh.
      shortest <- (mesh$steps - 2) * mesh$step
      finer <- if (shortest > 0) ceiling(1.25 * steps * mesh$span / shortest) + 1 else Inf
      if (finer * mesh$components <= max_cells) {
        mesh <- .bef_normal_mixture_hpd_mesh(mixture, level, as.integer(finer))
      }
    }
    if (mesh$steps >= steps) {
      return(mixture$origin + mixture$sd * mesh$limits)
    }
  }
  .bef_normal_mixture_hpd_exact(mixture, level)
}

# The shortest interval read from the distribution function on `n` equally
# spaced values, in units of the component SD. The means are increasing, as
# those of a site are. The values cover the components that hold all but
# 1e-10 of the probability at each side, and 6.5 SDs beyond them.
#
# `steps` is the number of steps that the shortest interval between two of
# the values spans; that interval holds `level` or a little more. The
# limits that are returned hold `level`. They are found in two ways, and the
# shorter interval is kept. One moves an end of that interval in, by less
# than a step. The other starts from the value at which an interval with
# probability `level` is shortest when its upper limit is read between the
# values, and moves its lower limit until no interval nearby is shorter. The
# first is never more than two steps longer than the shortest interval; the
# second is the one that puts the limits where they belong, since many
# intervals between values of the mesh have the same length.
.bef_normal_mixture_hpd_mesh <- function(mixture, level, n) {
  w <- mixture$weight
  keep <- which(cumsum(w) > 1e-10 & rev(cumsum(rev(w))) > 1e-10)
  mu <- mixture$mean[keep]
  w <- w[keep] / sum(w[keep])
  x <- seq(mu[1L] - 6.5, mu[length(mu)] + 6.5, length.out = n)
  # The sum can come out lower than at the value before by one unit in the
  # last place, depending on how the platform adds; the search below needs a
  # distribution function that does not decrease.
  cdf <- cummax(as.vector(stats::pnorm(outer(x, mu, `-`)) %*% w))
  limits <- .bef_shortest_interval(x, cdf, level)
  step <- x[2L] - x[1L]
  found <- list(step = step, steps = diff(limits) / step, span = x[n] - x[1L],
                components = length(mu))
  i <- match(limits[1L], x)
  j <- match(limits[2L], x)
  if (j - i < 2L) {
    return(c(found, list(limits = limits)))
  }
  cdf_at <- function(q) sum(w * stats::pnorm(q - mu))
  # Where the distribution function takes the value `p` between x[k] and
  # x[k + 1]. The values at the two ends are computed as the search
  # computes them, so that their signs are those it will see.
  crossing <- function(p, k) {
    f_a <- cdf_at(x[k]) - p
    f_b <- cdf_at(x[k + 1L]) - p
    if (f_a >= 0) return(x[k])
    if (f_b <= 0) return(x[k + 1L])
    x[k] + stats::uniroot(function(delta) cdf_at(x[k] + delta) - p, c(0, step),
                          f.lower = f_a, f.upper = f_b, tol = 1e-10)$root
  }
  upper_in <- c(limits[1L], crossing(cdf[i] + level, j - 1L))
  lower_in <- c(crossing(cdf[j] - level, i), limits[2L])
  limits <- if (diff(upper_in) <= diff(lower_in)) upper_in else lower_in

  # For every lower limit on the mesh, the upper limit by interpolation.
  target <- cdf + level
  last <- findInterval(target, cdf, left.open = TRUE) + 1L
  start <- which(last <= n)
  below <- last[start] - 1L
  upper <- x[below] + step * (target[start] - cdf[below]) / (cdf[last[start]] - cdf[below])
  best <- start[which.min(upper - x[start])]
  # The lower limit moves between the neighbors of that value, and not past
  # the last value from which an interval with probability `level` fits on
  # the mesh.
  quantile_at <- function(p) {
    k <- findInterval(p, cdf)
    if (k >= n) return(Inf)
    crossing(p, max(k, 1L))
  }
  length_from <- function(lower) quantile_at(cdf_at(lower) + level) - lower
  from <- x[max(1L, best - 1L)]
  to <- x[min(max(start), best + 1L)]
  if (from < to) {
    lower <- stats::optimize(length_from, c(from, to), tol = 1e-9)$minimum
    refined <- c(lower, quantile_at(cdf_at(lower) + level))
    holds <- is.finite(refined[2L]) &&
      cdf_at(refined[2L]) - cdf_at(refined[1L]) >= level - 1e-9
    if (holds && diff(refined) < diff(limits)) {
      limits <- refined
    }
  }
  c(found, list(limits = limits))
}

# A shortest connected interval has ends Q(u), Q(u + level). Scan starts
# at both a regular probability mesh and the component neighborhoods,
# refine each candidate minimum, and compare their lengths globally.
.bef_normal_mixture_hpd_exact <- function(mixture, level) {
  inverse <- .bef_normal_mixture_quantiles(mixture)
  q <- inverse$quantile
  u <- sort(unique(c(seq(0, 1 - level, length.out = 129L),
                     inverse$cdf, inverse$cdf - level)))
  u <- u[u >= 0 & u <= 1 - level]
  interval <- function(start) c(q(start), q(start + level))
  width <- function(start) diff(interval(start))
  widths <- vapply(u, width, numeric(1L))
  candidates <- u[which.min(widths)]
  minima <- which(widths <= c(Inf, utils::head(widths, -1L)) &
                    widths <= c(utils::tail(widths, -1L), Inf) & is.finite(widths))
  for (i in minima) {
    a <- u[max(1L, i - 1L)]
    b <- u[min(length(u), i + 1L)]
    if (a == b) next
    result <- stats::optimize(width, c(a, b), tol = 1e-10 * (1 - level))
    candidates <- c(candidates, result$minimum)
  }
  best <- candidates[which.min(vapply(candidates, width, numeric(1L)))]
  mixture$origin + mixture$sd * interval(best)
}

# The mode of a normal mixture. A maximum lies between the component means
# and within one SD of one of them: further out every component is convex.
# The density and its slope are computed at values around every mean, the
# places where the slope changes sign and the largest values are refined,
# and the candidates are compared.
.bef_normal_mixture_mode <- function(mixture, n = 1001L) {
  mu <- mixture$mean
  w <- mixture$weight
  if (length(mu) == 1L) return(mixture$origin)
  density <- function(x) sum(w * stats::dnorm(x - mu))
  slope <- function(x) sum(w * (mu - x) * stats::dnorm(x - mu))
  x <- sort(unique(c(as.vector(outer(mu, seq(-1, 1, length.out = 17L), `+`)),
                     seq(min(mu), max(mu), length.out = n))))
  x <- x[x >= min(mu) & x <= max(mu)]
  # One pass over the components for all values of x.
  d <- numeric(length(x))
  s <- numeric(length(x))
  for (j in seq_along(mu)) {
    z <- x - mu[j]
    part <- w[j] * stats::dnorm(z)
    d <- d + part
    s <- s - z * part
  }
  candidates <- x[which.max(d)]
  crossings <- which(s[-length(s)] > 0 & s[-1L] <= 0)
  for (i in crossings) {
    origin <- x[i]
    f_a <- slope(origin)
    f_b <- slope(x[i + 1L])
    # The slopes at the two ends, computed as the search computes them.
    if (!(f_a > 0 && f_b < 0)) {
      candidates <- c(candidates, if (f_a <= 0) origin else x[i + 1L])
      next
    }
    root <- stats::uniroot(function(delta) slope(origin + delta),
                           c(0, x[i + 1L] - origin), f.lower = f_a, f.upper = f_b,
                           tol = 1e-10)$root
    candidates <- c(candidates, origin + root)
  }
  # Refining sampled peaks also covers flat peaks whose derivative has
  # rounded to zero and modes very close to a change in modality.
  peaks <- which(d > 0 & d >= c(-Inf, utils::head(d, -1L)) &
                   d >= c(utils::tail(d, -1L), -Inf))
  for (i in peaks) {
    a <- x[max(1L, i - 1L)]
    b <- x[min(length(x), i + 1L)]
    if (a == b) next
    peak <- stats::optimize(function(delta) density(a + delta), c(0, b - a),
                            maximum = TRUE, tol = 1e-10)$maximum
    candidates <- c(candidates, a + peak)
  }
  candidates <- sort(unique(candidates))
  best <- candidates[which.max(vapply(candidates, density, numeric(1L)))]
  mixture$origin + mixture$sd * best
}

# How quantiles and probabilities are read from point masses when the caller
# does not say. "interpolated" spreads the mass of a grid point evenly over
# the cell of one spacing around it; "step" leaves it on the point.
.bef_point_mass_reading <- function() {
  "interpolated"
}

.bef_match_reading <- function(reading) {
  if (is.null(reading)) {
    return(.bef_point_mass_reading())
  }
  match.arg(reading, c("interpolated", "step"))
}

# P(theta[i] <= q) under the averaged posterior: one row per site and one
# column per value of `q`. `weights` are those of .bef_site_weights() for
# the same sites.
.bef_site_cdf <- function(engine, q, weights = NULL, sites = NULL, reading = NULL) {
  sites <- .bef_engine_sites(engine, sites)
  if (is.null(weights)) {
    weights <- .bef_site_weights(engine, sites = sites)
  }
  q <- as.numeric(q)
  if (engine$tau_c > 0) {
    comp <- .bef_engine_components(engine, sites)
    out <- vapply(q, function(value) {
      rowSums(weights * stats::pnorm((value - comp$mean) / comp$sd))
    }, numeric(length(sites)))
  } else {
    below <- .bef_point_mass_below(engine, q, .bef_match_reading(reading))
    out <- weights %*% below
  }
  matrix(pmin(pmax(out, 0), 1), nrow = length(sites), ncol = length(q))
}

# For point masses: the share of the mass of each grid point that lies at or
# below each value of `q`. One row per grid point and one column per value.
.bef_point_mass_below <- function(engine, q, reading) {
  grid <- engine$grid
  if (identical(reading, "step")) {
    return(outer(grid, q, function(point, value) as.numeric(point <= value)))
  }
  d <- engine$spacing
  outer(grid, q, function(point, value) pmin(pmax((value - point) / d + 0.5, 0), 1))
}

# Quantiles of the averaged posterior: one row per site and one column per
# probability.
.bef_site_quantile <- function(engine, probs, weights = NULL, sites = NULL,
                               reading = NULL) {
  sites <- .bef_engine_sites(engine, sites)
  if (is.null(weights)) {
    weights <- .bef_site_weights(engine, sites = sites)
  }
  K <- length(sites)
  L <- engine$L
  if (engine$tau_c > 0) {
    # The CDF is continuous and increasing, so bisection finds the quantile.
    comp <- .bef_engine_components(engine, sites)
    out <- vapply(probs, function(p) {
      lower <- apply(comp$mean, 1L, min) - 10 * comp$sd
      upper <- apply(comp$mean, 1L, max) + 10 * comp$sd
      for (step in seq_len(60L)) {
        middle <- (lower + upper) / 2
        below <- rowSums(weights * stats::pnorm((middle - comp$mean) / comp$sd)) < p
        lower[below] <- middle[below]
        upper[!below] <- middle[!below]
      }
      (lower + upper) / 2
    }, numeric(K))
    return(matrix(out, nrow = K, ncol = length(probs)))
  }

  reading <- .bef_match_reading(reading)
  cumulative <- t(apply(weights, 1L, cumsum))
  if (L == 1L) cumulative <- matrix(cumulative, ncol = 1L)
  rows <- seq_len(K)
  out <- vapply(probs, function(p) {
    # The first grid point at which the cumulative probability reaches p.
    # The tolerance keeps a probability that should equal p from falling on
    # either side by rounding.
    index <- pmin(rowSums(cumulative < p - 1e-12) + 1L, L)
    if (identical(reading, "step")) {
      return(engine$grid[index])
    }
    before <- ifelse(index > 1L, cumulative[cbind(rows, pmax(index - 1L, 1L))], 0)
    mass <- weights[cbind(rows, index)]
    share <- ifelse(mass > 0, pmin(pmax((p - before) / mass, 0), 1), 0.5)
    engine$grid[index] + engine$spacing * (share - 0.5)
  }, numeric(K))
  matrix(out, nrow = K, ncol = length(probs))
}

# P(theta[i] > threshold) under the averaged posterior, one value per site.
# With point masses read as steps, a grid point equal to the threshold does
# not count as above it.
.bef_site_exceedance <- function(engine, threshold, weights = NULL, sites = NULL,
                                 reading = NULL) {
  1 - .bef_site_cdf(engine, threshold, weights = weights, sites = sites,
                    reading = reading)[, 1L]
}

# The shortest interval that holds probability `level` of the averaged
# posterior of each site: a matrix with one row per site and the columns
# lower and upper. With normal components its length is within one
# hundredth of the shortest (.bef_normal_mixture_hpd()). For point masses
# it is the shortest interval of the reading: between two grid points for
# "step", and for "interpolated" the shortest one under the distribution
# function that rises evenly within each cell.
.bef_site_hpd <- function(engine, level, weights = NULL, sites = NULL,
                          reading = NULL) {
  sites <- .bef_engine_sites(engine, sites)
  if (is.null(weights)) {
    weights <- .bef_site_weights(engine, sites = sites)
  }
  if (engine$tau_c > 0) {
    comp <- .bef_engine_components(engine, sites)
    out <- t(vapply(seq_along(sites), function(k) {
      mixture <- .bef_normal_mixture_parts(comp$mean[k, ], comp$sd[k], weights[k, ])
      .bef_normal_mixture_hpd(mixture, level)
    }, numeric(2L)))
    colnames(out) <- c("lower", "upper")
    return(out)
  }
  shortest <- if (identical(.bef_match_reading(reading), "step")) {
    function(mass) .bef_point_mass_hpd_step(engine$grid, mass, level)
  } else {
    function(mass) .bef_point_mass_hpd_cells(engine$grid, engine$spacing, mass, level)
  }
  out <- t(vapply(seq_along(sites), function(k) shortest(weights[k, ]), numeric(2L)))
  colnames(out) <- c("lower", "upper")
  out
}

# The shortest run of grid points whose probabilities add up to `level` or
# more: its first and last point. Of the runs with that number of points
# the one with the largest probability is taken, and of those that are
# equal in that too, the middle one (the lower of two).
.bef_point_mass_hpd_step <- function(grid, mass, level) {
  L <- length(grid)
  cumulative <- c(0, cumsum(mass)) / sum(mass)
  # For every first point, the last point that is needed. The tolerance
  # keeps rounding from asking for one point more.
  last <- findInterval(cumulative[seq_len(L)] + level - 1e-12, cumulative[-1L],
                       left.open = TRUE) + 1L
  first <- which(last <= L)
  if (length(first) == 0L) {
    return(c(grid[1L], grid[L]))
  }
  points <- last[first] - first
  first <- first[points == min(points)]
  held <- cumulative[last[first] + 1L] - cumulative[first]
  first <- first[held >= max(held) - 1e-12]
  first <- first[ceiling(length(first) / 2)]
  c(grid[first], grid[last[first]])
}

# The shortest interval with probability `level` when the probability of
# each grid point is spread evenly over its cell of one spacing. The
# distribution function is then a broken line, and the length of the
# interval that starts at a value is a broken line too, with its corners
# where the lower or the upper limit is at the edge of a cell. So the
# shortest interval is found among those. Where the intervals between two
# corners are all equally short, as they are when the two limits move
# through cells of equal probability, the middle one is taken, and a
# symmetric distribution has a symmetric interval. Of two stretches of
# shortest intervals that lie apart, the lower one is taken. Where a cell
# has no probability the length jumps, and a shortest corner is taken.
#
# Positions are counted in cells from the lower edge of the first cell, so
# that the search does not depend on where the grid lies or on its spacing.
# A probability within 1e-12 of the value of the distribution function at an
# edge is taken to be that value: sums such as 0.1 + 0.8 need not be 0.9.
.bef_point_mass_hpd_cells <- function(grid, spacing, mass, level) {
  L <- length(grid)
  if (L == 1L || !(spacing > 0)) {
    return(rep(grid[1L], 2L))
  }
  cumulative <- c(0, cumsum(mass)) / sum(mass)
  near <- 1e-12
  share <- function(p, from, to) pmin(pmax((p - from) / (to - from), 0), 1)
  # The distribution function at a position, the first position at which
  # it reaches a probability, and the last position at which it still has it.
  at <- function(u) {
    j <- pmin(pmax(floor(u) + 1L, 1L), L)
    cumulative[j] + (cumulative[j + 1L] - cumulative[j]) * (u - (j - 1L))
  }
  reaches <- function(p) {
    j <- pmin(pmax(findInterval(p - near, cumulative, left.open = TRUE) + 1L, 2L), L + 1L)
    (j - 2L) + share(p, cumulative[j - 1L], cumulative[j])
  }
  leaves <- function(p) {
    j <- pmin(pmax(findInterval(p + near, cumulative), 1L), L)
    (j - 1L) + share(p, cumulative[j], cumulative[j + 1L])
  }
  from_edge <- which(cumulative + level <= 1 + near)
  to_edge <- which(cumulative - level >= -near)
  lower <- c(from_edge - 1L, leaves(pmax(cumulative[to_edge] - level, 0)))
  upper <- c(reaches(pmin(cumulative[from_edge] + level, 1)), to_edge - 1L)
  usable <- is.finite(lower) & is.finite(upper)
  lower <- lower[usable]
  upper <- upper[usable]
  by_lower <- order(lower)
  length_from <- (upper - lower)[by_lower]
  upper <- upper[by_lower]
  lower <- lower[by_lower]
  slack <- 1e-9
  shortest <- length_from <= min(length_from) + slack
  # The first stretch of corners at which the interval is shortest.
  first <- which(shortest)[1L]
  beyond <- which(!shortest & seq_along(shortest) > first)
  last <- if (length(beyond) == 0L) length(shortest) else beyond[1L] - 1L
  start <- (lower[first] + lower[last]) / 2
  limits <- c(start, reaches(min(at(start) + level, 1)))
  if (!(diff(limits) <= min(length_from) + slack)) {
    corner <- which(shortest)
    corner <- corner[ceiling(length(corner) / 2)]
    limits <- c(lower[corner], upper[corner])
  }
  grid[1L] - spacing / 2 + spacing * limits
}

# The shortest interval [x[a], x[b]] with cdf[b] - cdf[a] >= level, for a
# distribution function given at increasing values x.
.bef_shortest_interval <- function(x, cdf, level) {
  n <- length(x)
  # For each start a, the first b with cdf[b] >= cdf[a] + level. The
  # tolerance keeps rounding from pushing b one step further.
  b <- findInterval(cdf + level - 1e-12, cdf, left.open = TRUE) + 1L
  ok <- b <= n
  if (!any(ok)) {
    return(c(x[1L], x[n]))
  }
  width <- x[b[ok]] - x[which(ok)]
  best <- which(ok)[which.min(width)]
  c(x[best], x[b[best]])
}

# The shortest interval that holds a share `level` of a sample of draws.
.bef_draws_hpd <- function(draws, level) {
  sorted <- sort(as.numeric(draws))
  n <- length(sorted)
  keep <- max(1L, min(n, ceiling(level * n)))
  if (keep >= n) {
    return(c(sorted[1L], sorted[n]))
  }
  starts <- seq_len(n - keep + 1L)
  width <- sorted[starts + keep - 1L] - sorted[starts]
  best <- which.min(width)
  c(sorted[best], sorted[best + keep - 1L])
}

# Joint draws of the site effects: one row per draw and one column per
# site. Each row takes one posterior draw of g and, given it, draws the
# component of every site and then the effect within the component.
.bef_simulate_sites <- function(engine, nsim) {
  K <- engine$K
  L <- engine$L
  n_draws <- nrow(engine$log_g)
  like <- .bef_engine_likelihood(engine, seq_len(K))
  picked <- sample.int(n_draws, nsim, replace = nsim > n_draws)
  out <- matrix(NA_real_, nrow = nsim, ncol = K)
  shrink <- if (engine$tau_c > 0) {
    engine$tau_c^2 / (engine$tau_c^2 + engine$sigma^2)
  } else {
    NULL
  }
  for (r in seq_len(nsim)) {
    log_post <- sweep(like$log_lik, 2L, engine$log_g[picked[r], ], "+")
    weights <- exp(log_post - apply(log_post, 1L, max))
    cumulative <- weights
    for (j in seq_len(L)[-1L]) {
      cumulative[, j] <- cumulative[, j - 1L] + weights[, j]
    }
    u <- stats::runif(K) * cumulative[, L]
    index <- pmin(rowSums(cumulative < u) + 1L, L)
    center <- engine$grid[index]
    out[r, ] <- if (is.null(shrink)) {
      center
    } else {
      stats::rnorm(K, center + shrink * (engine$theta_hat - center), sqrt(shrink) * engine$sigma)
    }
  }
  out
}
