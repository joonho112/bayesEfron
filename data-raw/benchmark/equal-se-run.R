#!/usr/bin/env Rscript
# Equal standard errors, repeated: the true effects of the independent
# panels of 1,500 sites, each estimated with a normal error of standard
# deviation 1. For every panel the model is fitted and Efron's estimator is
# computed on the grid and with the number of spline functions of that fit.
#
#   Rscript data-raw/benchmark/equal-se-run.R <M> [first] [last] [outdir]
#
#   M            the number of spline functions, 6 or 8
#   first, last  the panels; default 1 to 20
#   outdir       default data-raw/benchmark/results-v040
#
# Run from the package root; needs CmdStan and deconvolveR. The data of
# panel p are made from the seed 20261200 + p and the fit has the seed
# 20261300 + p, four chains of 1,000 warmup and 3,000 sampling iterations,
# one draw in fifteen kept. A panel that already has rows in the output is
# skipped. Output: equal-se-M<M>.csv, three rows for each panel (the fit,
# its estimate of g taken as known, and Efron's estimator).

args <- commandArgs(trailingOnly = TRUE)
M <- as.integer(args[[1L]])
first <- if (length(args) >= 2L) as.integer(args[[2L]]) else 1L
last <- if (length(args) >= 3L) as.integer(args[[3L]]) else 20L
outdir <- if (length(args) >= 4L) args[[4L]] else file.path("data-raw", "benchmark", "results-v040")
suppressMessages(devtools::load_all(".", quiet = TRUE))
stopifnot(requireNamespace("deconvolveR", quietly = TRUE))
panels <- readRDS(file.path("data-raw", "benchmark", "panels_K1500.rds"))$panels
csv <- file.path(outdir, sprintf("equal-se-M%d.csv", M))
done <- if (file.exists(csv)) unique(utils::read.csv(csv)$panel) else integer()
probs <- c(0.05, 0.95)

under_one_g <- function(grid, g, x, s) {
  engine <- .bef_engine_from_parts(log_g = matrix(log(g / sum(g)), nrow = 1L), grid = grid,
                                   theta_hat = x, sigma = s)
  weights <- .bef_site_weights(engine)
  limits <- .bef_site_quantile(engine, probs, weights = weights)
  list(mean = .bef_site_moments(engine)$mean, lower = limits[, 1L], upper = limits[, 2L],
       engine = engine)
}
population_cdf <- function(x) {
  (1 / 3) * stats::punif(x, -1.7, -0.7) + (2 / 3) * stats::punif(x, 0.7, 2.7)
}

for (p in setdiff(seq.int(first, last), done)) {
  truth <- panels[[p]]$theta_true
  K <- length(truth)
  set.seed(20261200L + p)
  z <- truth + stats::rnorm(K)
  ones <- rep(1, K)
  fit <- suppressWarnings(suppressMessages(
    bayes_efron_fit(z, ones, M = M, seed = 20261300L + p, thin = 15L)
  ))
  d <- attr(fit$metadata, "diagnostics", exact = TRUE)
  core <- c("lambda", "mean_g", "sd_g", grep("^alpha\\[", names(d$rhat), value = TRUE))
  tau <- fit$metadata$data_list$grid
  sites <- site_summary(fit)
  g_fit <- as_normalmix(fit)$pi
  Q <- scale(splines::ns(tau, M), center = TRUE, scale = FALSE)
  Q <- apply(Q, 2L, function(w) w / sqrt(sum(w * w)))
  P <- outer(z, tau, function(value, point) stats::dnorm(value - point))
  g_efron <- deconvolveR::deconv(tau = tau, P = P, Q = Q, y = 1, family = "Binomial",
                                 c0 = 1)$stats[, "g"]
  rows <- list(
    "bayesEfron" = list(mean = sites$mean, lower = sites$lower, upper = sites$upper, g = g_fit),
    "g of bayesEfron, plugged in" = c(under_one_g(tau, g_fit, z, ones), list(g = g_fit)),
    "Efron's estimator" = c(under_one_g(tau, g_efron, z, ones), list(g = g_efron))
  )
  x <- seq(min(tau[1L], -12), max(tau[length(tau)], 12), length.out = 40001L)
  engine <- rows[[2L]]$engine
  out <- do.call(rbind, lapply(names(rows), function(method) {
    row <- rows[[method]]
    cdf <- as.numeric((row$g / sum(row$g)) %*% .bef_point_mass_below(engine, x, "interpolated"))
    inside <- truth >= row$lower & truth <= row$upper
    data.frame(
      panel = p, M = M, method = method, K = K,
      lower = tau[1L], upper = tau[length(tau)],
      rmse = sqrt(mean((row$mean - truth)^2)), inclusion = mean(inside),
      mean_length = mean(row$upper - row$lower),
      w1_population = sum(abs(cdf - population_cdf(x))) * (x[2L] - x[1L]),
      max_rhat_parameters = max(d$rhat[core]), min_ess_bulk_parameters = min(d$ess_bulk[core]),
      divergences = sum(d$divergences),
      treedepth_share = sum(d$max_treedepth) / (4 * 3000),
      sampling_seconds = fit$metadata$runtime_seconds, stringsAsFactors = FALSE
    )
  }))
  utils::write.table(out, csv, sep = ",", row.names = FALSE, col.names = !file.exists(csv),
                     append = file.exists(csv), qmethod = "double")
  cat(sprintf("[%s] M = %d panel %d: inclusion %.3f (fit), %.3f (plug-in), %.3f (Efron); R-hat %.3f\n",
              format(Sys.time(), "%H:%M:%S"), M, p, out$inclusion[1L], out$inclusion[2L],
              out$inclusion[3L], out$max_rhat_parameters[1L]))
  rm(fit)
  invisible(gc(verbose = FALSE))
}
