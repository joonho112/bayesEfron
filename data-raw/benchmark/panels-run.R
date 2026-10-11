#!/usr/bin/env Rscript
# Fit the independent panels of panels-make.R with one of three settings and
# keep, for every fit, its diagnostics, its cost, how well it recovers the
# true effects, and the summaries that a later comparison needs.
#
#   Rscript data-raw/benchmark/panels-run.R <arm> [outdir] [K] [first] [last]
#
#   arm     "default"  bayes_efron_fit(theta_hat, sigma)
#           "M8"       bayes_efron_fit(theta_hat, sigma, M = 8)
#           "v030"     bayes_efron_fit(theta_hat, sigma, settings = "0.3.0",
#                                      site_draws = FALSE)
#   outdir  default data-raw/benchmark/results-v040
#   K       default 1500 (the file panels_K<K>.rds)
#   first, last  the panels to fit; default all
#
# Run from the package root; needs CmdStan. Every fit has four chains of
# 1,000 warmup and 3,000 sampling iterations and the seed 20261109 + panel.
# A panel that already has a row in the output of its arm is skipped, so the
# script can be stopped and started again. Outputs, per arm:
#   panels-K<K>-<arm>.csv        one row for each fit
#   panels-K<K>-<arm>-fits.rds   for each fit, the estimate of g on its grid
#                                and the table of the sites
# The times are those of this machine while it ran; the numbers of leapfrog
# steps do not depend on what else was running.

args <- commandArgs(trailingOnly = TRUE)
arm <- args[[1L]]
stopifnot(arm %in% c("default", "M8", "v030"))
outdir <- if (length(args) >= 2L) args[[2L]] else file.path("data-raw", "benchmark", "results-v040")
K <- if (length(args) >= 3L) as.integer(args[[3L]]) else 1500L
suppressMessages(devtools::load_all(".", quiet = TRUE))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
panels <- readRDS(file.path("data-raw", "benchmark", sprintf("panels_K%d.rds", K)))$panels
first <- if (length(args) >= 4L) as.integer(args[[4L]]) else 1L
last <- if (length(args) >= 5L) as.integer(args[[5L]]) else length(panels)
csv <- file.path(outdir, sprintf("panels-K%d-%s.csv", K, arm))
rds <- file.path(outdir, sprintf("panels-K%d-%s-fits.rds", K, arm))
done <- if (file.exists(csv)) utils::read.csv(csv)$panel else integer()
kept <- if (file.exists(rds)) readRDS(rds) else list()

arm_args <- switch(
  arm,
  default = list(),
  M8 = list(M = 8L),
  v030 = list(settings = "0.3.0", site_draws = FALSE)
)

# The distribution from which the true effects are drawn.
population_cdf <- function(x) {
  (1 / 3) * stats::punif(x, -1.7, -0.7) + (2 / 3) * stats::punif(x, 0.7, 2.7)
}

for (p in setdiff(seq.int(first, last), done)) {
  panel <- panels[[p]]
  caught <- character()
  started <- Sys.time()
  fit <- tryCatch(
    withCallingHandlers(
      do.call(bayes_efron_fit, c(
        list(panel$theta_hat, panel$sigma, seed = 20261109L + p), arm_args
      )),
      warning = function(w) {
        caught <<- c(caught, class(w)[[1L]])
        invokeRestart("muffleWarning")
      },
      message = function(m) invokeRestart("muffleMessage")
    ),
    error = function(err) err
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (inherits(fit, "error")) {
    row <- data.frame(panel = p, arm = arm, ok = FALSE, error = conditionMessage(fit),
                      stringsAsFactors = FALSE)
    failed <- sub("[.]csv$", "-failed.csv", csv)
    utils::write.table(row, failed, sep = ",", row.names = FALSE,
                       col.names = !file.exists(failed), append = file.exists(failed))
    next
  }
  d <- attr(fit$metadata, "diagnostics", exact = TRUE)
  core <- c("lambda", "mean_g", "sd_g", grep("^alpha\\[", names(d$rhat), value = TRUE))
  nuts <- attr(fit$metadata, "nuts", exact = TRUE)
  grid_record <- attr(fit$metadata, "grid", exact = TRUE)
  engine <- .bef_engine(fit)
  n_draws <- nrow(engine$log_g)
  g <- colMeans(exp(engine$log_g))
  sites <- site_summary(fit)
  truth <- panel$theta_true
  x <- seq(min(engine$grid[1L], -12), max(engine$grid[engine$L], 12), length.out = 40001L)
  cdf <- as.numeric(g %*% .bef_point_mass_below(engine, x, "interpolated"))
  step <- x[2L] - x[1L]
  inside <- truth >= sites$lower & truth <= sites$upper
  ends <- diagnose(fit)$grid$ends
  row <- data.frame(
    panel = p, arm = arm, ok = TRUE, K = panel$K, seed = fit$metadata$seed,
    program = attr(fit$metadata, "program", exact = TRUE),
    rule = fit$metadata$grid_method, L = engine$L, M = fit$metadata$data_list$M,
    lower = engine$grid[1L], upper = engine$grid[engine$L],
    share_of_wide_range = if (is.null(grid_record$pilot)) NA_real_ else grid_record$pilot$share,
    max_rhat_parameters = max(d$rhat[core]), max_rhat = max(d$rhat, na.rm = TRUE),
    min_ess_bulk_parameters = min(d$ess_bulk[core]),
    min_ess_tail_parameters = min(d$ess_tail[core]),
    divergences = sum(d$divergences), treedepth_hits = sum(d$max_treedepth),
    treedepth_share = sum(d$max_treedepth) / n_draws, min_ebfmi = min(d$ebfmi),
    valid = max(d$rhat[core]) < 1.01 && min(d$ess_bulk[core]) >= 400 &&
      sum(d$divergences) == 0,
    rhat_above_1.05 = max(d$rhat[core]) > 1.05,
    leapfrog = if (is.null(nuts)) NA_real_ else sum(nuts[, , "n_leapfrog__"]),
    pilot_seconds = if (is.null(grid_record$pilot)) 0 else grid_record$pilot$seconds,
    sampling_seconds = fit$metadata$runtime_seconds, elapsed_seconds = elapsed,
    mean_g = mean(fit$posterior$mean_g), sd_g = mean(fit$posterior$sd_g),
    w1_panel = sum(abs(cdf - stats::ecdf(truth)(x))) * step,
    w1_population = sum(abs(cdf - population_cdf(x))) * step,
    rmse = sqrt(mean((sites$mean - truth)^2)),
    coverage_90 = mean(inside), mean_length_90 = mean(sites$upper - sites$lower),
    end_mass = max(ends$mass), farthest_estimate_beyond_se = max(ends$estimate_beyond),
    warnings = paste(unique(caught), collapse = ";"),
    stringsAsFactors = FALSE
  )
  utils::write.table(row, csv, sep = ",", row.names = FALSE, col.names = !file.exists(csv),
                     append = file.exists(csv), qmethod = "double")
  kept[[as.character(p)]] <- list(
    panel = p, arm = arm, grid = engine$grid, g = g,
    sites = sites[, c("site", "mean", "sd", "lower", "upper")]
  )
  saveRDS(kept, rds, version = 3, compress = "xz")
  cat(sprintf("[%s] %s panel %d: R-hat %.3f, tree depth %.1f%%, %.0f s sampling, W1 %.3f\n",
              format(Sys.time(), "%H:%M:%S"), arm, p, row$max_rhat_parameters,
              100 * row$treedepth_share, row$sampling_seconds, row$w1_panel))
  rm(fit, engine)
  invisible(gc(verbose = FALSE))
}
