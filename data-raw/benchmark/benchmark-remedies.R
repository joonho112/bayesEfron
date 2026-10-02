#!/usr/bin/env Rscript
# The benchmark fits with 1,500 sites that did not converge with the default
# settings, refitted with a narrower grid and with a larger maximum tree
# depth. Run benchmark-coverage.R first, from the package root: this script
# reads its per-replication results to find the fits with an R-hat above
# 1.05.
#
# Usage:
#   Rscript --vanilla data-raw/benchmark/benchmark-remedies.R [outdir]
#     outdir  default data-raw/benchmark/results
#
# Each fit with the larger tree depth takes about a quarter of an hour.
# Afterwards copy benchmark-remedies.csv from the output directory to
# inst/verification/, where vignette("implementation") reads it.
suppressMessages(devtools::load_all(".", quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[[1L]] else file.path("data-raw", "benchmark", "results")
K <- 1500L
SEED <- 20260509L               # as in benchmark-coverage.R: replication r uses SEED + r

per_replication <- utils::read.csv(file.path(outdir, "benchmark-per-replication.csv"))
default_rows <- per_replication[per_replication$K == K & per_replication$max_rhat > 1.05, ]
fixture <- readRDS(file.path("data-raw", "benchmark", sprintf("lee_sui_K%d.rds", K)))

refit <- function(r, setting, ...) {
  rp <- fixture$replications[[r]]
  arguments <- utils::modifyList(
    list(theta_hat = rp$theta_hat, sigma = rp$sigma, L = rp$L, M = rp$M,
         grid_method = rp$grid_method, expansion = rp$expansion, seed = SEED + r),
    list(...)
  )
  fit <- suppressWarnings(do.call(bayes_efron_fit, arguments))
  intervals <- confint(fit)
  diagnostics <- diagnose(fit)
  data.frame(
    replication = r,
    setting = setting,
    max_rhat = max(diagnostics$rhat, na.rm = TRUE),
    min_ess_bulk = min(diagnostics$ess_bulk, na.rm = TRUE),
    treedepth_hits = sum(diagnostics$max_treedepth),
    coverage = mean(intervals$lower <= rp$theta_true & rp$theta_true <= intervals$upper)
  )
}

rows <- list(data.frame(
  replication = default_rows$replication,
  setting = "default",
  max_rhat = default_rows$max_rhat,
  min_ess_bulk = default_rows$min_ess_bulk,
  treedepth_hits = default_rows$treedepth_hits,
  coverage = default_rows$coverage
))
for (r in default_rows$replication) {
  rows[[length(rows) + 1L]] <- refit(r, "expansion = 0.05", expansion = 0.05)
  rows[[length(rows) + 1L]] <- refit(r, "max_treedepth = 12", max_treedepth = 12L)
  message("replication ", r, " done")
}
out <- do.call(rbind, rows)
out <- out[order(out$replication, match(out$setting,
                 c("default", "expansion = 0.05", "max_treedepth = 12"))), ]
rownames(out) <- NULL
out$generated_on <- format(Sys.Date(), "%Y-%m-%d")
out$cmdstan_version <- as.character(cmdstanr::cmdstan_version())
utils::write.csv(out, file.path(outdir, "benchmark-remedies.csv"), row.names = FALSE)
print(out)
