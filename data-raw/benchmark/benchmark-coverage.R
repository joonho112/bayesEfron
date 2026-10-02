#!/usr/bin/env Rscript
# Coverage of the 90 percent posterior intervals on the benchmark of Lee and
# Sui (2025): 50, 100, 200, 500 and 1,500 sites, 20 data sets of each size.
# The result is the table in inst/verification/benchmark-coverage.csv.
#
# Run this from the package root when the Stan model changes. The run takes
# an hour to two hours at the default settings. Afterwards copy
# benchmark-coverage.csv from the output directory to inst/verification/,
# where vignette("implementation") and vignette("diagnostics") read it.
#
# Usage:
#   Rscript --vanilla data-raw/benchmark/benchmark-coverage.R [outdir] [reps] [config]
#     outdir  default data-raw/benchmark/results
#     reps    default 20
#     config  "default" (4 chains, 1000 warmup, 3000 sampling) or
#             "reduced" (2 chains, 500 warmup, 500 sampling; a quick rehearsal)
#
# Input: tests/testthat/_fixtures/lee_sui_K50.rds and
#        data-raw/benchmark/lee_sui_K{100,200,500,1500}.rds
suppressMessages(devtools::load_all(".", quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[[1L]] else file.path("data-raw", "benchmark", "results")
reps <- if (length(args) >= 2L) as.integer(args[[2L]]) else 20L
config <- if (length(args) >= 3L) args[[3L]] else "default"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cfg <- switch(
  config,
  default = list(chains = 4L, warmup = 1000L, sampling = 3000L),
  reduced = list(chains = 2L, warmup = 500L, sampling = 500L),
  stop("config must be \"default\" or \"reduced\"")
)

Ks <- c(50L, 100L, 200L, 500L, 1500L)
PROBS <- c(0.05, 0.95)          # nominal 90 percent
SEED <- 20260509L               # sampler seed; replication r uses SEED + r

message("Benchmark coverage: config=", config, " reps=", reps,
        " K=", paste(Ks, collapse = ","))
started <- Sys.time()

rows <- list()
for (K in Ks) {
  fixture_dir <- if (K == 50L) {
    file.path("tests", "testthat", "_fixtures")
  } else {
    file.path("data-raw", "benchmark")
  }
  fx <- readRDS(file.path(fixture_dir, sprintf("lee_sui_K%d.rds", K)))
  n_rep <- min(reps, length(fx$replications))
  for (r in seq_len(n_rep)) {
    rp <- fx$replications[[r]]
    t0 <- Sys.time()
    fit <- tryCatch(
      suppressWarnings(bayes_efron_fit(
        theta_hat = rp$theta_hat, sigma = rp$sigma,
        L = rp$L, M = rp$M, grid_method = rp$grid_method,
        expansion = rp$expansion,
        chains = cfg$chains, iter_warmup = cfg$warmup,
        iter_sampling = cfg$sampling, seed = SEED + r
      )),
      error = function(err) {
        message("  FAIL K=", K, " rep=", r, ": ", conditionMessage(err))
        NULL
      }
    )
    if (is.null(fit)) next

    theta_rep <- .bef_site_draws(fit, "theta_rep")
    q <- apply(theta_rep, 2, stats::quantile, probs = PROBS, names = FALSE)
    covered <- rp$theta_true >= q[1, ] & rp$theta_true <= q[2, ]
    d <- diagnose(fit)

    rows[[length(rows) + 1L]] <- data.frame(
      K = K, replication = r,
      coverage = mean(covered),
      n_sites = length(covered),
      mean_width = mean(q[2, ] - q[1, ]),
      max_rhat = suppressWarnings(max(d$rhat, na.rm = TRUE)),
      min_ess_bulk = suppressWarnings(min(d$ess_bulk, na.rm = TRUE)),
      divergences = sum(d$divergences, na.rm = TRUE),
      treedepth_hits = sum(d$max_treedepth, na.rm = TRUE),
      sampler_seconds = fit$metadata$runtime_seconds,
      wall_seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
      stringsAsFactors = FALSE
    )
    rm(fit, theta_rep); invisible(gc(verbose = FALSE))
  }
  done <- do.call(rbind, rows)
  sub <- done[done$K == K, ]
  message(sprintf("  K=%4d done: %d reps, aggregate coverage %.6f, %.1f min elapsed",
                  K, nrow(sub),
                  sum(sub$coverage * sub$n_sites) / sum(sub$n_sites),
                  as.numeric(difftime(Sys.time(), started, units = "mins"))))
  utils::write.csv(done, file.path(outdir, "benchmark-per-replication.csv"),
                   row.names = FALSE)
}

per_rep <- do.call(rbind, rows)
rownames(per_rep) <- NULL

agg <- do.call(rbind, lapply(split(per_rep, per_rep$K), function(d) {
  # site-weighted mean over replications
  cov <- sum(d$coverage * d$n_sites) / sum(d$n_sites)
  data.frame(
    K = d$K[1], n_replications = nrow(d), n_sites_total = sum(d$n_sites),
    coverage = cov,
    coverage_se = stats::sd(d$coverage) / sqrt(nrow(d)),
    mean_width = mean(d$mean_width),
    max_rhat = max(d$max_rhat),
    n_rhat_above_1.05 = sum(d$max_rhat > 1.05),
    min_ess_bulk = min(d$min_ess_bulk),
    total_divergences = sum(d$divergences),
    total_treedepth_hits = sum(d$treedepth_hits),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}))
rownames(agg) <- NULL
agg$generated_on <- format(Sys.Date(), "%Y-%m-%d")
agg$cmdstan_version <- as.character(cmdstanr::cmdstan_version())

utils::write.csv(per_rep, file.path(outdir, "benchmark-per-replication.csv"), row.names = FALSE)
utils::write.csv(agg, file.path(outdir, "benchmark-coverage.csv"), row.names = FALSE)
saveRDS(
  list(per_replication = per_rep, aggregate = agg, config = cfg,
       probs = PROBS, seed = SEED,
       stan_sha256 = digest::digest(file = file.path("inst", "stan", "efron_re.stan"),
                                    algo = "sha256")),
  file.path(outdir, "benchmark-run.rds")
)

cat("\nCoverage of the 90 percent intervals\n")
print(agg[, c("K", "n_replications", "coverage", "coverage_se", "max_rhat",
              "n_rhat_above_1.05", "min_ess_bulk", "total_divergences",
              "total_treedepth_hits")],
      row.names = FALSE, digits = 6)
cat("elapsed: ", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1),
    " min\n", sep = "")
cat("written to ", outdir, "\n", sep = "")
