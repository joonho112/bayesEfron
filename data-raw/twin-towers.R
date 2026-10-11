# Package one simulation panel of Lee and Sui (2025).
# Rscript data-raw/twin-towers.R <public-panel.rds> [audit-dir]
# Public source: https://github.com/joonho112/fully-bayes-efron-prior
# The source effects are deconvolveR::disjointTheta. I = .7 and R = 9.
# With `audit-dir` the script writes a one-row check of the panel there.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || !file.exists(args[1L])) stop("Supply the public simulation panel RDS.")
panels <- readRDS(args[1L])
selected <- which(panels$I == 0.7 & panels$R == 9)
stopifnot(length(selected) == 1L)
d <- as.data.frame(panels$data[[selected]])
twin_towers <- data.frame(site = seq_len(nrow(d)), theta_true = d$theta_true,
                          yi = d$theta_hat, sei = d$se, vi = d$se2,
                          I = 0.7, R = 9)
stopifnot(nrow(twin_towers) == 1500L, !anyNA(twin_towers),
          all(is.finite(twin_towers$yi)), all(twin_towers$sei > 0),
          max(abs(twin_towers$sei^2 - twin_towers$vi)) < 1e-12,
          abs(max(twin_towers$sei) / min(twin_towers$sei) - 9) < 1e-12)
truth_difference <- NA_real_
if (requireNamespace("deconvolveR", quietly = TRUE)) {
  e <- new.env(parent = emptyenv())
  utils::data("disjointTheta", package = "deconvolveR", envir = e)
  truth_difference <- max(abs(twin_towers$theta_true - e$disjointTheta))
  stopifnot(truth_difference < 1e-12)
}
attr(twin_towers, "source") <- list(
  reference = "Lee and Sui (2025), doi:10.3390/math13162639",
  true_effects = "deconvolveR::disjointTheta (Twin Towers)",
  source_url = "https://github.com/joonho112/fully-bayes-efron-prior",
  source_commit = "f0a57a0235397a420fb5d1562934abbc3f725886",
  source_filename = "data/simulated_data_varied_by_I and R.rds",
  source_sha256 = digest::digest(file = args[1L], algo = "sha256"),
  generation_seed = 12345L, panel_index = selected,
  generation = "Historical public panel; normal errors; shuffled log-spaced variances",
  reliability = "I = 1/(1 + geometric mean(sei^2))",
  heterogeneity = "R = largest/smallest SE (variance ratio R^2)",
  I = 0.7, R = 9
)
dir.create("data", showWarnings = FALSE)
save(twin_towers, file = "data/twin_towers.rda", compress = "xz", version = 3)
if (length(args) >= 2L && nzchar(args[2L])) {
  dir.create(args[2L], recursive = TRUE, showWarnings = FALSE)
  write.csv(data.frame(panel = selected, K = nrow(twin_towers), I = 0.7, R = 9,
                       true_effect_max_difference = truth_difference,
                       geometric_mean_variance = exp(mean(log(twin_towers$vi))),
                       se_ratio = max(twin_towers$sei) / min(twin_towers$sei)),
            file.path(args[2L], "twin-towers-parity.csv"), row.names = FALSE)
}
cat(sprintf("Twin Towers: panel %d, %d sites, I=.7, SE ratio=9, truth difference %.3g\n",
            selected, nrow(twin_towers), truth_difference))
