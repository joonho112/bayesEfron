#!/usr/bin/env Rscript
# Independent panels for the benchmark: every panel draws its own true
# effects. The 20 data sets with 1,500 sites in lee_sui_K1500.rds share one
# set of true effects, so they are not independent replications.
#
#   Rscript data-raw/benchmark/panels-make.R [K] [panels] [first-seed]
#
# Defaults: K = 1500, 20 panels, first seed 20261009; panel p is generated
# with the seed first-seed + p. Run from the package root. Writes
# data-raw/benchmark/panels_K<K>.rds; running it again gives the same file.
#
# The data-generating rule is that of the Twin Towers example of Lee and
# Sui (2025) with I = 0.7 and R = 9:
#   - a true effect is uniform on (-1.7, -0.7) with probability 1/3 and
#     uniform on (0.7, 2.7) with probability 2/3, the two populations from
#     which deconvolveR::disjointTheta was drawn (500 and 1,000 values);
#   - the sampling variances are equally spaced on the log scale, with
#     geometric mean (1 - I) / I and the largest R^2 times the smallest, so
#     that the largest standard error is R times the smallest, and they are
#     assigned to the sites at random;
#   - an estimate is its true effect plus a normal error with its standard
#     error.

args <- commandArgs(trailingOnly = TRUE)
K <- if (length(args) >= 1L) as.integer(args[[1L]]) else 1500L
n_panels <- if (length(args) >= 2L) as.integer(args[[2L]]) else 20L
first_seed <- if (length(args) >= 3L) as.integer(args[[3L]]) else 20261009L
reliability <- 0.7
ratio <- 9

make_panel <- function(p) {
  seed <- first_seed + p
  set.seed(seed)
  right <- stats::runif(K) < 2 / 3
  theta_true <- ifelse(right, stats::runif(K, 0.7, 2.7), stats::runif(K, -1.7, -0.7))
  center <- log((1 - reliability) / reliability)
  variances <- exp(seq(center - log(ratio), center + log(ratio), length.out = K))
  sigma <- sqrt(sample(variances))
  theta_hat <- theta_true + stats::rnorm(K, 0, sigma)
  list(panel = p, seed = seed, K = K, theta_hat = theta_hat, sigma = sigma,
       theta_true = theta_true)
}

panels <- lapply(seq_len(n_panels), make_panel)
out <- list(
  description = paste(
    "Independent panels: true effects drawn anew for each panel from the two uniform",
    "populations of the Twin Towers example; I = 0.7, R = 9."
  ),
  generator = "data-raw/benchmark/panels-make.R",
  K = K, n_panels = n_panels, first_seed = first_seed,
  reliability = reliability, ratio = ratio,
  r_version = R.version.string,
  rng = RNGkind(),
  panels = panels
)
file <- file.path("data-raw", "benchmark", sprintf("panels_K%d.rds", K))
saveRDS(out, file, version = 3, compress = "xz")
cat(sprintf("%s: %d panels of %d sites, %d bytes, sha256 %s\n", file, n_panels, K,
            file.size(file),
            digest::digest(lapply(panels, `[`, c("theta_hat", "sigma", "theta_true")),
                           algo = "sha256")))
