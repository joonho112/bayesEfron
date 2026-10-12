# ==============================================================================
# tail-coverage.R: coverage of the 90% site intervals for the sites in the
# tails of g and for the sites whose estimate lies in the empty middle, in
# the two-group designs (100 data sets of 100 sites, 20 of 1,500), for the
# fit, its posterior mean of g plugged in, the true distribution of the
# population plugged in, and the ebnm normal and NPMLE fits. Reads the
# per-site results of the package's benchmark (data-raw/benchmark/ in the
# source repository) and writes data/tail-coverage.csv and
# data/tail-coverage-numbers.csv; data/tail-coverage-panels.csv retains the
# panel-level results underlying the Monte Carlo standard errors.
# tables.R draws Table E.2 from the aggregate CSV.
#
# Run from the paper folder root:
#   Rscript scripts/tail-coverage.R [path/to/package-checkout] [--check-inputs]
# The default first checks the parent checkout (checkout/paper layout),
# then the two sibling checkout names used for the manuscript source.
# --check-inputs validates paths and the package version without fitting.
# ==============================================================================
root <- normalizePath(".")
.libPaths(c(file.path(root, "r-library"), .libPaths()))
suppressPackageStartupMessages({ library(bayesEfron); library(ebnm) })
print(packageVersion("bayesEfron"))
stopifnot(packageVersion("bayesEfron") == "0.4.0")
args <- commandArgs(trailingOnly = TRUE)
check_inputs <- "--check-inputs" %in% args
args <- args[args != "--check-inputs"]
if (length(args) > 1L || any(startsWith(args, "--")))
  stop("Usage: Rscript scripts/tail-coverage.R [path/to/package-checkout] [--check-inputs]", call. = FALSE)
checkouts <- if (length(args)) path.expand(args[[1L]]) else
  file.path(root, "..", c(".", "bayesEfron-R-package-public", "bayesEfron-R-package"))
candidates <- file.path(checkouts, "data-raw", "benchmark")
required <- c(sprintf("panels_K%d.rds", c(100L, 1500L)),
              file.path("results-v040", sprintf("panels-K%d-default-fits.rds", c(100L, 1500L))))
complete <- vapply(candidates, function(path) all(file.exists(file.path(path, required))), logical(1))
if (!any(complete)) {
  missing <- vapply(candidates, function(path) paste0("  ", path, "\n    missing: ",
    paste(required[!file.exists(file.path(path, required))], collapse = ", ")), character(1))
  stop(paste0("Benchmark inputs not found. Supply the package checkout as the first argument.\n",
              paste(missing, collapse = "\n")), call. = FALSE)
}
bench <- normalizePath(candidates[which(complete)[1L]])
message("benchmark folder: ", bench)
if (check_inputs) {
  message("Input check passed: all four benchmark files exist; no fits or outputs were generated.")
  quit(save = "no", status = 0L)
}
everything <- ebnm::ebnm_output_all()
quietly <- function(expr) withCallingHandlers(expr, warning = function(w) {
  if (grepl("local false sign rates", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
})
limits_ebnm <- function(fit) {
  q <- stats::quantile(fit, probs = c(0.05, 0.95), nsim = 4000L, seed = 1L)
  list(lower = unname(q[, 1]), upper = unname(q[, 2]))
}
plug_in <- function(grid, g, x, s) {
  engine <- bayesEfron:::.bef_engine_from_parts(log_g = matrix(log(g / sum(g)), nrow = 1L),
                                                grid = grid, theta_hat = x, sigma = s, tau_c = 0)
  w <- bayesEfron:::.bef_site_weights(engine)
  q <- bayesEfron:::.bef_site_quantile(engine, c(0.05, 0.95), weights = w)
  list(lower = q[, 1], upper = q[, 2])
}
# The generating distribution in benchmark/panels-make.R is continuous:
# (1/3) Uniform(-1.7, -0.7) + (2/3) Uniform(0.7, 2.7). Its posterior is a
# mixture of two truncated normals. Integrate the likelihood exactly on
# each interval, then invert the component CDF; no empirical population,
# grid approximation or Monte Carlo sampling is used for this comparator.
# Log probabilities avoid cancellation when an estimate lies far from g.
known_g <- function(x, s) {
  stopifnot(length(x) == length(s), all(is.finite(x)),
            all(is.finite(s)), all(s > 0))
  a <- c(-1.7, 0.7); b <- c(-0.7, 2.7)
  log_density <- log(c(1 / 3, 2 / 3) / (b - a))
  log_add <- function(u, v) {
    m <- max(u, v)
    if (is.infinite(m)) return(m)
    m + log1p(exp(min(u, v) - m))
  }
  lower <- upper <- numeric(length(x))
  for (i in seq_along(x)) {
    za <- (a - x[i]) / s[i]; zb <- (b - x[i]) / s[i]
    # Use the upper tail when both endpoints are positive.
    upper_tail <- za > 0
    log_big <- ifelse(upper_tail, stats::pnorm(za, lower.tail = FALSE, log.p = TRUE),
                      stats::pnorm(zb, log.p = TRUE))
    log_small <- ifelse(upper_tail, stats::pnorm(zb, lower.tail = FALSE, log.p = TRUE),
                        stats::pnorm(za, log.p = TRUE))
    log_mass <- log_big + log(-expm1(log_small - log_big))
    log_weight <- log_density + log_mass
    w <- exp(log_weight - max(log_weight)); w <- w / sum(w)
    quantile_at <- function(p) {
      j <- if (p <= w[1L]) 1L else 2L
      before <- if (j == 1L) 0 else w[1L]
      fraction <- (p - before) / w[j]
      if (fraction <= 0) return(a[j])
      if (fraction >= 1) return(b[j])
      lp <- if (upper_tail[j])
        log_add(stats::pnorm(zb[j], lower.tail = FALSE, log.p = TRUE),
                log1p(-fraction) + log_mass[j]) else
        log_add(stats::pnorm(za[j], log.p = TRUE), log(fraction) + log_mass[j])
      x[i] + s[i] * stats::qnorm(lp, lower.tail = !upper_tail[j], log.p = TRUE)
    }
    lower[i] <- quantile_at(0.05)
    upper[i] <- quantile_at(0.95)
  }
  stopifnot(all(is.finite(lower)), all(is.finite(upper)), all(lower <= upper))
  list(lower = lower, upper = upper)
}
# the groups are (-1.7, -0.7) and (0.7, 2.7); the outer tenth at either end
# of a group is a tenth of its width at each end
subsets <- function(x, truth) {
  cut <- stats::quantile(x, c(0.025, 0.975))
  list("All sites" = rep(TRUE, length(x)),
       "Estimate between the groups (-0.7, 0.7)" = abs(x) < 0.7,
       "Estimates outside the sample 2.5th and 97.5th percentiles" = x < cut[1] | x > cut[2],
       "True effect in the outer tenth at either end of its group" =
         truth < -1.6 | (truth > -0.8 & truth < -0.7) | (truth > 0.7 & truth < 0.9) | truth > 2.5)
}
rows <- list()
for (K in c(100L, 1500L)) {
  panels <- readRDS(file.path(bench, sprintf("panels_K%d.rds", K)))$panels
  kept <- readRDS(file.path(bench, "results-v040", sprintf("panels-K%d-default-fits.rds", K)))
  stopifnot(length(kept) == length(panels))
  for (p in seq_along(panels)) {
    pan <- panels[[p]]; st <- kept[[as.character(p)]]
    stopifnot(identical(st$panel, pan$panel), nrow(st$sites) == pan$K)
    x <- pan$theta_hat; s <- pan$sigma; truth <- pan$theta_true
    methods <- list(
      "bayesEfron" = list(lower = st$sites$lower, upper = st$sites$upper),
      "bayesEfron g, plugged in" = plug_in(st$grid, unname(st$g), x, s),
      "True g, plugged in" = known_g(x, s),
      "Normal (ebnm)" = limits_ebnm(quietly(ebnm_normal(x, s, mode = "estimate", output = everything))),
      "NPMLE (ebnm)" = limits_ebnm(quietly(ebnm_npmle(x, s, output = everything)))
    )
    sub <- subsets(x, truth)
    for (m in names(methods)) for (nm in names(sub)) {
      idx <- sub[[nm]]
      inside <- truth[idx] >= methods[[m]]$lower[idx] & truth[idx] <= methods[[m]]$upper[idx]
      rows[[length(rows) + 1L]] <- data.frame(K = K, panel = p, method = m, subset = nm,
                                              n = sum(idx), inclusion = mean(inside))
    }
  }
  message("K = ", K, ": ", length(panels), " panels done")
}
rows <- do.call(rbind, rows)
write.csv(rows, "data/tail-coverage-panels.csv", row.names = FALSE)
agg <- do.call(rbind, lapply(split(rows, list(rows$K, rows$method, rows$subset), drop = TRUE), function(d) {
  data.frame(K = d$K[1], method = d$method[1], subset = d$subset[1], panels = nrow(d),
             sites_per_panel = mean(d$n), inclusion = mean(d$inclusion),
             inclusion_se = stats::sd(d$inclusion) / sqrt(nrow(d)))
}))
method_order <- c("bayesEfron", "bayesEfron g, plugged in", "True g, plugged in", "Normal (ebnm)", "NPMLE (ebnm)")
subset_order <- names(subsets(0, 0))
agg <- agg[order(agg$K, match(agg$subset, subset_order), match(agg$method, method_order)), ]
rownames(agg) <- NULL
write.csv(agg, "data/tail-coverage.csv", row.names = FALSE)
print(agg, digits = 3)


# ---- numbers quoted in the text ----------------------------------------------
get <- function(K, m, nm) agg[agg$K == K & agg$method == m & agg$subset == nm, ]
reg <- list()
for (K in c(100L, 1500L)) {
  for (m in method_order) {
    key <- c("bayesEfron" = "be", "bayesEfron g, plugged in" = "plug", "True g, plugged in" = "true",
             "Normal (ebnm)" = "normal", "NPMLE (ebnm)" = "npmle")[[m]]
    reg[[sprintf("tail_K%d_%s_middle", K, key)]] <- get(K, m, subset_order[2])$inclusion
    reg[[sprintf("tail_K%d_%s_extreme", K, key)]] <- get(K, m, subset_order[3])$inclusion
    reg[[sprintf("tail_K%d_%s_outer", K, key)]] <- get(K, m, subset_order[4])$inclusion
    reg[[sprintf("tail_K%d_%s_middle_se", K, key)]] <- get(K, m, subset_order[2])$inclusion_se
    reg[[sprintf("tail_K%d_%s_extreme_se", K, key)]] <- get(K, m, subset_order[3])$inclusion_se
    reg[[sprintf("tail_K%d_%s_outer_se", K, key)]] <- get(K, m, subset_order[4])$inclusion_se
  }
  reg[[sprintf("tail_K%d_middle_n", K)]] <- get(K, "bayesEfron", subset_order[2])$sites_per_panel
  reg[[sprintf("tail_K%d_extreme_n", K)]] <- get(K, "bayesEfron", subset_order[3])$sites_per_panel
  reg[[sprintf("tail_K%d_outer_n", K)]] <- get(K, "bayesEfron", subset_order[4])$sites_per_panel
}
for (K in c(100L, 1500L)) {
  d <- agg[agg$K == K & agg$method %in% c("bayesEfron", "bayesEfron g, plugged in"), ]
  gap <- abs(d$inclusion[d$method == "bayesEfron"] - d$inclusion[d$method == "bayesEfron g, plugged in"])
  reg[[sprintf("tail_K%d_plug_max_gap_pp", K)]] <- 100 * max(gap)
}
extra <- data.frame(name = names(reg), value = unlist(reg), stringsAsFactors = FALSE)
write.csv(extra, "data/tail-coverage-numbers.csv", row.names = FALSE)
message("done")
