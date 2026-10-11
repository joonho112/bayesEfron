# Rebuild the firm-level racial callback gaps from the public Walters
# replication of Kline, Rose and Walters (2022, 2024). Run at package root:
# Rscript data-raw/krw-discrimination.R <krw_data.dta> [reference.rds] [audit-dir]
# The raw data are available through:
# https://sites.google.com/view/christopher-walters/research
# Uses haven only to read Stata data; calculations below use base R.
# With `reference.rds`, the firm estimates of the replication repository of
# Lee and Sui (2025), data/firm_discrimination_estimates.rds, the script
# checks that it reproduces them. With `audit-dir` it writes there the table
# of the firms that shows which were included, and that comparison.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || !file.exists(args[1L])) {
  stop("Supply the path to the public krw_data.dta.")
}
if (!requireNamespace("haven", quietly = TRUE)) stop("Install haven to read the Stata data.")
raw <- as.data.frame(haven::read_dta(args[1L]))
required <- c("firm_id", "firm_name", "job_id", "callback", "white")
stopifnot(all(required %in% names(raw)), !anyNA(raw[required]),
          all(raw$callback %in% 0:1), all(raw$white %in% 0:1))
job_firms <- tapply(raw$firm_id, raw$job_id, function(x) length(unique(x)))
stopifnot(all(job_firms == 1L))
firm_ids <- sort(unique(raw$firm_id))
sample_audit <- do.call(rbind, lapply(firm_ids, function(id) {
  d <- raw[raw$firm_id == id, ]
  data.frame(firm_id = as.numeric(id), callback_rate = mean(d$callback),
             n_jobs = length(unique(d$job_id)), n_obs = nrow(d))
}))
sample_audit$included <- sample_audit$callback_rate >= 0.03 & sample_audit$n_jobs >= 40L
ids <- sample_audit$firm_id[sample_audit$included]
krw_discrimination <- do.call(rbind, lapply(seq_along(ids), function(site) {
  d <- raw[raw$firm_id == ids[site], ]
  X <- cbind(intercept = 1, white = as.numeric(d$white))
  ols <- lm.fit(X, as.numeric(d$callback))
  N <- nrow(X)
  G <- length(unique(d$job_id))
  stopifnot(ols$rank == 2L, N > 2L, G > 1L, length(unique(d$firm_name)) == 1L)
  scores <- rowsum(X * ols$residuals, group = d$job_id, reorder = FALSE)
  bread <- solve(crossprod(X))
  # Stata reg ..., cluster(job_id): CR1 with two fitted coefficients.
  covariance <- (G / (G - 1)) * ((N - 1) / (N - 2)) *
    bread %*% crossprod(scores) %*% bread
  data.frame(site = site, firm_id = ids[site], firm_name = as.character(d$firm_name[1L]),
             yi = unname(ols$coefficients[2L]), sei = sqrt(covariance[2L, 2L]),
             n_jobs = G, n_obs = N)
}))
krw_discrimination$vi <- krw_discrimination$sei^2
krw_discrimination <- krw_discrimination[
  c("site", "firm_id", "firm_name", "yi", "sei", "vi", "n_jobs", "n_obs")
]
rownames(krw_discrimination) <- NULL
stopifnot(nrow(krw_discrimination) == 97L,
          sum(krw_discrimination$n_obs) == 78910L,
          sum(krw_discrimination$n_jobs) == 10453L,
          !anyDuplicated(krw_discrimination$firm_id),
          all(is.finite(krw_discrimination$yi)), all(krw_discrimination$sei > 0))
parity <- NULL
if (length(args) >= 2L && nzchar(args[2L])) {
  reference <- as.data.frame(readRDS(args[2L]))
  ref <- reference[match(krw_discrimination$firm_id, reference$firm_id), ]
  stopifnot(!anyNA(ref$firm_id), identical(as.numeric(ref$firm_id), krw_discrimination$firm_id))
  parity <- data.frame(firm_id = krw_discrimination$firm_id,
                      estimate_difference = krw_discrimination$yi - ref$theta_white,
                      se_difference = krw_discrimination$sei - ref$s_white)
  stopifnot(max(abs(parity$estimate_difference)) < 1e-10,
            max(abs(parity$se_difference)) < 1e-10)
}
attr(krw_discrimination, "source") <- list(
  original_study = "Kline, Rose and Walters (2022), QJE 137(4), 1963-2036",
  analysis_sample = "Kline, Rose and Walters (2024), AER 114(8), 2472-2525",
  replication = "Walters (2024), Empirical Bayes methods in labor economics",
  source_url = "https://sites.google.com/view/christopher-walters/research",
  raw_filename = "krw_data.dta",
  raw_sha256 = digest::digest(file = args[1L], algo = "sha256"),
  estimator = "OLS with intercept; job-clustered CR1 SE, G/(G-1)*(N-1)/(N-2)",
  units = "probability difference, white minus Black",
  inclusion = "mean(callback) >= 0.03 and at least 40 distinct job_id",
  reanalysis_method = "Lee and Sui (2025), doi:10.3390/math13162639"
)
dir.create("data", showWarnings = FALSE)
save(krw_discrimination, file = "data/krw_discrimination.rda", compress = "xz", version = 3)
if (length(args) >= 3L && nzchar(args[3L])) {
  dir.create(args[3L], recursive = TRUE, showWarnings = FALSE)
  write.csv(sample_audit, file.path(args[3L], "krw-sample.csv"), row.names = FALSE)
  if (!is.null(parity)) write.csv(parity, file.path(args[3L], "krw-parity.csv"), row.names = FALSE)
}
cat(sprintf("97 firms; %d applications; %d vacancies; estimate mean %.9f; positive SE: %s\n",
            sum(krw_discrimination$n_obs), sum(krw_discrimination$n_jobs),
            mean(krw_discrimination$yi), all(krw_discrimination$sei > 0)))
if (!is.null(parity)) cat(sprintf("Maximum absolute parity difference: estimate %.3g; SE %.3g\n",
                                 max(abs(parity$estimate_difference)), max(abs(parity$se_difference))))
