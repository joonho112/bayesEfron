# Builds inst/examples/lee_sui_examples.rds: numbers that the documentation
# reads for the callback gaps and the simulated sites, in the same units as
# the data. Needs CmdStan and the metafor package. Run from the package root
# after data-raw/example-fits.R:
#
#   Rscript data-raw/example-evidence.R [audit-dir]
#
# The table `sensitivity` (three fits of the callback gaps and three rows
# for the simulated sites) is also written as a CSV file to `audit-dir`, or
# to a temporary directory if none is given.
#
# For the callback gaps the default fit is set beside two others at full
# length (12,000 draws): a narrower range and eight spline functions. For
# the simulated sites the rows come from the comparison that
# data-raw/example-fits.R stored. Each case also holds the normal
# random-effects model fitted by REML to the same estimates, and the
# probability of a positive effect for every site.

args <- commandArgs(trailingOnly = TRUE)
audit_dir <- if (length(args) >= 1L) args[[1L]] else file.path(tempdir(), "bayesEfron-example-fits")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
suppressMessages(devtools::load_all(".", quiet = TRUE))
stopifnot(requireNamespace("metafor", quietly = TRUE))

quiet_fit <- function(...) {
  suppressWarnings(suppressMessages(bayes_efron_fit(...), classes = "bef_default_changed"))
}
diagnostics_of <- function(fit) {
  d <- attr(fit$metadata, "diagnostics", exact = TRUE)
  data.frame(
    K = nobs(fit),
    max_rhat = max(d$rhat, na.rm = TRUE), min_ess_bulk = min(d$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(d$ess_tail, na.rm = TRUE), divergences = sum(d$divergences),
    treedepth_hits = sum(d$max_treedepth), min_ebfmi = min(d$ebfmi),
    sampler_seconds = fit$metadata$runtime_seconds,
    mean_g = mean(fit$posterior$mean_g), sd_g = mean(fit$posterior$sd_g)
  )
}
settings_of <- function(fit, case, spec, data) {
  settings <- attr(fit$metadata, "sampler_settings", exact = TRUE)
  grid <- fit$metadata$data_list$grid
  data.frame(
    case = case, spec = spec, data_sha256 = digest::digest(data, algo = "sha256"),
    seed = fit$metadata$seed, L = fit$metadata$data_list$L, M = fit$metadata$data_list$M,
    grid_method = fit$metadata$grid_method, lower = min(grid), upper = max(grid),
    chains = settings$chains, warmup = settings$iter_warmup, sampling = settings$iter_sampling,
    adapt_delta = settings$adapt_delta, max_treedepth = settings$max_treedepth,
    stan_sha256 = fit$metadata$stan_file_sha256, stringsAsFactors = FALSE
  )
}
# The normal random-effects model with its mean and variance held at their
# REML estimates.
normal_of <- function(data) {
  # rma() stops when the variance of the effects changes by less than 1e-5
  # from one iteration to the next. For the callback gaps, whose variance is
  # near 2e-4, that is too early, so the threshold is tightened. The result
  # must agree with the estimate that plot_g(normal = TRUE) draws.
  model <- metafor::rma(yi = data$yi, vi = data$vi, method = "REML",
                        control = list(threshold = 1e-12, maxiter = 1000))
  mu <- as.numeric(stats::coef(model))
  own <- .bef_reml_normal(data$yi, sqrt(data$vi))
  stopifnot(abs(mu - own$mean) < 1e-6 * stats::sd(data$yi),
            abs(sqrt(model$tau2) - own$sd) < 1e-6 * stats::sd(data$yi))
  weight <- model$tau2 / (model$tau2 + data$vi)
  mean <- mu + weight * (data$yi - mu)
  sd <- sqrt(weight * data$vi)
  list(mean = mu, sd = sqrt(model$tau2), method = "REML plug-in normal random effects",
       sites = data.frame(site = data$site, mean = mean, sd = sd,
                          lower = mean - stats::qnorm(0.95) * sd,
                          upper = mean + stats::qnorm(0.95) * sd))
}
selected_of <- function(case, data) {
  if (identical(case, "krw")) {
    c(which.min(data$sei), order(data$sei)[[49L]], which.max(data$sei))
  } else {
    c(1L, 750L, 1500L)
  }
}
metrics_of <- function(fit, case, data) {
  grid <- fit$metadata$data_list$grid
  g <- g_summary(fit, type = "mass")$mean
  limits <- confint(fit)
  positive <- exceedance_prob(fit, 0)
  selected <- selected_of(case, data)
  truth <- data$theta_true
  data.frame(
    g_positive_mass = g_summary(fit, type = "prob", lower = 0)$mean,
    g_peak = grid[which.max(g)],
    observed_sd = stats::sd(data$yi), posterior_means_sd = stats::sd(coef(fit)),
    positive_intervals = sum(limits$lower > 0),
    selected_sites = paste(selected, collapse = ";"),
    selected_means = paste(signif(coef(fit)[selected], 7), collapse = ";"),
    selected_positive_probabilities = paste(signif(positive[selected], 7), collapse = ";"),
    rmse = if (is.null(truth)) NA_real_ else sqrt(mean((coef(fit) - truth)^2)),
    coverage90 = if (is.null(truth)) NA_real_ else mean(truth >= limits$lower & truth <= limits$upper),
    stringsAsFactors = FALSE
  )
}
curve_of <- function(fit, case, spec) {
  mass <- g_summary(fit, type = "mass")
  data.frame(case = case, spec = spec, grid = mass$value, probability = mass$mean,
             density = mass$mean / diff(mass$value[1:2]))
}
case_of <- function(stored, full, case, data) {
  normal <- normal_of(data)
  limits <- confint(stored)
  truth <- data$theta_true
  inside <- function(lower, upper) {
    if (is.null(truth)) NA_real_ else mean(truth >= lower & truth <= upper)
  }
  error <- function(x) if (is.null(truth)) NA_real_ else sqrt(mean((x - truth)^2))
  list(
    normal = normal,
    # With `reading = "step"` the probabilities are those of the grid points
    # above zero; `vignette("implementation")` compares them between the
    # stored draws and all the draws. The interpolated ones follow.
    positive_probability = unname(exceedance_prob(stored, 0, reading = "step")),
    full_positive_probability = unname(exceedance_prob(full, 0, reading = "step")),
    positive_probability_interpolated = unname(exceedance_prob(stored, 0)),
    compact_metrics = data.frame(
      observed_sd = stats::sd(data$yi), posterior_means_sd = stats::sd(coef(stored)),
      positive_intervals = sum(limits$lower > 0),
      rmse = error(coef(stored)), coverage90 = inside(limits$lower, limits$upper)
    ),
    normal_metrics = data.frame(
      rmse = error(normal$sites$mean),
      coverage90 = inside(normal$sites$lower, normal$sites$upper)
    )
  )
}

# ---- the callback gaps: three fits at full length ----------------------------
krw_specs <- list(
  baseline = list(),
  narrower = list(grid_method = "paper_realdata", expansion = 0.25),
  basis8 = list(M = 8L)
)
sensitivity <- list()
settings <- list()
curves <- list()
krw_full <- NULL
for (spec in names(krw_specs)) {
  fit <- do.call(quiet_fit, c(list(krw_discrimination, seed = 20250807L), krw_specs[[spec]]))
  if (identical(spec, "baseline")) {
    krw_full <- fit
  }
  sensitivity[[paste0("krw-", spec)]] <- cbind(
    data.frame(case = "krw", spec = spec, stringsAsFactors = FALSE),
    diagnostics_of(fit), metrics_of(fit, "krw", krw_discrimination)
  )
  settings[[paste0("krw-", spec)]] <- settings_of(fit, "krw", spec, krw_discrimination)
  curves[[paste0("krw-", spec)]] <- curve_of(fit, "krw", spec)
  message("fitted: krw ", spec)
}
# The stored fit holds one draw in fifteen of the baseline.
stopifnot(identical(
  unclass(krw_fit$draws)[, , "lambda"],
  unclass(posterior::thin_draws(krw_full$draws, 15L))[, , "lambda"]
))

# ---- the simulated sites -----------------------------------------------------
# The fits were made at full length by data-raw/example-fits.R; their rows
# are read from its table. The stored fit gives the curve and the case.
comparison <- readRDS(file.path("inst", "examples", "twin_towers_comparison.rds"))$table
twin_rows <- c(
  baseline = "example: rule \"pilot\", M = 8",
  default_M6 = "defaults: rule \"pilot\", M = 6",
  true_range = "range around the true effects, M = 6"
)
for (spec in names(twin_rows)) {
  row <- comparison[comparison$fit == twin_rows[[spec]], ]
  stopifnot(nrow(row) == 1L)
  sensitivity[[paste0("twin-", spec)]] <- data.frame(
    case = "twin", spec = spec, K = nrow(twin_towers),
    max_rhat = row$max_rhat, min_ess_bulk = row$min_ess_bulk, min_ess_tail = row$min_ess_tail,
    divergences = row$divergences, treedepth_hits = row$treedepth_hits,
    min_ebfmi = row$min_ebfmi, sampler_seconds = row$sampling_seconds,
    mean_g = row$mean_g, sd_g = row$sd_g,
    g_positive_mass = NA_real_, g_peak = NA_real_,
    observed_sd = stats::sd(twin_towers$yi), posterior_means_sd = NA_real_,
    positive_intervals = NA_integer_, selected_sites = NA_character_,
    selected_means = NA_character_, selected_positive_probabilities = NA_character_,
    rmse = row$rmse, coverage90 = row$coverage_90, stringsAsFactors = FALSE
  )
}
settings[["twin-baseline"]] <- settings_of(twin_towers_fit, "twin", "baseline", twin_towers)
curves[["twin-baseline"]] <- curve_of(twin_towers_fit, "twin", "baseline")
twin_full <- quiet_fit(twin_towers, M = 8L, seed = 20250627L)
stopifnot(identical(
  unclass(twin_towers_fit$draws)[, , "lambda"],
  unclass(posterior::thin_draws(twin_full$draws, 15L))[, , "lambda"]
))

result <- list(
  sensitivity = do.call(rbind, sensitivity),
  settings = do.call(rbind, settings),
  curves = do.call(rbind, curves),
  cases = list(
    krw = case_of(krw_fit, krw_full, "krw", krw_discrimination),
    twin = case_of(twin_towers_fit, twin_full, "twin", twin_towers)
  ),
  about = list(
    analysis_date = format(Sys.Date()), interval_level = 0.9,
    package_version = as.character(utils::packageVersion("bayesEfron")),
    cmdstan_version = krw_fit$metadata$cmdstan_version,
    fitting_script = "data-raw/example-fits.R",
    evidence_script = "data-raw/example-evidence.R",
    warning = "One simulated panel; the normal intervals leave out the uncertainty of the mean and the variance."
  )
)
rownames(result$sensitivity) <- rownames(result$settings) <- rownames(result$curves) <- NULL
saveRDS(result, file.path("inst", "examples", "lee_sui_examples.rds"), compress = "xz",
        version = 3)
utils::write.csv(result$sensitivity, file.path(audit_dir, "example-evidence-sensitivity.csv"),
                 row.names = FALSE)
options(width = 200)
print(result$sensitivity[, c("case", "spec", "max_rhat", "min_ess_bulk", "divergences",
                             "treedepth_hits", "mean_g", "sd_g", "positive_intervals", "rmse",
                             "coverage90")], digits = 4, row.names = FALSE)
message(sprintf("inst/examples/lee_sui_examples.rds: %d bytes",
                file.size(file.path("inst", "examples", "lee_sui_examples.rds"))))
