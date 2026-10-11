# Builds the two files that the article "Small classes in 79 schools" reads,
# in vignettes/articles/:
#
#   class-size-fit.rds       the fit of the 79 schools of `star_schools`:
#                            bayes_efron_fit(star_schools, seed = 1985,
#                              adapt_delta = 0.99, max_treedepth = 12,
#                              iter_warmup = 2000, iter_sampling = 6000,
#                              thin = 30)
#   class-size-evidence.rds  the runs that the article sets beside it
#
# Needs CmdStan, and class-size-first-stage.rds, which data-raw/star.R
# writes. Run from the package root:
#
#   Rscript data-raw/class-size.R [audit-dir]
#
# Four tables are also written as CSV files to `audit-dir`, or to a
# temporary directory if none is given: the runs, the stored draws beside
# all the draws, the normal model, and the fits under the other assumptions
# about the variance of the classes.
#
# The fit runs four chains of 2,000 warmup and 6,000 sampling iterations and
# keeps one draw in thirty, 800 in all. The same call is also run with all
# its draws kept, and the script stops if the summaries of the 800 draws are
# far from those of the 24,000, or if the run misses the targets below.
# Unlike the fits in data/, this one keeps a warning: the energy diagnostic
# E-BFMI is below 0.3 in each of the twelve runs of `runs` below, and the
# article reports it.
#
# The evidence file holds
#
#   runs         the diagnostics and the summaries of g for six seeds under
#                the defaults and under the settings of the stored fit;
#   full_sites   the summaries of the schools from all the draws of the
#                stored call;
#   ols          the fit of the same schools with the standard errors that
#                take the students of a school as independent (`yi_ols`,
#                `sei_ols`);
#   normal       the normal model fitted by REML to both sets of estimates;
#   first_stage  four rows: the stored call with all its draws, and the same
#                call with the contrasts and standard errors under three
#                other assumptions about the variance of the classes (that
#                variance at the two ends of its 95% interval, and one
#                variance for each kind of class).

args <- commandArgs(trailingOnly = TRUE)
audit_dir <- if (length(args) >= 1L) args[[1L]] else file.path(tempdir(), "bayesEfron-class-size")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
suppressMessages(devtools::load_all(".", quiet = TRUE))
out_dir <- file.path("vignettes", "articles")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

quiet_fit <- function(...) {
  withCallingHandlers(
    bayes_efron_fit(...),
    warning = function(w) {
      if (inherits(w, "bef_sampler_diagnostics_failed")) invokeRestart("muffleWarning")
    },
    message = function(m) {
      if (inherits(m, "bef_default_changed")) invokeRestart("muffleMessage")
    }
  )
}
level <- 0.9
run_row <- function(fit, setting, seed) {
  d <- diagnose(fit)
  g <- confint(fit, type = "g", level = level)
  rownames(g) <- g$parameter
  sampler <- attr(fit$metadata, "sampler_settings", exact = TRUE)
  data.frame(
    setting = setting, seed = seed, chains = sampler$chains,
    iter_warmup = sampler$iter_warmup, iter_sampling = sampler$iter_sampling,
    adapt_delta = sampler$adapt_delta, max_treedepth = sampler$max_treedepth,
    ebfmi_min = min(d$ebfmi), ebfmi_max = max(d$ebfmi),
    chains_below_0.3 = sum(d$ebfmi < 0.3),
    max_rhat = max(d$rhat, na.rm = TRUE),
    min_ess_bulk = min(d$ess_bulk, na.rm = TRUE), min_ess_tail = min(d$ess_tail, na.rm = TRUE),
    divergences = sum(d$divergences), treedepth_hits = sum(d$max_treedepth),
    failed = paste(d$sampler_diagnostics_failed, collapse = ", "),
    warned = paste(d$sampler_diagnostics_warned, collapse = ", "),
    mean_g = g["mean_g", "mean"], mean_g_lower = g["mean_g", "lower"],
    mean_g_upper = g["mean_g", "upper"],
    sd_g = g["sd_g", "mean"], sd_g_lower = g["sd_g", "lower"], sd_g_upper = g["sd_g", "upper"],
    sampling_seconds = fit$metadata$runtime_seconds,
    stringsAsFactors = FALSE
  )
}

# ---- the stored fit ------------------------------------------------------------
stored_settings <- list(adapt_delta = 0.99, max_treedepth = 12L,
                        iter_warmup = 2000L, iter_sampling = 6000L)
stored_seed <- 1985L
full <- do.call(quiet_fit, c(list(star_schools, seed = stored_seed), stored_settings))
fit <- do.call(quiet_fit, c(list(star_schools, seed = stored_seed), stored_settings,
                            list(thin = 30L)))
d <- diagnose(fit)
stopifnot(
  identical(dim(fit$draws)[1:2], c(200L, 4L)),
  identical(unclass(fit$draws)[, , "lambda"],
            unclass(posterior::thin_draws(full$draws, 30L))[, , "lambda"]),
  identical(attr(fit$metadata, "diagnostics", exact = TRUE),
            attr(full$metadata, "diagnostics", exact = TRUE)),
  is.null(fit$cmdstan_fit),
  # The targets: the warning about E-BFMI and no other, no failed check.
  identical(d$sampler_diagnostics_failed, character()),
  identical(d$sampler_diagnostics_warned, "ebfmi"),
  max(d$rhat, na.rm = TRUE) < 1.01,
  min(d$ess_bulk, na.rm = TRUE) >= 400, min(d$ess_tail, na.rm = TRUE) >= 400,
  !any(d$grid$ends$status == "warning")
)
# The summaries of the 800 draws beside those of the 24,000.
spacing <- diff(fit$metadata$data_list$grid[1:2])
a <- site_summary(full, level = level)
b <- site_summary(fit, level = level)
limits <- abs(as.matrix(b[, c("lower", "upper")]) - as.matrix(a[, c("lower", "upper")])) / spacing
means <- abs(b$mean - a$mean) / b$se
thinning <- data.frame(
  largest_limit_difference_in_spacings = max(limits),
  largest_mean_difference_in_standard_errors = max(means),
  mean_g_difference = mean(fit$posterior$mean_g) - mean(full$posterior$mean_g),
  sd_g_difference = mean(fit$posterior$sd_g) - mean(full$posterior$sd_g)
)
stopifnot(max(limits) < 2, max(means) < 0.05)
saveRDS(fit, file.path(out_dir, "class-size-fit.rds"), compress = "xz")

# ---- the same data under other seeds and the default settings -------------------
seeds <- c(stored_seed, 1:5)
runs <- list(run_row(full, "stored", stored_seed))
for (seed in seeds[-1L]) {
  runs[[length(runs) + 1L]] <- run_row(
    do.call(quiet_fit, c(list(star_schools, seed = seed), stored_settings)), "stored", seed
  )
}
default_fits <- list()
for (seed in seeds) {
  default_fits[[as.character(seed)]] <- quiet_fit(star_schools, seed = seed)
  runs[[length(runs) + 1L]] <- run_row(default_fits[[as.character(seed)]], "defaults", seed)
}
runs <- do.call(rbind, runs)
runs <- runs[order(runs$setting != "defaults", match(runs$seed, seeds)), ]
rownames(runs) <- NULL

# ---- the standard errors that ignore the classes --------------------------------
ols_data <- as_bef_data(star_schools, theta_hat = "yi_ols", sigma = "sei_ols")
ols_fit <- quiet_fit(ols_data, seed = stored_seed)
ols_sites <- site_summary(ols_fit, level = level)
ols <- list(
  call = "bayes_efron_fit(as_bef_data(star_schools, theta_hat = \"yi_ols\", sigma = \"sei_ols\"), seed = 1985)",
  run = run_row(ols_fit, "defaults, OLS standard errors", stored_seed),
  g = confint(ols_fit, type = "g", level = level),
  density = g_summary(ols_fit, type = "density", level = level),
  sites = ols_sites[, c("site", "estimate", "se", "mean", "sd", "lower", "upper")],
  intervals_excluding_zero = sum(ols_sites$lower > 0 | ols_sites$upper < 0)
)

# ---- the normal model -----------------------------------------------------------
reml_row <- function(y, s, label) {
  estimate <- .bef_reml_normal(y, s)
  w <- 1 / s^2
  q <- sum(w * (y - sum(w * y) / sum(w))^2)
  data.frame(standard_errors = label, mean = estimate$mean, sd = estimate$sd,
             q = q, df = length(y) - 1L,
             p = stats::pchisq(q, length(y) - 1L, lower.tail = FALSE),
             stringsAsFactors = FALSE)
}
normal <- rbind(
  reml_row(star_schools$yi, star_schools$sei, "classes"),
  reml_row(star_schools$yi_ols, star_schools$sei_ols, "OLS")
)

# ---- other assumptions about the variance of the classes ------------------------
# The contrasts and standard errors of data-raw/star.R under each assumption,
# fitted with the call of the stored fit. The first row is the stored call
# itself, with all its draws.
first_stage <- readRDS(file.path(out_dir, "class-size-first-stage.rds"))
variant_names <- c("common", "lower", "upper", "by_kind")
stopifnot(identical(unique(first_stage$variants$variant), variant_names))
variant_rows <- lapply(variant_names, function(name) {
  v <- first_stage$variants[first_stage$variants$variant == name, ]
  stopifnot(identical(v$site, star_schools$site))
  if (identical(name, "common")) {
    stopifnot(max(abs(v$yi - star_schools$yi)) < 1e-10, max(abs(v$sei - star_schools$sei)) < 1e-10)
    fitted <- full
  } else {
    fitted <- do.call(quiet_fit, c(list(v$yi, v$sei, seed = stored_seed), stored_settings))
  }
  sites <- site_summary(fitted, level = level)
  reml <- reml_row(v$yi, v$sei, name)
  run <- run_row(fitted, name, stored_seed)
  data.frame(
    variant = name,
    se_ratio_median = stats::median(v$sei / star_schools$sei),
    largest_estimate_change = max(abs(v$yi - star_schools$yi)),
    reml_sd = reml$sd, q = reml$q, p = reml$p,
    mean_g = run$mean_g, sd_g = run$sd_g, sd_g_lower = run$sd_g_lower,
    sd_g_upper = run$sd_g_upper,
    ebfmi_min = run$ebfmi_min, max_rhat = run$max_rhat, min_ess_bulk = run$min_ess_bulk,
    divergences = run$divergences, treedepth_hits = run$treedepth_hits,
    failed = run$failed, warned = run$warned,
    intervals_excluding_zero = sum(sites$lower > 0 | sites$upper < 0),
    site_5_lower = sites$lower[5L], site_5_upper = sites$upper[5L],
    site_60_lower = sites$lower[60L], site_60_upper = sites$upper[60L],
    stringsAsFactors = FALSE
  )
})
first_stage_fits <- do.call(rbind, variant_rows)
full_sites <- site_summary(full, level = level)[, c("site", "estimate", "se", "mean", "sd",
                                                    "lower", "upper")]

evidence <- list(
  runs = runs, full_sites = full_sites, ols = ols, normal = normal,
  first_stage = first_stage_fits, thinning = thinning, level = level,
  stored_call = paste0("bayes_efron_fit(star_schools, seed = 1985, adapt_delta = 0.99, ",
                       "max_treedepth = 12, iter_warmup = 2000, iter_sampling = 6000, thin = 30)"),
  seeds = seeds,
  created = format(Sys.Date()),
  bayesEfron_version = as.character(utils::packageVersion("bayesEfron")),
  cmdstan_version = fit$metadata$cmdstan_version
)
saveRDS(evidence, file.path(out_dir, "class-size-evidence.rds"), compress = "xz")

utils::write.csv(runs, file.path(audit_dir, "class-size-runs.csv"), row.names = FALSE)
utils::write.csv(thinning, file.path(audit_dir, "class-size-thinning.csv"), row.names = FALSE)
utils::write.csv(normal, file.path(audit_dir, "class-size-normal.csv"), row.names = FALSE)
utils::write.csv(first_stage_fits, file.path(audit_dir, "class-size-first-stage.csv"),
                 row.names = FALSE)
options(width = 200)
print(runs[, c("setting", "seed", "ebfmi_min", "chains_below_0.3", "max_rhat", "min_ess_bulk",
               "min_ess_tail", "divergences", "treedepth_hits", "mean_g", "sd_g",
               "sd_g_lower", "sd_g_upper")], digits = 3)
print(ols$run[, c("setting", "ebfmi_min", "max_rhat", "divergences", "mean_g", "sd_g",
                  "sd_g_lower", "sd_g_upper")], digits = 3)
print(normal, digits = 4)
print(first_stage_fits[, c("variant", "se_ratio_median", "reml_sd", "p", "sd_g", "sd_g_lower",
                           "sd_g_upper", "ebfmi_min", "divergences", "warned",
                           "intervals_excluding_zero", "site_60_lower", "site_60_upper")],
      digits = 3)
print(thinning, digits = 3)
for (file in c("class-size-fit.rds", "class-size-evidence.rds")) {
  cat(sprintf("%s %d bytes\n", file.path(out_dir, file), file.size(file.path(out_dir, file))))
}
