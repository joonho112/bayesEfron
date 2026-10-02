# Builds the two data sets in data/: `raudenbush1985`, the teacher expectancy
# studies, and `raudenbush_fit`, the model fitted to them with the default
# settings. Needs CmdStan. Run from the package root:
#
#   Rscript data-raw/raudenbush.R
#
# The fit takes a few seconds once the model is compiled. To keep the package
# small, one draw in fifteen is kept (800 of 12,000). The convergence
# diagnostics stored in the object are those of the full run.

devtools::load_all(".", quiet = TRUE)

# ---- the data ----------------------------------------------------------------
# Raudenbush and Bryk (1985, Table 1); `setting` and `tester` from Raudenbush
# (1984). The same data are distributed as `dat.raudenbush1985` in the metadat
# package, which is where the values were checked.
raudenbush1985 <- data.frame(
  study = 1:19,
  author = c(
    "Rosenthal et al.", "Conn et al.", "Jose & Cody", "Pellegrini & Hicks",
    "Pellegrini & Hicks", "Evans & Rosenthal", "Fielder et al.", "Claiborn",
    "Kester", "Maxwell", "Carter", "Flowers", "Keshock", "Henrikson", "Fine",
    "Grieger", "Rosenthal & Jacobson", "Fleming & Anttonen", "Ginsburg"
  ),
  year = c(
    1974L, 1968L, 1971L, 1972L, 1972L, 1969L, 1971L, 1969L, 1969L, 1970L,
    1970L, 1966L, 1970L, 1970L, 1972L, 1970L, 1968L, 1971L, 1970L
  ),
  weeks = c(2L, 21L, 19L, 0L, 0L, 3L, 17L, 24L, 0L, 1L, 0L, 0L, 1L, 2L, 17L, 5L, 1L, 2L, 7L),
  setting = c(
    "group", "group", "group", "group", "group", "group", "group", "group",
    "group", "indiv", "group", "group", "indiv", "indiv", "group", "group",
    "group", "group", "group"
  ),
  tester = c(
    "aware", "aware", "aware", "aware", "blind", "aware", "blind", "aware",
    "aware", "blind", "blind", "blind", "blind", "blind", "aware", "blind",
    "aware", "blind", "aware"
  ),
  n1i = c(77L, 60L, 72L, 11L, 11L, 129L, 110L, 26L, 75L, 32L, 22L, 43L, 24L, 19L, 80L, 72L, 65L, 233L, 65L),
  n2i = c(339L, 198L, 72L, 22L, 22L, 348L, 636L, 99L, 74L, 32L, 22L, 38L, 24L, 32L, 79L, 72L, 255L, 224L, 67L),
  yi = c(
    0.03, 0.12, -0.14, 1.18, 0.26, -0.06, -0.02, -0.32, 0.27, 0.80,
    0.54, 0.18, -0.02, 0.23, -0.18, -0.06, 0.30, 0.07, -0.07
  ),
  vi = c(
    0.0156, 0.0216, 0.0279, 0.1391, 0.1362, 0.0106, 0.0106, 0.0484, 0.0269,
    0.0630, 0.0912, 0.0497, 0.0835, 0.0841, 0.0253, 0.0279, 0.0193, 0.0088,
    0.0303
  ),
  stringsAsFactors = FALSE
)
raudenbush1985$sei <- sqrt(raudenbush1985$vi)

if (requireNamespace("metadat", quietly = TRUE)) {
  reference <- metadat::dat.raudenbush1985
  shared <- setdiff(names(raudenbush1985), "sei")
  for (column in shared) {
    stopifnot(isTRUE(all.equal(
      raudenbush1985[[column]], as.vector(reference[[column]]),
      check.attributes = FALSE
    )))
  }
  message("raudenbush1985: ", length(shared), " columns equal to metadat::dat.raudenbush1985")
}

# ---- the fit -----------------------------------------------------------------
full <- bayes_efron_fit(
  theta_hat = raudenbush1985$yi,
  sigma = raudenbush1985$sei,
  seed = 1985L,
  keep_cmdstan_fit = TRUE
)
# The run must have converged: no failed check, every Rhat below 1.01, every
# bulk effective sample size above 400 and no divergent transition. A few
# transitions that reach the maximum tree depth are common with these data
# and are reported by diagnose(); they are left as they are.
full_diagnostic <- diagnose(full)
print(summary(full_diagnostic))
stopifnot(
  length(full_diagnostic$sampler_diagnostics_failed) == 0L,
  length(full_diagnostic$diagnostic_skipped) == 0L,
  max(full_diagnostic$rhat) < 1.01,
  min(full_diagnostic$ess_bulk) > 400,
  full_diagnostic$divergences == 0
)

# One draw in fifteen. The summaries are computed again from the kept draws,
# so that they agree with what is stored; the diagnostics are carried over
# from the full run.
kept <- posterior::thin_draws(full$draws, thin = 15L)
stand_in <- list(
  draws = function(format = "draws_array", ...) kept,
  diagnostic_summary = function(...) full$cmdstan_fit$diagnostic_summary(...)
)
processed <- postprocess_stan_draws(
  cmdstan_fit = stand_in,
  stan_data = full$metadata$data_list,
  model_family = "RE"
)
metadata <- full$metadata
for (field in c(
  "mean_g_summary", "var_g_summary", "theta_summary",
  "effective_params_summary", "log_marginal_likelihood_summary"
)) {
  metadata[[field]] <- processed[[field]]
}
attr(metadata, "sd_g_summary") <- processed$sd_g_summary
raudenbush_fit <- validate_bef_fit_re(new_bef_fit_re(
  draws = processed$draws,
  metadata = metadata,
  posterior = processed$posterior
))
# The posterior means from the kept draws agree with those of the full run to
# within Monte Carlo error.
mc_error <- apply(.bef_site_draws(raudenbush_fit, "theta_mean"), 2L, stats::sd) / sqrt(800)
shift <- abs(coef(raudenbush_fit) - coef(full))
message(sprintf(
  "largest change in a posterior mean after thinning: %.4f (%.1f Monte Carlo standard errors)",
  max(shift), max(shift / mc_error)
))
stopifnot(
  identical(attr(raudenbush_fit$metadata, "diagnostics"), attr(full$metadata, "diagnostics")),
  posterior::ndraws(raudenbush_fit$draws) == 800L,
  all(shift < 4 * mc_error)
)

# ---- save --------------------------------------------------------------------
dir.create("data", showWarnings = FALSE)
save(raudenbush1985, file = "data/raudenbush1985.rda", compress = "xz", version = 3)
save(raudenbush_fit, file = "data/raudenbush_fit.rda", compress = "xz", version = 3)
message(sprintf(
  "data/raudenbush1985.rda %d bytes; data/raudenbush_fit.rda %d bytes (%s in memory)",
  file.size("data/raudenbush1985.rda"), file.size("data/raudenbush_fit.rda"),
  format(utils::object.size(raudenbush_fit), units = "MB")
))
message(sprintf(
  "CmdStan %s, cmdstanr %s, %s, seed %d, runtime %.1f s",
  raudenbush_fit$metadata$cmdstan_version, utils::packageVersion("cmdstanr"),
  format(Sys.Date()), raudenbush_fit$metadata$seed, raudenbush_fit$metadata$runtime_seconds
))
