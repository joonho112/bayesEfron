# ==============================================================================
# analysis.R: the fits and registered numerical results for the SoftwareX article on
# bayesEfron 0.4.0.
#
# Run from the paper folder root:   Rscript scripts/analysis.R
#
# Needs CmdStan (2.38.0 here) and the packages bayesEfron 0.4.0 (installed in
# r-library/ by scripts/install-package.sh), metafor, ebnm and loo. It writes
#
#   data/fits.rds            the fits of the two examples and of the 19 studies
#                            (about 50 MB; not needed by the figure script)
#   data/krw.rds             tables and curves for Figure 3 and Section 3.1
#   data/twin.rds            tables and curves for Figure 4 and Section 3.2
#   data/evidence.rds        the simulation tables shipped with the package
#   data/grid-comparison.rds compare_efron_grids() of the 19 studies (Figure 5)
#   data/numbers.csv         the registry of selected numerical quotations
#   data/console/*.txt       the console output printed in the article
#
# The simulation tables in inst/verification/ and inst/examples/ of the
# package are read, not recomputed (the fits behind them take hours).
# ==============================================================================

root <- normalizePath(".")
stopifnot(file.exists(file.path(root, "scripts", "analysis.R")))
.libPaths(c(file.path(root, "r-library"), .libPaths()))
suppressPackageStartupMessages({
  library(bayesEfron)
  library(metafor)
  library(ebnm)
  library(loo)
})
print(packageVersion("bayesEfron"))
stopifnot(packageVersion("bayesEfron") == "0.4.0")
options(bayesEfron.use_cli = FALSE, width = 76)
dir.create("data/console", showWarnings = FALSE, recursive = TRUE)

registry <- list()
reg <- function(name, value, digits = NA) {
  registry[[name]] <<- list(name = name, value = value, digits = digits)
  invisible(value)
}
console <- function(name, expr) {
  text <- utils::capture.output(expr)
  writeLines(text, file.path("data", "console", paste0(name, ".txt")))
  invisible(text)
}
pkg_file <- function(...) system.file(..., package = "bayesEfron")
started_all <- Sys.time()

# ---- 1. The callback gaps of 97 firms ---------------------------------------
message("fit: 97 firms")
dat <- as_bef_data(krw_discrimination, site = "firm_name")
krw_timer <- system.time(
  fit_krw <- bayes_efron_fit(dat, seed = 20250807)
)
grid_krw <- attr(fit_krw$metadata, "grid")
checks_krw <- as.data.frame(diagnose(fit_krw))
stopifnot(all(checks_krw$status == "ok"), !grid_krw$pilot$narrowed,
          grid_krw$L == 101L, grid_krw$M == 6L)

# Console output for Section 3.1
console("krw-data", print(head(krw_discrimination[, c("firm_name", "yi", "sei", "n_jobs")], 3)))
console("krw-fit", print(fit_krw))
console("krw-diagnose", print(diagnose(fit_krw)))
console("krw-g-moments", print(g_summary(fit_krw, type = "moments", scale = 100)))
console("krw-g-prob", print(g_summary(fit_krw, type = "prob", lower = 0.02)))
sites_krw <- site_summary(fit_krw, threshold = 0, scale = 100)
top3 <- sites_krw[order(-sites_krw$mean)[1:3],
                  c("label", "estimate", "se", "mean", "lower", "upper", "prob_above")]
rownames(top3) <- NULL
console("krw-sites", print(top3, digits = 3))
p_pos <- exceedance_prob(fit_krw)
console("krw-exceedance", print(sum(p_pos >= 0.95)))
console("krw-predict", print(predict(fit_krw, newdata = data.frame(estimate = c(0.05, 0.05),
                                                                  se = c(0.01, 0.03)))))

# The normal random-effects model fitted by REML to the same estimates.
reml_krw <- metafor::rma(yi = krw_discrimination$yi, vi = krw_discrimination$vi,
                         method = "REML", control = list(threshold = 1e-12, maxiter = 1000))
normal_krw <- list(mean = as.numeric(coef(reml_krw)), sd = sqrt(reml_krw$tau2))

g_ci <- confint(fit_krw, type = "g")
rownames(g_ci) <- g_ci$parameter
above_two <- g_summary(fit_krw, type = "prob", lower = 0.02)
below_zero <- g_summary(fit_krw, type = "prob", upper = 0)
above_five <- g_summary(fit_krw, type = "prob", lower = 0.05)
share_about_g <- sites_krw$sd_about_g^2 / sites_krw$sd^2

reg("krw_K", nrow(krw_discrimination))
reg("krw_applications", sum(krw_discrimination$n_obs))
reg("krw_vacancies", sum(krw_discrimination$n_jobs))
reg("krw_se_ratio", max(krw_discrimination$sei) / min(krw_discrimination$sei), 1)
reg("krw_seconds_total", as.numeric(krw_timer[["elapsed"]]), 0)
reg("krw_seconds_sampling", fit_krw$metadata$runtime_seconds, 1)
reg("krw_seconds_pilot", grid_krw$pilot$seconds, 1)
reg("krw_L", grid_krw$L); reg("krw_M", grid_krw$M)
reg("krw_grid_lower_pp", 100 * grid_krw$support[1], 2)
reg("krw_grid_upper_pp", 100 * grid_krw$support[2], 2)
reg("krw_mean_g_pp", 100 * g_ci["mean_g", "mean"], 2)
reg("krw_mean_g_lower_pp", 100 * g_ci["mean_g", "lower"], 2)
reg("krw_mean_g_upper_pp", 100 * g_ci["mean_g", "upper"], 2)
reg("krw_sd_g_pp", 100 * g_ci["sd_g", "mean"], 2)
reg("krw_sd_g_lower_pp", 100 * g_ci["sd_g", "lower"], 2)
reg("krw_sd_g_upper_pp", 100 * g_ci["sd_g", "upper"], 2)
reg("krw_sd_estimates_pp", 100 * sd(krw_discrimination$yi), 2)
reg("krw_sd_posterior_means_pp", 100 * sd(sites_krw$mean / 100), 2)
reg("krw_reml_mean_pp", 100 * normal_krw$mean, 2)
reg("krw_reml_sd_pp", 100 * normal_krw$sd, 2)
reg("krw_prob_above_2", above_two$mean, 2)
reg("krw_prob_above_2_lower", above_two$lower, 2)
reg("krw_prob_above_2_upper", above_two$upper, 2)
reg("krw_prob_below_0", below_zero$mean, 2)
reg("krw_prob_above_5", above_five$mean, 3)
reg("krw_prob_above_5_lower", above_five$lower, 3)
reg("krw_prob_above_5_upper", above_five$upper, 3)
reg("krw_normal_above_5", 1 - pnorm(0.05, normal_krw$mean, normal_krw$sd), 3)
reg("krw_n_prob_ge_95", sum(p_pos >= 0.95))
reg("krw_n_interval_above_0", sum(sites_krw$lower > 0))
reg("krw_top_firm", top3$label[1])
reg("krw_top_estimate_pp", top3$estimate[1], 2)
reg("krw_top_mean_pp", top3$mean[1], 2)
reg("krw_share_about_g_median_pct", 100 * median(share_about_g), 1)
reg("krw_share_about_g_max_pct", 100 * max(share_about_g), 0)
reg("krw_share_about_g_max_firm", sites_krw$label[which.max(share_about_g)])
reg("krw_rhat_max", max(checks_krw$value[checks_krw$diagnostic == "rhat"]), 3)
reg("krw_divergences", sum(diagnose(fit_krw)$divergences))

# Empirical Bayes fits of the same data with ebnm, for the comparison of the
# estimates of g (Figure 2a inset / OSM) and of the intervals. The protocol
# is that of data-raw/ebnm-protocol.md in the package: mode estimated,
# 90% equal-tailed intervals from 4,000 posterior draws.
message("ebnm: 97 firms")
everything <- ebnm::ebnm_output_all()
quietly <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("local false sign rates", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}
fit_ebnm <- function(fun, x, s) {
  fit <- quietly(fun(x, s))
  limits <- stats::quantile(fit, probs = c(0.05, 0.95), nsim = 4000L, seed = 1L)
  list(mean = unname(coef(fit)), sd = fitted(fit)$sd, lower = unname(limits[, 1]),
       upper = unname(limits[, 2]), loglik = as.numeric(logLik(fit)), g = fit$fitted_g)
}
cdf_of <- function(g, x) {
  if (inherits(g, "normalmix")) {
    smooth <- g$sd > 0
    return(vapply(x, function(v) {
      sum(g$pi[smooth] * pnorm(v, g$mean[smooth], g$sd[smooth])) +
        sum(g$pi[!smooth] * (v >= g$mean[!smooth]))
    }, numeric(1)))
  }
  if (inherits(g, "unimix")) {
    wide <- g$b > g$a
    return(vapply(x, function(v) {
      sum(g$pi[wide] * punif(v, g$a[wide], g$b[wide])) + sum(g$pi[!wide] * (v >= g$a[!wide]))
    }, numeric(1)))
  }
  rep(NA_real_, length(x))
}
eb_methods <- list(
  Normal = function(x, s) ebnm_normal(x, s, mode = "estimate", output = everything),
  Unimodal = function(x, s) ebnm_unimodal(x, s, mode = "estimate", output = everything),
  NPMLE = function(x, s) ebnm_npmle(x, s, output = everything)
)
compare_with_ebnm <- function(fit, x, s, truth = NULL) {
  sites <- site_summary(fit)
  rows <- lapply(eb_methods, fit_ebnm, x = x, s = s)
  g_hat <- as_normalmix(fit)
  plugged <- fit_ebnm(function(x, s) ebnm_npmle(x, s, g_init = g_hat, fix_g = TRUE,
                                                output = everything), x, s)
  # The plug-in intervals are computed from the same g in the way the package
  # computes its own, so that the two rows differ by the uncertainty about g
  # alone: a one-draw engine with the posterior mean of g.
  engine <- bayesEfron:::.bef_engine_from_parts(
    log_g = matrix(log(g_hat$pi / sum(g_hat$pi)), nrow = 1L),
    grid = fit$metadata$data_list$grid, theta_hat = x, sigma = s, tau_c = 0)
  weights <- bayesEfron:::.bef_site_weights(engine)
  limits <- bayesEfron:::.bef_site_quantile(engine, c(0.05, 0.95), weights = weights)
  plugged$lower <- limits[, 1]; plugged$upper <- limits[, 2]
  rows[["Plug-in g"]] <- plugged
  rows[["bayesEfron"]] <- list(mean = sites$mean, sd = sites$sd, lower = sites$lower,
                               upper = sites$upper, loglik = NA_real_, g = NULL)
  table <- do.call(rbind, lapply(names(rows), function(m) {
    r <- rows[[m]]
    out <- data.frame(
      method = m, correlation = cor(r$mean, sites$mean),
      largest_difference_se = max(abs(r$mean - sites$mean) / s),
      mean_sd = mean(r$sd), mean_length = mean(r$upper - r$lower),
      length_ratio = mean(r$upper - r$lower) / mean(sites$upper - sites$lower),
      loglik = r$loglik, stringsAsFactors = FALSE)
    if (!is.null(truth)) {
      inside <- truth >= r$lower & truth <= r$upper
      out$rmse <- sqrt(mean((r$mean - truth)^2))
      out$inclusion <- mean(inside)
      out$inclusion_se <- sqrt(mean(inside) * (1 - mean(inside)) / length(truth))
    }
    out
  }))
  grid <- fit$metadata$data_list$grid
  at <- seq(grid[1], grid[length(grid)], length.out = 401)
  band <- g_summary(fit, type = "cdf", at = at)
  cdf <- rbind(
    do.call(rbind, lapply(names(eb_methods), function(m) {
      data.frame(method = m, x = at, cdf = cdf_of(rows[[m]]$g, at), lower = NA, upper = NA)
    })),
    data.frame(method = "bayesEfron", x = at, cdf = band$mean, lower = band$lower,
               upper = band$upper)
  )
  list(table = table, cdf = cdf, rows = rows)
}
eb_krw <- compare_with_ebnm(fit_krw, krw_discrimination$yi, krw_discrimination$sei)
kt <- eb_krw$table; rownames(kt) <- kt$method
reg("krw_cor_npmle", kt["NPMLE", "correlation"], 3)
reg("krw_cor_normal", kt["Normal", "correlation"], 3)
reg("krw_length_ratio_npmle", kt["NPMLE", "length_ratio"], 2)
reg("krw_length_ratio_normal", kt["Normal", "length_ratio"], 2)
reg("krw_length_ratio_plugin", kt["Plug-in g", "length_ratio"], 3)
reg("krw_plugin_shorter_pct", 100 * (1 - kt["Plug-in g", "length_ratio"]), 1)

# The curves of Figure 2
dens_krw <- g_summary(fit_krw, type = "density")
sites50_krw <- site_summary(fit_krw, level = 0.9, inner_level = 0.5, scale = 100)
krw <- list(
  data = krw_discrimination, sites = sites_krw, sites50 = sites50_krw, p_pos = p_pos, g_ci = g_ci,
  density = dens_krw, normal = normal_krw, grid = grid_krw, checks = checks_krw,
  ebnm = eb_krw[c("table", "cdf")], share_about_g = share_about_g,
  loo = tryCatch(loo(fit_krw), error = function(e) NULL)
)
saveRDS(krw, "data/krw.rds")

# ---- 2. Twin Towers: 1,500 simulated sites with known effects ---------------
message("fit: 1,500 simulated sites")
sim <- twin_towers
twin_timer <- system.time(
  fit_twin <- bayes_efron_fit(twin_towers, M = 8, seed = 20250627)
)
grid_twin <- attr(fit_twin$metadata, "grid")
diag_twin <- diagnose(fit_twin)
checks_twin <- as.data.frame(diag_twin)
stopifnot(grid_twin$pilot$narrowed, grid_twin$M == 8L,
          all(checks_twin$status[checks_twin$diagnostic != "grid_ends"] == "ok"))
console("twin-fit-call", cat('fit <- bayes_efron_fit(twin_towers, M = 8, seed = 20250627)\n'))
console("twin-fit", print(fit_twin))
console("twin-diagnose", print(diag_twin))
console("twin-g", {
  print(g_summary(fit_twin, type = "prob", lower = -0.7, upper = 0.7))
  print(g_summary(fit_twin, type = "moments"))
})
sites_twin <- site_summary(fit_twin)
inside <- sim$theta_true >= sites_twin$lower & sim$theta_true <= sites_twin$upper
rmse <- function(x) sqrt(mean((x - sim$theta_true)^2))
console("twin-recovery", {
  cat("mean(theta_true inside 90% interval):", round(mean(inside), 3), "\n")
  print(c(posterior_means = rmse(sites_twin$mean), estimates = rmse(sim$yi)))
})
reml_twin <- metafor::rma(yi = sim$yi, vi = sim$vi, method = "REML",
                          control = list(threshold = 1e-12, maxiter = 1000))
normal_twin <- list(mean = as.numeric(coef(reml_twin)), sd = sqrt(reml_twin$tau2))
w_normal <- reml_twin$tau2 / (reml_twin$tau2 + sim$vi)
normal_means <- normal_twin$mean + w_normal * (sim$yi - normal_twin$mean)
normal_sd <- sqrt(w_normal * sim$vi)
normal_inside <- sim$theta_true >= normal_means - qnorm(0.95) * normal_sd &
  sim$theta_true <= normal_means + qnorm(0.95) * normal_sd
between <- g_summary(fit_twin, type = "prob", lower = -0.7, upper = 0.7)
normal_between <- pnorm(0.7, normal_twin$mean, normal_twin$sd) -
  pnorm(-0.7, normal_twin$mean, normal_twin$sd)
g_ci_twin <- confint(fit_twin, type = "g"); rownames(g_ci_twin) <- g_ci_twin$parameter
# Calibration: the share of true effects inside the interval at each level.
levels <- seq(0.05, 0.95, by = 0.05)
calibration <- data.frame(level = levels, inside = vapply(levels, function(l) {
  tab <- site_summary(fit_twin, level = l)
  mean(sim$theta_true >= tab$lower & sim$theta_true <= tab$upper)
}, numeric(1)))
# The posterior distribution of an imprecise site between the groups and of
# the most precise site, on the grid (interpolated reading not needed here).
engine_twin <- bayesEfron:::.bef_engine(fit_twin)
site_posterior <- function(i) {
  w <- bayesEfron:::.bef_site_weights(engine_twin, sites = i)
  data.frame(site = i, value = engine_twin$grid, posterior = as.numeric(w) /
               diff(engine_twin$grid[1:2]))
}
between_groups <- c(268L, 959L)
most_precise <- which.min(sim$sei)
shown_sites <- c(between_groups, most_precise)
post_curves <- do.call(rbind, lapply(shown_sites, site_posterior))
three <- site_summary(fit_twin, sites = shown_sites)
ends <- diag_twin$grid$ends
leaning <- ends[which.max(ends$mass), ]

reg("twin_K", nrow(sim))
reg("twin_low_n", sum(sim$theta_true < 0))
reg("twin_high_n", sum(sim$theta_true > 0))
reg("twin_se_min", min(sim$sei), 2); reg("twin_se_max", max(sim$sei), 2)
reg("twin_se_ratio", max(sim$sei) / min(sim$sei), 0)
reg("twin_se_above_gap_pct", 100 * mean(sim$sei > 1.4), 0)
reg("twin_sd_estimates", sd(sim$yi), 2); reg("twin_sd_truth", sd(sim$theta_true), 2)
reg("twin_seconds_total", as.numeric(twin_timer[["elapsed"]]), 0)
reg("twin_seconds_sampling", fit_twin$metadata$runtime_seconds, 0)
reg("twin_seconds_pilot", grid_twin$pilot$seconds, 0)
reg("twin_grid_lower", grid_twin$support[1], 2); reg("twin_grid_upper", grid_twin$support[2], 2)
reg("twin_wide_lower", grid_twin$pilot$wide[1], 1); reg("twin_wide_upper", grid_twin$pilot$wide[2], 1)
reg("twin_share_of_wide_pct", 100 * grid_twin$pilot$share, 0)
reg("twin_L", grid_twin$L); reg("twin_M", grid_twin$M)
reg("twin_prob_between", between$mean, 2)
reg("twin_normal_between", normal_between, 2)
reg("twin_mean_g", g_ci_twin["mean_g", "mean"], 2)
reg("twin_sd_g", g_ci_twin["sd_g", "mean"], 2)
reg("twin_mean_truth", mean(sim$theta_true), 2); reg("twin_sd_truth2", sd(sim$theta_true), 2)
reg("twin_inside_90_pct", 100 * mean(inside), 1)
reg("twin_normal_inside_90_pct", 100 * mean(normal_inside), 1)
reg("twin_rmse_posterior", rmse(sites_twin$mean), 3)
reg("twin_rmse_estimates", rmse(sim$yi), 3)
reg("twin_rmse_normal", rmse(normal_means), 3)
reg("twin_leaning_site", leaning$label); reg("twin_leaning_mass_pct", 100 * leaning$mass, 1)
reg("twin_site_a", shown_sites[1]); reg("twin_site_b", shown_sites[2])
reg("twin_site_a_se", three$se[1], 2); reg("twin_site_b_se", three$se[2], 2)
reg("twin_site_a_mean", three$mean[1], 2); reg("twin_site_b_mean", three$mean[2], 2)
reg("twin_site_a_width", three$upper[1] - three$lower[1], 1)
reg("twin_rhat_max", max(checks_twin$value[checks_twin$diagnostic == "rhat"]), 3)
reg("twin_divergences", sum(diag_twin$divergences))
reg("twin_treedepth_hits", sum(diag_twin$max_treedepth))

message("ebnm: 1,500 simulated sites")
eb_twin <- compare_with_ebnm(fit_twin, sim$yi, sim$sei, truth = sim$theta_true)
tt <- eb_twin$table; rownames(tt) <- tt$method
reg("twin_eb_rmse_npmle", tt["NPMLE", "rmse"], 3)
reg("twin_eb_rmse_normal", tt["Normal", "rmse"], 3)
reg("twin_eb_inclusion_npmle_pct", 100 * tt["NPMLE", "inclusion"], 1)
reg("twin_eb_inclusion_normal_pct", 100 * tt["Normal", "inclusion"], 0)
reg("twin_eb_inclusion_plugin_pct", 100 * tt["Plug-in g", "inclusion"], 0)
reg("twin_eb_length_ratio_plugin", tt["Plug-in g", "length_ratio"], 3)

twin <- list(
  data = sim, sites = sites_twin, g_ci = g_ci_twin, density = g_summary(fit_twin, type = "density"),
  normal = normal_twin, normal_means = normal_means, grid = grid_twin, checks = checks_twin,
  calibration = calibration, posterior_curves = post_curves, three = three,
  ebnm = eb_twin[c("table", "cdf")], inside = inside
)
saveRDS(twin, "data/twin.rds")

# ---- 3. The 19 teacher expectancy studies: sensitivity to the range ----------
message("compare_efron_grids: 19 studies")
grid_timer <- system.time(
  comparison <- compare_efron_grids(raudenbush1985, L = 101, M = 6,
                                    expansion = c(0.25, 0.5, 1), sites = c(4, 8),
                                    seed = 1985, keep_fits = TRUE)
)
console("raudenbush-comparison", print(comparison))
cmp <- as.data.frame(comparison)
cmp <- cmp[order(cmp$expansion), ]
for (a in c("fits", "g", "sites")) attr(cmp, a) <- NULL   # keep the table small
reg("raud_site4_narrow", cmp$site_4[1], 2); reg("raud_site4_default", cmp$site_4[2], 2)
reg("raud_site4_wide", cmp$site_4[3], 2)
reg("raud_sd_g_min", min(cmp$sd_g), 2); reg("raud_sd_g_max", max(cmp$sd_g), 2)
reg("raud_seconds_total", as.numeric(grid_timer[["elapsed"]]), 0)
fits_raud <- attr(comparison, "fits")
attr(comparison, "fits") <- NULL   # the fits go to data/fits.rds; the comparison stays small
reml_raud <- metafor::rma(yi = raudenbush1985$yi, vi = raudenbush1985$sei^2, method = "REML",
                          control = list(threshold = 1e-12, maxiter = 1000))
reg("raud_reml_sd", sqrt(reml_raud$tau2), 2)
reg("raud_reml_mean", as.numeric(coef(reml_raud)), 2)
saveRDS(list(comparison = comparison, table = cmp, g = attr(comparison, "g"),
             data = raudenbush1985,
             reml = list(mean = as.numeric(coef(reml_raud)), sd = sqrt(reml_raud$tau2))),
        "data/grid-comparison.rds")

# ---- 4. The simulation evidence shipped with the package --------------------
message("evidence tables")
verification <- function(name) read.csv(pkg_file("verification", name), check.names = FALSE)
panels_methods <- verification("panels-methods.csv")
panels_convergence <- verification("panels-convergence.csv")
few <- verification("few-studies.csv")
stored <- readRDS(pkg_file("examples", "ebnm_comparison.rds"))
equal <- stored$equal_se$repeated
evidence <- list(panels_methods = panels_methods, panels_convergence = panels_convergence,
                 few = few, equal = equal, stored_krw = stored$krw$table,
                 stored_twin = stored$twin$table, about = stored$about)
saveRDS(evidence, "data/evidence.rds")
be <- panels_methods[panels_methods$method == "bayesEfron", ]
reg("cov_K100", be$inclusion[be$K == 100], 3); reg("cov_K100_se", be$inclusion_se[be$K == 100], 3)
reg("cov_K1500", be$inclusion[be$K == 1500], 3); reg("cov_K1500_se", be$inclusion_se[be$K == 1500], 3)
np <- panels_methods[panels_methods$method == "NPMLE", ]
reg("cov_npmle_K100", np$inclusion[np$K == 100], 2); reg("cov_npmle_K1500", np$inclusion[np$K == 1500], 2)
fb <- few[few$method == "bayesEfron", ]
reg("cov_few_min", min(fb$inclusion), 3); reg("cov_few_max", max(fb$inclusion), 3)
fn <- few[few$method == "normal", ]
reg("cov_few_normal_min", min(fn$inclusion), 2); reg("cov_few_normal_max", max(fn$inclusion), 2)
eq <- equal[equal$method == "bayesEfron" & equal$M == 6, ]
reg("cov_equal_M6", eq$inclusion, 3); reg("cov_equal_M6_se", eq$inclusion_se, 3)
reg("cov_equal_efron_M6", equal$inclusion[equal$method == "Efron's estimator" & equal$M == 6], 3)
reg("few_sd_factor_min", min(fb$sd_mean / fb$true_sd), 1)
reg("few_sd_factor_max", max(fb$sd_mean / fb$true_sd), 1)
plug <- panels_methods[panels_methods$method == "g of bayesEfron, plugged in", ]
reg("extra_length_K100_pct", 100 * (be$mean_length[be$K == 100] / plug$mean_length[plug$K == 100] - 1), 1)
reg("extra_length_K1500_pct", 100 * (be$mean_length[be$K == 1500] / plug$mean_length[plug$K == 1500] - 1), 2)
reg("panels_default_valid", panels_convergence$valid[panels_convergence$setting == "defaults"])
reg("panels_default_seconds", panels_convergence$sampling_seconds_median[panels_convergence$setting == "defaults"], 0)
reg("panels_v030_valid", panels_convergence$valid[grepl("0.3.0", panels_convergence$setting)])
reg("panels_v030_seconds", panels_convergence$sampling_seconds_median[grepl("0.3.0", panels_convergence$setting)], 0)

# ---- 5. Package facts ---------------------------------------------------------
desc <- read.dcf(pkg_file("DESCRIPTION"))
reg("pkg_version", desc[, "Version"])
reg("pkg_r_min", "4.1.0")
reg("n_exported", length(getNamespaceExports("bayesEfron")))
reg("n_plot_functions", sum(startsWith(getNamespaceExports("bayesEfron"), "plot_")))
reg("cmdstan_version", as.character(cmdstanr::cmdstan_version()))
reg("r_version", paste(R.version$major, R.version$minor, sep = "."))
reg("cmdstanr_version", as.character(packageVersion("cmdstanr")))
reg("ebnm_version", as.character(packageVersion("ebnm")))
reg("metafor_version", as.character(packageVersion("metafor")))
reg("platform", R.version$platform)
reg("analysis_date", format(Sys.Date()))
reg("analysis_minutes", as.numeric(difftime(Sys.time(), started_all, units = "mins")), 1)

saveRDS(list(krw = fit_krw, twin = fit_twin, raudenbush = fits_raud), "data/fits.rds")
numbers <- do.call(rbind, lapply(registry, function(r) {
  value <- r$value
  shown <- if (is.numeric(value) && !is.na(r$digits)) formatC(value, format = "f", digits = r$digits)
           else as.character(value)
  data.frame(name = r$name, value = as.character(value), shown = shown, stringsAsFactors = FALSE)
}))
write.csv(numbers, "data/numbers.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "qa/sessionInfo-analysis.txt")
message("done: ", nrow(numbers), " numbers in data/numbers.csv")
