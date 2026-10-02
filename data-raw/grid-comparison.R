# Builds inst/examples/grid_comparison_example.rds, which
# vignette("choosing-a-grid") and vignette("diagnostics") read. Needs CmdStan
# and the loo package. Run from the package root:
#
#   Rscript data-raw/grid-comparison.R
#
# The file holds
#
#   raudenbush              compare_efron_grids() with its default arguments
#                           for the teacher expectancy studies (9 fits);
#   raudenbush_sensitivity  the same studies fitted for three widths of the
#                           grid and three numbers of spline functions, with
#                           the summaries that the vignette compares (9 fits);
#   two_groups              compare_efron_grids() with its default arguments
#                           for simulated data in which the true effects form
#                           two separated groups (9 fits);
#   two_groups_runs         three fits of those data: four spline functions
#                           with the default length of run and with a run four
#                           times as long, and six spline functions;
#   two_groups_sensitivity  those data fitted for three widths of the grid and
#                           two numbers of spline functions (6 fits);
#   two_groups_data         the simulated effects, estimates and standard
#                           errors.
#
# The simulated data follow the design of Lee and Sui (2025, Section 5.1) with
# reliability 0.7 and R = 9, here with 200 sites. The sampling variances run
# from (1 - I) / I / R to (1 - I) / I * R, so the largest standard error is
# R = 9 times the smallest.

devtools::load_all(".", quiet = TRUE)

# Summaries of one fit that compare_efron_grids() does not return.
fit_summaries <- function(fit) {
  diagnostics <- diagnose(fit)
  list(
    max_rhat = max(diagnostics$rhat, na.rm = TRUE),
    min_ess_bulk = min(diagnostics$ess_bulk, na.rm = TRUE),
    divergences = sum(diagnostics$divergences),
    treedepth_hits = sum(diagnostics$max_treedepth),
    draws = posterior::ndraws(fit$draws)
  )
}

# --- teacher expectancy studies ------------------------------------------------

raudenbush <- compare_efron_grids(
  raudenbush1985$yi, raudenbush1985$sei,
  seed = 1985L
)

largest <- which.max(raudenbush1985$yi)      # study 4
precise <- which.min(raudenbush1985$sei)     # study 18
stopifnot(largest == 4L, precise == 18L)

specifications <- expand.grid(expansion = c(0.25, 0.5, 1), M = c(6L, 8L, 10L))
sensitivity_fits <- lapply(
  seq_len(nrow(specifications)),
  function(row) {
    expansion <- specifications$expansion[row]
    M <- specifications$M[row]
    fit <- suppressWarnings(bayes_efron_fit(
      raudenbush1985$yi, raudenbush1985$sei,
      expansion = expansion, M = M, seed = 1985L
    ))
    grid <- fit$metadata$data_list$grid
    g <- colMeans(.bef_grid_draws(fit, "g"))
    summaries_g <- confint(fit, type = "g")
    intervals <- confint(fit)
    checks <- fit_summaries(fit)
    fit_loo <- loo::loo(fit)
    list(
      loo = fit_loo,
      row = data.frame(
        expansion = expansion,
        M = M,
        grid_lower = min(grid),
        grid_upper = max(grid),
        mean_g = summaries_g$point[summaries_g$site == "mean_g"],
        sd_g = summaries_g$point[summaries_g$site == "sd_g"],
        above_0.5 = sum(g[grid > 0.5]),
        study4_mean = unname(coef(fit)[largest]),
        study4_lower = intervals$lower[largest],
        study4_upper = intervals$upper[largest],
        study18_mean = unname(coef(fit)[precise]),
        elpd_loo = fit_loo$estimates["elpd_loo", "Estimate"],
        max_rhat = checks$max_rhat,
        divergences = checks$divergences,
        check.names = FALSE
      )
    )
  }
)
raudenbush_sensitivity <- do.call(rbind, lapply(sensitivity_fits, `[[`, "row"))
# Differences in elpd_loo from the best of the nine fits, with standard errors
# from the pointwise values, as in compare_efron_grids().
differences <- .bef_elpd_differences(lapply(sensitivity_fits, `[[`, "loo"))
raudenbush_sensitivity$elpd_diff <- differences$elpd_diff
raudenbush_sensitivity$se_diff <- differences$se_diff

# --- simulated data with two groups of effects -----------------------------------

set.seed(2025)
K <- 200L
n_left <- round(K / 3)
theta <- c(runif(n_left, -1.7, -0.7), runif(K - n_left, 0.7, 2.7))
reliability <- 0.7
R <- 9
variance_range <- (1 - reliability) / reliability * c(1 / R, R)
sigma <- sqrt(exp(seq(log(variance_range[1]), log(variance_range[2]), length.out = K)))
sigma <- sigma[sample.int(K)]
theta_hat <- rnorm(K, mean = theta, sd = sigma)

two_groups <- compare_efron_grids(theta_hat, sigma, seed = 2025L)

one_run <- function(M, expansion = 0.5, iter_warmup = 1000L, iter_sampling = 3000L) {
  fit <- suppressWarnings(bayes_efron_fit(
    theta_hat, sigma, M = M, expansion = expansion, seed = 2025L,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling
  ))
  intervals <- confint(fit)
  summaries_g <- confint(fit, type = "g")
  checks <- fit_summaries(fit)
  data.frame(
    M = M,
    expansion = expansion,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    draws = checks$draws,
    max_rhat = checks$max_rhat,
    min_ess_bulk = checks$min_ess_bulk,
    treedepth_hits = checks$treedepth_hits,
    mean_g = summaries_g$point[summaries_g$site == "mean_g"],
    sd_g = summaries_g$point[summaries_g$site == "sd_g"],
    coverage = mean(intervals$lower <= theta & theta <= intervals$upper),
    rmse = sqrt(mean((coef(fit) - theta)^2)),
    mean_width = mean(intervals$upper - intervals$lower),
    sampler_seconds = fit$metadata$runtime_seconds
  )
}

# Four spline functions with the default length of run and with a run four
# times as long, and six spline functions.
two_groups_runs <- rbind(
  one_run(M = 4L),
  one_run(M = 4L, iter_warmup = 4000L, iter_sampling = 12000L),
  one_run(M = 6L)
)

# Three widths of the grid and two numbers of spline functions.
two_groups_specifications <- expand.grid(expansion = c(0.05, 0.25, 0.5), M = c(6L, 8L))
two_groups_sensitivity <- do.call(rbind, lapply(
  seq_len(nrow(two_groups_specifications)),
  function(row) {
    one_run(
      M = two_groups_specifications$M[row],
      expansion = two_groups_specifications$expansion[row]
    )
  }
))

print(raudenbush)
print(raudenbush_sensitivity)
print(two_groups)
print(two_groups_runs)
print(two_groups_sensitivity)

out <- list(
  raudenbush = raudenbush,
  raudenbush_sensitivity = raudenbush_sensitivity,
  two_groups = two_groups,
  two_groups_runs = two_groups_runs,
  two_groups_sensitivity = two_groups_sensitivity,
  two_groups_data = data.frame(theta = theta, theta_hat = theta_hat, sigma = sigma),
  cmdstan_version = as.character(cmdstanr::cmdstan_version()),
  created = format(Sys.Date())
)
saveRDS(out, "inst/examples/grid_comparison_example.rds", version = 3)
message("inst/examples/grid_comparison_example.rds: ",
        file.size("inst/examples/grid_comparison_example.rds"), " bytes")
