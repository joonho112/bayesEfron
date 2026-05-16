# Maya smoke path: escalc -> as_bef_data -> bayes_efron_fit -> confint -> plot.

## ---- maya-packages
if (!requireNamespace("metafor", quietly = TRUE)) {
  stop("The Maya example requires the suggested package `metafor`.", call. = FALSE)
}
if (!requireNamespace("bayesEfron", quietly = TRUE)) {
  stop("The Maya example requires the package `bayesEfron`.", call. = FALSE)
}

bayes_efron_fit <- getExportedValue("bayesEfron", "bayes_efron_fit")
as_bef_data <- getExportedValue("bayesEfron", "as_bef_data")
diagnose <- getExportedValue("bayesEfron", "diagnose")

## ---- maya-data
my_multisite_data <- data.frame(
  site_name = paste0("site_", 1:10),
  m_treatment = c(0.42, 0.18, 0.05, 0.31, 0.54, 0.12, 0.39, 0.27, 0.48, 0.16),
  sd_treatment = c(1.05, 1.12, 0.98, 1.08, 1.15, 1.03, 1.10, 0.96, 1.20, 1.07),
  n_treatment = c(42L, 39L, 44L, 41L, 40L, 43L, 37L, 45L, 38L, 46L),
  m_control = c(0.05, 0.02, -0.08, 0.04, 0.10, -0.03, 0.11, 0.00, 0.08, -0.02),
  sd_control = c(1.01, 1.09, 1.04, 1.02, 1.12, 1.00, 1.06, 0.99, 1.16, 1.05),
  n_control = c(40L, 38L, 42L, 39L, 41L, 44L, 36L, 43L, 39L, 45L)
)

## ---- maya-escalc
dat <- metafor::escalc(
  measure = "MD",
  m1i = m_treatment,
  sd1i = sd_treatment,
  n1i = n_treatment,
  m2i = m_control,
  sd2i = sd_control,
  n2i = n_control,
  data = my_multisite_data,
  slab = site_name
)

bef_dat <- as_bef_data(dat)

## ---- maya-fit
# Smoke-scale controls keep this example fast. Increase chains and iterations
# for a real analysis.
fit <- bayes_efron_fit(
  theta_hat = bef_dat$theta_hat,
  sigma = bef_dat$sigma,
  L = 51L,
  M = 3L,
  chains = 1L,
  iter_warmup = 150L,
  iter_sampling = 4L,
  seed = 1234L
)

## ---- maya-outputs
fit_summary <- summary(fit)
theta_ci <- confint(fit, level = 0.9, type = "theta")
fit_diagnostic <- diagnose(fit)
fit_plot <- plot(fit, type = "caterpillar")

maya_result <- list(
  data = bef_dat,
  fit = fit,
  summary = fit_summary,
  theta_ci = theta_ci,
  diagnostic = fit_diagnostic,
  plot = fit_plot
)

invisible(maya_result)
