# Rebuilds the figure of the simulated sites in the README from the stored
# fit with the plot function of the package. No model is fitted. Needs the
# ragg package. The other figure of the README is made when README.Rmd is
# knitted (`devtools::build_readme()`). Run from the package root:
#
#   Rscript data-raw/readme-figures.R
suppressMessages(devtools::load_all(".", quiet = TRUE))

save_figure <- function(plot, name) {
  ggplot2::ggsave(file.path("man", "figures", name), plot, device = ragg::agg_png,
                  width = 7, height = 4.2, units = "in", dpi = 144)
  message(sprintf("man/figures/%s: %d bytes", name,
                  file.size(file.path("man", "figures", name))))
}

# The 1,500 simulated sites: the fit, the normal model and the true effects.
save_figure(
  plot_g(twin_towers_fit, truth = twin_towers$theta_true, normal = TRUE),
  "twin-towers.png"
)
