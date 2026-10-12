# ==============================================================================
# figures.R: Figures 2 to 5 of the SoftwareX article on bayesEfron 0.4.0, and
# the figures of the supplement. Reads only what scripts/analysis.R wrote in
# data/, and asserts every headline number that a figure displays against
# data/numbers.csv before a file is written.
#
# Run from the paper folder root:   Rscript scripts/figures.R
#
# House style (Applied Geography / Geometry papers): Helvetica, theme_minimal
# with a grey90 grid, bold (a)(b)(c) panel tags, the Okabe-Ito palette with
# fixed semantic assignments, portable PDF + ragg PNG. The figures are drawn at
# the text width of the elsarticle preprint (137 mm), so that they are
# included at their natural size.
# ==============================================================================

root <- normalizePath(".")
stopifnot(file.exists(file.path(root, "scripts", "figures.R")))
.libPaths(c(file.path(root, "r-library"), .libPaths()))
suppressPackageStartupMessages({
  library(bayesEfron)
  library(ggplot2)
  library(patchwork)
  library(scales)
})
stopifnot(packageVersion("bayesEfron") == "0.4.0")
print(packageVersion("bayesEfron"))

FIGDIR <- "figures"; PNGDIR <- "figures/png"; OSMDIR <- "osm/figures"
for (d in c(FIGDIR, PNGDIR, OSMDIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

numbers <- read.csv("data/numbers.csv", stringsAsFactors = FALSE)
num <- function(name) {
  v <- numbers$value[numbers$name == name]
  stopifnot(length(v) == 1L)
  suppressWarnings(as.numeric(v))
}
assert <- function(cond, msg) if (!isTRUE(cond)) stop(msg, call. = FALSE)
near <- function(x, y, tol = 1e-8) is.finite(x) && is.finite(y) && abs(x - y) <= tol

krw <- readRDS("data/krw.rds")
twin <- readRDS("data/twin.rds")
evidence <- readRDS("data/evidence.rds")
gridcmp <- readRDS("data/grid-comparison.rds")

# ---- style -------------------------------------------------------------------
BASE_SIZE <- 8
BASE_FAMILY <- "Helvetica"
PAL <- list(
  ink = "#1A1A1A", ink2 = "#4D4D4D", ink3 = "#767676", grid = "grey90",
  blue = "#0072B2", sky = "#56B4E9", vermillion = "#D55E00", green = "#009E73",
  orange = "#E69F00", purple = "#CC79A7", bar = "#C9CDD2", band = "#0072B2"
)
# Fixed semantic colors across all figures.
METHOD_COL <- c("bayesEfron" = PAL$blue, "bayesEfron g, plugged in" = PAL$sky,
                "Normal (plug-in)" = PAL$vermillion, "NPMLE" = PAL$green,
                "Unimodal" = PAL$purple, "Efron's estimator" = PAL$orange, "Known g" = PAL$ink,
                "Normal" = PAL$vermillion)
theme_paper <- function(base_size = BASE_SIZE) {
  theme_minimal(base_size = base_size, base_family = BASE_FAMILY) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = PAL$grid, linewidth = 0.28),
      axis.ticks = element_blank(),
      axis.title = element_text(size = base_size, color = PAL$ink),
      axis.text = element_text(size = base_size - 0.6, color = PAL$ink2),
      strip.text = element_text(size = base_size + 0.2, face = "bold", color = PAL$ink,
                                hjust = 0, margin = margin(b = 3)),
      legend.title = element_text(size = base_size - 0.4, color = PAL$ink),
      legend.text = element_text(size = base_size - 0.8, color = PAL$ink2),
      legend.key.height = unit(3, "mm"), legend.key.width = unit(3.2, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.title = element_blank(), plot.subtitle = element_blank(),
      plot.tag = element_text(size = base_size + 1.2, face = "bold", family = BASE_FAMILY),
      plot.margin = margin(2, 4, 2, 2)
    )
}
# A build can advertise Cairo support while its runtime shared library is
# missing. Require a loadable Cairo version, then fall back to base PDF.
cairo_version <- if (capabilities("cairo")) {
  tryCatch(suppressWarnings(grDevices::cairoVersion()), error = function(e) "")
} else ""
use_cairo <- length(cairo_version) == 1L && nzchar(cairo_version)
dev_paper_pdf <- function(filename, width, height, ...) {
  if (use_cairo) {
    grDevices::cairo_pdf(filename = filename, width = width, height = height,
                        family = BASE_FAMILY, ...)
  } else {
    grDevices::pdf(file = filename, width = width, height = height,
                  family = BASE_FAMILY, useDingbats = FALSE, ...)
  }
}
message("PDF device: ", if (use_cairo) "cairo_pdf" else "pdf")
# Text in the PDF figures is drawn as vector outlines with showtext, so that
# the PDFs carry no font to embed (neither the quartz nor the base pdf device
# embeds Helvetica, and the journal asks for embedded fonts). The outlines come
# from the two Helvetica faces of macOS when they have been extracted to
# qa/fonts/ (see the README), else from Arial, which has the same metrics.
# Without showtext the PDFs are written as before. The PNG proofs keep ragg's
# own text rendering.
OUTLINE_FONTS <- list(
  c(regular = file.path("qa", "fonts", "Helvetica-Regular.ttf"),
    bold = file.path("qa", "fonts", "Helvetica-Bold.ttf")),
  c(regular = "/System/Library/Fonts/Supplemental/Arial.ttf",
    bold = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"))
use_showtext <- FALSE
if (requireNamespace("showtext", quietly = TRUE) && requireNamespace("sysfonts", quietly = TRUE)) {
  for (faces in OUTLINE_FONTS) if (all(file.exists(faces))) {
    sysfonts::font_add(BASE_FAMILY, regular = faces[["regular"]], bold = faces[["bold"]])
    use_showtext <- TRUE
    message("PDF text drawn as outlines from ", basename(faces[["regular"]]))
    break
  }
}
if (!use_showtext) message("showtext or its font files not found: PDF text keeps the device font")
save_fig <- function(plot, name, width_mm, height_mm, dir = FIGDIR) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  if (use_showtext) { showtext::showtext_auto(TRUE); showtext::showtext_opts(dpi = 72) }
  ggsave(file.path(dir, paste0(name, ".pdf")), plot, width = w, height = h,
         device = dev_paper_pdf)
  if (use_showtext) showtext::showtext_auto(FALSE)
  ggsave(file.path(PNGDIR, paste0(name, ".png")), plot, width = w, height = h,
         dpi = 400, device = ragg::agg_png, bg = "white")
  message("wrote ", name, " (", width_mm, " x ", height_mm, " mm)")
}
tags <- function(p) p + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
MM_FULL <- 137
label_text <- function(..., size = 2.4, color = PAL$ink2)
  annotate("text", ..., size = size, color = color, family = BASE_FAMILY, lineheight = 0.95)

# ==============================================================================
# Figure 2: the callback gaps of 97 firms
# ==============================================================================
pp <- 100
d <- krw$density
dens <- data.frame(x = d$value * pp, mean = d$mean / pp, lower = d$lower / pp, upper = d$upper / pp)
normal_curve <- data.frame(x = dens$x,
                           y = dnorm(dens$x, krw$normal$mean * pp, krw$normal$sd * pp))
est <- data.frame(x = krw$data$yi * pp)
assert(near(round(num("krw_mean_g_pp"), 2), round(100 * krw$g_ci["mean_g", "mean"], 2)),
       "fig2: mean of g drifted")
assert(near(round(num("krw_reml_sd_pp"), 2), round(100 * krw$normal$sd, 2)), "fig2: REML sd drifted")
top_y <- max(dens$upper)
p2a <- ggplot() +
  geom_histogram(data = est, aes(x = x, y = after_stat(density)), binwidth = 1,
                 boundary = 0, fill = PAL$bar, color = "white", linewidth = 0.2) +
  geom_ribbon(data = dens, aes(x = x, ymin = lower, ymax = upper), fill = PAL$band, alpha = 0.18) +
  geom_line(data = normal_curve, aes(x = x, y = y), color = PAL$vermillion, linewidth = 0.55,
            linetype = "22") +
  geom_line(data = dens, aes(x = x, y = mean), color = PAL$blue, linewidth = 0.75) +
  geom_vline(xintercept = 0, color = PAL$ink3, linewidth = 0.3) +
  label_text(x = 4.6, y = 0.97 * top_y, hjust = 0, vjust = 1, color = PAL$blue, size = 2.1,
             label = sprintf("g, with 90%% band\nmean %.2f, SD %.2f",
                             num("krw_mean_g_pp"), num("krw_sd_g_pp"))) +
  label_text(x = 4.6, y = 0.74 * top_y, hjust = 0, vjust = 1, color = PAL$vermillion, size = 2.1,
             label = sprintf("Normal (REML)\nmean %.2f, SD %.2f",
                             num("krw_reml_mean_pp"), num("krw_reml_sd_pp"))) +
  label_text(x = 4.6, y = 0.51 * top_y, hjust = 0, vjust = 1, color = PAL$ink2, size = 2.1,
             label = sprintf("bars: 97 estimates\nSD %.2f", num("krw_sd_estimates_pp"))) +
  scale_x_continuous(breaks = seq(-8, 16, 4), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  coord_cartesian(xlim = c(-8.3, 15.9)) +
  labs(x = "Callback gap, white minus Black (percentage points)", y = "Density") +
  theme_paper()

# (b) caterpillar of the firms with the eight largest and eight smallest posterior means
sites50 <- krw$sites50
ord <- order(sites50$mean)
pick <- c(head(ord, 8), tail(ord, 8))
cat_df <- sites50[pick, ]
cat_df$group <- factor(rep(c("Smallest posterior means", "Largest posterior means"), each = 8),
                       levels = c("Largest posterior means", "Smallest posterior means"))
cat_df$label <- factor(cat_df$label, levels = cat_df$label[order(cat_df$mean)])
assert(identical(as.character(cat_df$label[which.max(cat_df$mean)]),
                 numbers$value[numbers$name == "krw_top_firm"]), "fig2b: top firm drifted")
p2b <- ggplot(cat_df, aes(y = label)) +
  geom_vline(xintercept = 0, color = PAL$ink3, linewidth = 0.3) +
  geom_segment(aes(x = lower, xend = upper, yend = label), color = PAL$blue, linewidth = 0.45,
               alpha = 0.6) +
  geom_segment(aes(x = inner_lower, xend = inner_upper, yend = label), color = PAL$blue,
               linewidth = 1.4) +
  geom_point(aes(x = estimate), shape = 21, fill = "white", color = PAL$ink2, size = 1.6,
             stroke = 0.5) +
  geom_point(aes(x = mean), color = PAL$blue, size = 1.8) +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = seq(-4, 12, 4)) +
  labs(x = "Callback gap (percentage points)", y = NULL) +
  theme_paper() +
  theme(axis.text.y = element_text(size = BASE_SIZE - 1.2),
        strip.text = element_text(size = BASE_SIZE - 0.4),
        panel.spacing.y = unit(2, "mm"))

# (c) shrinkage: posterior mean against estimate, colored by the standard error
sh <- krw$sites
p2c <- ggplot(sh, aes(x = estimate, y = mean)) +
  geom_abline(slope = 1, intercept = 0, color = PAL$ink3, linewidth = 0.35, linetype = "22") +
  geom_hline(yintercept = num("krw_mean_g_pp"), color = PAL$ink3, linewidth = 0.25) +
  geom_point(aes(fill = se), shape = 21, color = "white", size = 2.0, stroke = 0.25) +
  scale_fill_gradient(low = "#084571", high = "#CBDEEC", name = "SE (pp)") +
  coord_equal(xlim = c(-5.5, 10.5), ylim = c(-5.5, 10.5)) +
  scale_x_continuous(breaks = seq(-4, 8, 4)) + scale_y_continuous(breaks = seq(-4, 8, 4)) +
  label_text(x = 9.8, y = 9.0, label = "no shrinkage", hjust = 1, size = 2.3, color = PAL$ink3) +
  label_text(x = -5.0, y = num("krw_mean_g_pp") + 0.8, label = "mean of g", hjust = 0,
             size = 2.3, color = PAL$ink3) +
  labs(x = "Estimate (percentage points)", y = "Posterior mean (percentage points)") +
  theme_paper() +
  theme(legend.position = c(0.84, 0.24), legend.background = element_rect(fill = "white", color = NA),
        legend.key.height = unit(2.4, "mm"), legend.key.width = unit(2.4, "mm"))

# (d) the probability of a positive gap, every firm, sorted
pos <- data.frame(firm = names(krw$p_pos), p = as.numeric(krw$p_pos))
pos <- pos[order(-pos$p), ]
pos$rank <- seq_len(nrow(pos))
pos$above <- pos$p >= 0.95
n_above <- sum(pos$above)
assert(n_above == num("krw_n_prob_ge_95"), "fig2d: count above 0.95 drifted")
p2d <- ggplot(pos, aes(x = rank, y = p)) +
  geom_hline(yintercept = 0.95, color = PAL$orange, linewidth = 0.45, linetype = "22") +
  geom_point(aes(color = above), size = 1.4) +
  scale_color_manual(values = c(`TRUE` = PAL$blue, `FALSE` = PAL$ink3), guide = "none") +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 97)) +
  label_text(x = n_above + 3, y = 0.58, hjust = 0, vjust = 1, color = PAL$blue, size = 2.3,
             label = sprintf("%d of 97 firms:\nP(gap > 0) of at least 0.95", n_above)) +
  label_text(x = 96, y = 0.985, hjust = 1, vjust = 1, label = "0.95", size = 2.3,
             color = PAL$orange) +
  labs(x = "Firms, sorted by the probability", y = "Posterior probability of a positive gap") +
  theme_paper()

fig2 <- tags((p2a | p2b) / (p2c | p2d) + plot_layout(widths = c(1.12, 1)))
save_fig(fig2, "fig2_firms", MM_FULL, 124)

# ==============================================================================
# Figure 3: the Twin Towers simulation, 1,500 sites with known effects
# ==============================================================================
sim <- twin$data
d <- twin$density
dens3 <- data.frame(x = d$value, mean = d$mean, lower = d$lower, upper = d$upper)
normal3 <- data.frame(x = seq(-3, 3.7, length.out = 300))
normal3$y <- dnorm(normal3$x, twin$normal$mean, twin$normal$sd)
truth_df <- data.frame(x = sim$theta_true)
est_density <- as.data.frame(density(sim$yi, bw = "SJ", from = -3, to = 3.7)[c("x", "y")])
assert(near(round(num("twin_prob_between"), 2), 0.06), "fig3: probability between groups drifted")
p3a <- ggplot() +
  annotate("rect", xmin = -0.7, xmax = 0.7, ymin = 0, ymax = Inf, fill = PAL$orange, alpha = 0.10) +
  geom_histogram(data = truth_df, aes(x = x, y = after_stat(density)), binwidth = 0.2,
                 boundary = -1.7, fill = PAL$bar, color = "white", linewidth = 0.2) +
  geom_line(data = est_density, aes(x = x, y = y), color = PAL$ink3, linewidth = 0.45,
            linetype = "dotted") +
  geom_ribbon(data = dens3, aes(x = x, ymin = lower, ymax = upper), fill = PAL$band, alpha = 0.18) +
  geom_line(data = normal3, aes(x = x, y = y), color = PAL$vermillion, linewidth = 0.55,
            linetype = "22") +
  geom_line(data = dens3, aes(x = x, y = mean), color = PAL$blue, linewidth = 0.75) +
  label_text(x = 0.95, y = 0.655, vjust = 1, size = 2.2, color = PAL$ink2,
             label = sprintf("no true effect here:\ng puts %.2f, normal %.2f",
                             num("twin_prob_between"), num("twin_normal_between"))) +
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.68), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-3.05, 3.75)) +
  labs(x = "Effect", y = "Density",
       caption = "Blue: g, 90% band; dashed: normal (REML)\nBars: true effects; dotted: estimates") +
  theme_paper() +
  theme(plot.caption = element_text(size = BASE_SIZE - 0.8, color = PAL$ink2,
                                     hjust = 0, margin = margin(t = 3)))

# (b) shrinkage under the two models, colored by the group of the true effect
sb <- rbind(
  data.frame(model = "bayesEfron", estimate = sim$yi, mean = twin$sites$mean,
             group = ifelse(sim$theta_true < 0, "lower group", "upper group")),
  data.frame(model = "Normal (REML)", estimate = sim$yi, mean = twin$normal_means,
             group = ifelse(sim$theta_true < 0, "lower group", "upper group"))
)
sb$model <- factor(sb$model, levels = c("bayesEfron", "Normal (REML)"))
p3b <- ggplot(sb, aes(x = estimate, y = mean, color = group)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -0.7, ymax = 0.7, fill = PAL$orange, alpha = 0.10) +
  geom_abline(slope = 1, intercept = 0, color = PAL$ink3, linewidth = 0.3, linetype = "22") +
  geom_point(size = 0.5, alpha = 0.55) +
  scale_color_manual(values = c("lower group" = PAL$vermillion, "upper group" = PAL$blue),
                     name = NULL) +
  facet_wrap(~ model, nrow = 1) +
  coord_cartesian(xlim = c(-6, 7), ylim = c(-2.5, 3.5)) +
  labs(x = "Estimate", y = "Posterior mean") +
  theme_paper() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.key.height = unit(2.5, "mm"), legend.box.margin = margin(-6, 0, 0, 0)) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

# (c) the posterior distribution of an imprecise site between the groups and of the most precise site
pc <- twin$posterior_curves
shown <- c(num("twin_site_a"), 304)
pc <- pc[pc$site %in% shown, ]
site_tab <- twin$three
site_tab <- site_tab[match(shown, site_tab$site), ]
lik <- do.call(rbind, lapply(seq_len(nrow(site_tab)), function(i) {
  x <- seq(-3, 3.7, length.out = 400)
  data.frame(site = site_tab$site[i], value = x,
             density = dnorm(x, site_tab$estimate[i], site_tab$se[i]))
}))
gcurve <- data.frame(site = rep(shown, each = nrow(dens3)), value = rep(dens3$x, 2),
                     density = rep(dens3$mean, 2))
panel_lab <- setNames(sprintf("Site %d (SE %.2f)", site_tab$site, site_tab$se), site_tab$site)
to_panel <- function(site) factor(panel_lab[as.character(site)], levels = panel_lab)
pc$panel <- to_panel(pc$site); lik$panel <- to_panel(lik$site); gcurve$panel <- to_panel(gcurve$site)
means <- data.frame(panel = to_panel(site_tab$site), mean = site_tab$mean,
                    lower = site_tab$lower, upper = site_tab$upper)
p3c <- ggplot() +
  geom_area(data = pc, aes(x = value, y = posterior), fill = PAL$blue, alpha = 0.25) +
  geom_line(data = pc, aes(x = value, y = posterior, linetype = "Posterior"),
            color = PAL$blue, linewidth = 0.6) +
  geom_line(data = lik, aes(x = value, y = density, linetype = "Likelihood"),
            color = PAL$ink2, linewidth = 0.45) +
  geom_line(data = gcurve, aes(x = value, y = density, linetype = "g"),
            color = PAL$ink3, linewidth = 0.4) +
  geom_segment(data = means, aes(x = lower, xend = upper, y = -0.08, yend = -0.08),
               color = PAL$blue, linewidth = 0.6) +
  geom_point(data = means, aes(x = mean, y = -0.08), color = PAL$blue, size = 1.4) +
  facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  scale_linetype_manual(values = c("Posterior" = "solid", "Likelihood" = "22", "g" = "dotted"),
                        breaks = c("Posterior", "Likelihood", "g"), name = NULL) +
  guides(linetype = guide_legend(nrow = 1,
    override.aes = list(color = c(PAL$blue, PAL$ink2, PAL$ink3)))) +
  scale_x_continuous(limits = c(-3.05, 3.75), breaks = seq(-3, 3, 1)) +
  labs(x = "Effect", y = "Density") +
  theme_paper() + theme(strip.text = element_text(size = BASE_SIZE - 0.4),
                         legend.position = "bottom", legend.direction = "horizontal",
                         legend.key.width = unit(4.5, "mm"),
                         legend.key.height = unit(2.5, "mm"),
                         legend.box.margin = margin(-5, 0, 0, 0))

# (d) calibration: the share of the true effects inside the interval at each level
cal <- twin$calibration
at90 <- abs(cal$level - 0.9) < 1e-9
assert(near(round(100 * cal$inside[at90], 1), round(num("twin_inside_90_pct"), 1)),
       "fig3d: 90% inclusion drifted")
p3d <- ggplot(cal, aes(x = level, y = inside)) +
  geom_abline(slope = 1, intercept = 0, color = PAL$ink3, linewidth = 0.3, linetype = "22") +
  geom_line(color = PAL$blue, linewidth = 0.55) +
  geom_point(color = PAL$blue, size = 1.5) +
  geom_point(data = cal[at90, ], color = PAL$orange, size = 2.6, shape = 21, stroke = 0.8,
             fill = NA) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent_format(accuracy = 1),
                     expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent_format(accuracy = 1),
                     expand = c(0.01, 0)) +
  coord_equal() +
  labs(x = "Level of the posterior interval", y = "True effects inside the interval",
       caption = sprintf("90%% intervals contain %.1f%%\nof the 1,500 true effects",
                         num("twin_inside_90_pct"))) +
  theme_paper() +
  theme(plot.caption = element_text(size = BASE_SIZE - 1.5, color = PAL$ink2,
                                     hjust = 1, margin = margin(t = 3)))

fig3 <- tags((p3a | p3b) / (p3c | p3d) + plot_layout(widths = c(1.12, 1)))
save_fig(fig3, "fig3_twin_towers", MM_FULL, 124)

# ==============================================================================
# Figure 4: coverage, length and error over simulated data sets, by design
# ==============================================================================
pm <- evidence$panels_methods
few <- evidence$few
eq <- evidence$equal
design_levels <- c(
  "Two groups, K = 100\n(100 data sets)", "Two groups, K = 1,500\n(20 data sets)",
  "19 studies, normal SD 0.14\n(200 data sets)", "19 studies, normal SD 0.31\n(200 data sets)",
  "19 studies, fitted g\n(200 data sets)", "Two groups, equal SE, K = 1,500\n(20 data sets)"
)
rows <- list()
add <- function(design, method, inclusion, inclusion_se, length, length_se, rmse, rmse_se) {
  rows[[length(rows) + 1L]] <<- data.frame(design = design, method = method, inclusion = inclusion,
                                           inclusion_se = inclusion_se, length = length,
                                           length_se = length_se, rmse = rmse, rmse_se = rmse_se)
}
rename_pm <- c("bayesEfron" = "bayesEfron", "g of bayesEfron, plugged in" = "bayesEfron g, plugged in",
               "Normal" = "Normal (plug-in)", "NPMLE" = "NPMLE")
for (K in c(100, 1500)) {
  design <- design_levels[if (K == 100) 1 else 2]
  for (m in names(rename_pm)) {
    r <- pm[pm$K == K & pm$method == m, ]
    add(design, rename_pm[[m]], r$inclusion, r$inclusion_se, r$mean_length, r$mean_length_se,
        r$rmse, r$rmse_se)
  }
}
few_methods <- c("bayesEfron" = "bayesEfron", "normal" = "Normal (plug-in)", "known g" = "Known g")
few_scen <- c("normal, SD 0.14" = design_levels[3], "normal, SD 0.31" = design_levels[4],
              "fitted g" = design_levels[5])
for (s in names(few_scen)) for (m in names(few_methods)) {
  r <- few[few$scenario == s & few$method == m, ]
  add(few_scen[[s]], few_methods[[m]], r$inclusion, r$inclusion_se, r$length, r$length_se,
      r$rmse, r$rmse_se)
}
eq_methods <- c("bayesEfron" = "bayesEfron", "g of bayesEfron, plugged in" = "bayesEfron g, plugged in",
                "Efron's estimator" = "Efron's estimator")
for (m in names(eq_methods)) {
  r <- eq[eq$M == 6 & eq$method == m, ]
  add(design_levels[6], eq_methods[[m]], r$inclusion, r$inclusion_se, r$mean_length,
      r$mean_length_se, r$rmse, r$rmse_se)
}
ev <- do.call(rbind, rows)
ev$design <- factor(ev$design, levels = rev(design_levels))
method_levels <- c("bayesEfron", "bayesEfron g, plugged in", "Normal (plug-in)", "NPMLE",
                   "Efron's estimator", "Known g")
ev$method <- factor(ev$method, levels = method_levels)
ref <- ev[ev$method == "bayesEfron", c("design", "length", "rmse")]
names(ref)[2:3] <- c("length_ref", "rmse_ref")
ev <- merge(ev, ref, by = "design")
ev$length_ratio <- ev$length / ev$length_ref
ev$rmse_ratio <- ev$rmse / ev$rmse_ref
# Paired delta-method errors use the same saved data sets for each method and
# its bayesEfron reference; scripts/ratio-mcse.R generates this portable CSV.
ratio_mcse <- read.csv("data/ratio-mcse.csv", stringsAsFactors = FALSE)
ratio_key <- function(x) paste(as.character(x$design), as.character(x$method), sep = " | ")
assert(!anyDuplicated(ratio_key(ratio_mcse)), "fig4: duplicate ratio MCSE keys")
j <- match(ratio_key(ev), ratio_key(ratio_mcse))
assert(!anyNA(j) && nrow(ev) == nrow(ratio_mcse), "fig4: ratio MCSE keys do not match")
assert(max(abs(ev$length_ratio - ratio_mcse$length_ratio[j])) < 1e-12 &&
         max(abs(ev$rmse_ratio - ratio_mcse$rmse_ratio[j])) < 1e-12,
       "fig4: paired-ratio point estimates drifted")
ev$length_ratio_se <- ratio_mcse$length_ratio_se[j]
ev$rmse_ratio_se <- ratio_mcse$rmse_ratio_se[j]
ratio_limits <- function(ratio, se) {
  bounds <- range(c(ratio - 2 * se, ratio + 2 * se, 1))
  assert(all(is.finite(bounds)), "fig4: non-finite ratio bounds")
  bounds + c(-1, 1) * max(0.015, 0.04 * diff(bounds))
}
length_limits <- ratio_limits(ev$length_ratio, ev$length_ratio_se)
rmse_limits <- ratio_limits(ev$rmse_ratio, ev$rmse_ratio_se)
assert(near(round(ev$inclusion[ev$method == "bayesEfron" & ev$design == design_levels[2]], 3),
            round(num("cov_K1500"), 3)), "fig4: coverage K = 1500 drifted")
assert(near(round(ev$inclusion[ev$method == "bayesEfron" & ev$design == design_levels[6]], 3),
            round(num("cov_equal_M6"), 3)), "fig4: coverage equal SE drifted")
METHOD_SHP <- c("bayesEfron" = 16, "bayesEfron g, plugged in" = 15, "Normal (plug-in)" = 17,
                "NPMLE" = 18, "Efron's estimator" = 15, "Known g" = 1)
dodge <- position_dodge(width = 0.7)
base4 <- function(p) p + theme_paper() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.box.margin = margin(-4, 0, 0, 0), panel.grid.major.y = element_blank(),
        axis.title.x = element_text(size = BASE_SIZE - 0.4))
p4a <- ggplot(ev, aes(x = inclusion, y = design, color = method, shape = method)) +
  geom_vline(xintercept = 0.90, color = PAL$ink2, linewidth = 0.4) +
  geom_errorbar(aes(xmin = inclusion - 2 * inclusion_se, xmax = inclusion + 2 * inclusion_se),
                width = 0, position = dodge, linewidth = 0.4, orientation = "y") +
  geom_point(position = dodge, size = 1.8) +
  scale_color_manual(values = METHOD_COL[method_levels], breaks = method_levels, name = NULL) +
  scale_shape_manual(values = METHOD_SHP, breaks = method_levels, name = NULL) +
  scale_x_continuous(breaks = seq(0.6, 0.95, 0.1), labels = percent_format(accuracy = 1)) +
  labs(x = "Coverage of the 90% intervals", y = NULL)
p4a <- base4(p4a) + theme(axis.text.y = element_text(size = BASE_SIZE - 1, lineheight = 0.9))
p4b <- ggplot(ev, aes(x = length_ratio, y = design, color = method, shape = method)) +
  geom_vline(xintercept = 1, color = PAL$ink2, linewidth = 0.4) +
  geom_errorbar(aes(xmin = length_ratio - 2 * length_ratio_se, xmax = length_ratio + 2 * length_ratio_se),
                width = 0, position = dodge, linewidth = 0.4, orientation = "y") +
  geom_point(position = dodge, size = 1.8) +
  scale_color_manual(values = METHOD_COL[method_levels], breaks = method_levels, guide = "none") +
  scale_shape_manual(values = METHOD_SHP, breaks = method_levels, guide = "none") +
  scale_x_continuous(breaks = c(0.8, 1.0, 1.2), limits = length_limits) +
  labs(x = "Length ratio", y = NULL)
p4b <- base4(p4b) + theme(axis.text.y = element_blank())
p4c <- ggplot(ev, aes(x = rmse_ratio, y = design, color = method, shape = method)) +
  geom_vline(xintercept = 1, color = PAL$ink2, linewidth = 0.4) +
  geom_errorbar(aes(xmin = rmse_ratio - 2 * rmse_ratio_se, xmax = rmse_ratio + 2 * rmse_ratio_se),
                width = 0, position = dodge, linewidth = 0.4, orientation = "y") +
  geom_point(position = dodge, size = 1.8) +
  scale_color_manual(values = METHOD_COL[method_levels], breaks = method_levels, guide = "none") +
  scale_shape_manual(values = METHOD_SHP, breaks = method_levels, guide = "none") +
  scale_x_continuous(breaks = c(0.9, 1.0, 1.1), limits = rmse_limits) +
  labs(x = "RMSE ratio", y = NULL)
p4c <- base4(p4c) + theme(axis.text.y = element_blank())
fig4 <- tags((p4a | p4b | p4c) + plot_layout(widths = c(1.6, 1, 1), guides = "collect") &
               theme(legend.position = "bottom", plot.margin = margin(2, 6, 2, 2)))
save_fig(fig4, "fig4_evidence", MM_FULL, 72)
write.csv(ev, "data/figure4-table.csv", row.names = FALSE)

# ==============================================================================
# Figure 5: the range of the grid with 19 studies, and what the uncertainty
# about g adds to a site
# ==============================================================================
tab <- gridcmp$table
gl <- gridcmp$g
range_lab <- c("1" = "0.25 (narrow)", "2" = "0.5 (default)", "3" = "1 (wide)")
specs <- as.character(tab$spec[order(tab$expansion)])
curves <- do.call(rbind, lapply(specs, function(id) cbind(gl[[id]], range = range_lab[[id]])))
curves$range <- factor(curves$range, levels = range_lab[specs])
RANGE_COL <- setNames(c("#9CC2DC", PAL$blue, "#084571"), range_lab[specs])
normal5 <- data.frame(x = seq(-1.0, 1.6, length.out = 300))
normal5$y <- dnorm(normal5$x, gridcmp$reml$mean, gridcmp$reml$sd)
studies <- gridcmp$data
p5a <- ggplot() +
  geom_ribbon(data = curves[curves$range == range_lab[["2"]], ],
              aes(x = value, ymin = lower, ymax = upper), fill = PAL$band, alpha = 0.14) +
  geom_line(data = normal5, aes(x = x, y = y), color = PAL$vermillion, linewidth = 0.5,
            linetype = "22") +
  geom_line(data = curves, aes(x = value, y = mean, color = range), linewidth = 0.65) +
  geom_rug(data = studies, aes(x = yi), sides = "b", color = PAL$ink2, length = unit(1.8, "mm"),
           linewidth = 0.35) +
  scale_color_manual(values = RANGE_COL, name = "Range: the estimates\nextended by") +
  scale_x_continuous(breaks = seq(-1, 1.5, 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(xlim = c(-1.0, 1.6)) +
  label_text(x = 1.55, y = 1.25, hjust = 1, vjust = 1, color = PAL$vermillion, size = 2.3,
             label = sprintf("Normal (REML)\nSD %.2f", num("raud_reml_sd"))) +
  labs(x = "Effect (standardized mean difference)", y = "Density of g") +
  theme_paper() +
  theme(legend.position = c(0.79, 0.80), legend.background = element_rect(fill = "white", color = NA),
        legend.key.height = unit(2.6, "mm"), legend.text = element_text(size = BASE_SIZE - 1.2),
        legend.title = element_text(size = BASE_SIZE - 1))

sites5 <- attr(gridcmp$comparison, "sites")
sites5$range <- factor(range_lab[as.character(sites5$spec)], levels = range_lab[specs])
sites5$study <- factor(paste("Study", sites5$label), levels = c("Study 4", "Study 8"))
assert(near(round(sites5$mean[sites5$label == "4" & sites5$spec == 3], 2),
            round(num("raud_site4_wide"), 2)), "fig5b: study 4 on the wide range drifted")
p5b <- ggplot(sites5, aes(y = range, color = range)) +
  geom_vline(xintercept = 0, color = PAL$ink3, linewidth = 0.3) +
  geom_segment(aes(x = lower, xend = upper, yend = range), linewidth = 0.9) +
  geom_point(aes(x = estimate), shape = 21, fill = "white", color = PAL$ink2, size = 1.8,
             stroke = 0.5) +
  geom_point(aes(x = mean), size = 2.1) +
  scale_color_manual(values = RANGE_COL, guide = "none") +
  scale_y_discrete(limits = rev(levels(sites5$range))) +
  facet_wrap(~ study, ncol = 1, scales = "free_x") +
  labs(x = "Effect", y = "Range") +
  theme_paper() + theme(panel.spacing.y = unit(3, "mm"))

share_df <- rbind(
  data.frame(case = "97 firms", z = (krw$sites$estimate - num("krw_mean_g_pp")) / num("krw_sd_g_pp"),
             share = 100 * krw$share_about_g, se = krw$sites$se / num("krw_sd_g_pp")),
  data.frame(case = "1,500 simulated sites",
             z = (twin$sites$estimate - num("twin_mean_g")) / num("twin_sd_g"),
             share = 100 * twin$sites$sd_about_g^2 / twin$sites$sd^2,
             se = twin$sites$se / num("twin_sd_g"))
)
share_df$case <- factor(share_df$case, levels = c("97 firms", "1,500 simulated sites"))
assert(near(round(max(share_df$share[share_df$case == "97 firms"])),
            round(num("krw_share_about_g_max_pct"))), "fig5c: largest share drifted")
share_panel <- function(case, legend) {
  ggplot(share_df[share_df$case == case, ], aes(x = z, y = share)) +
    geom_point(aes(fill = se), shape = 21, color = "white", size = 1.7, stroke = 0.2, alpha = 0.9) +
    scale_fill_gradient(low = "#084571", high = "#CBDEEC", name = "SE / SD of g") +
    scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 10),
                       labels = function(x) paste0(x, "%")) +
    labs(x = "Estimate, in SDs of g from its mean",
         y = "Share of posterior variance\nfrom uncertainty about g") +
    facet_wrap(~ case) +
    theme_paper() +
    theme(legend.position = if (legend) c(0.16, 0.74) else "none",
          legend.key.height = unit(2.6, "mm"), legend.key.width = unit(2.4, "mm"),
          legend.background = element_rect(fill = NA, color = NA))
}
p5c <- share_panel("97 firms", legend = TRUE)
p5d <- share_panel("1,500 simulated sites", legend = TRUE)
fig5 <- tags((p5a | p5b) / (p5c | p5d) + plot_layout(widths = c(1.15, 1)))
save_fig(fig5, "fig5_grid_and_uncertainty", MM_FULL, 112)

# ==============================================================================
# Supplement figures
# ==============================================================================
# S1: distribution functions of the estimates of g, 97 firms and 1,500 sites,
# against three empirical Bayes fits
cdf_plot <- function(cdf, scale = 1, xlab, truth = NULL) {
  be <- cdf[cdf$method == "bayesEfron", ]
  others <- cdf[cdf$method != "bayesEfron", ]
  others$method <- factor(others$method, levels = c("NPMLE", "Normal", "Unimodal"))
  p <- ggplot() +
    geom_ribbon(data = be, aes(x = x * scale, ymin = lower, ymax = upper), fill = PAL$band,
                alpha = 0.18) +
    geom_line(data = be, aes(x = x * scale, y = cdf, color = "bayesEfron"), linewidth = 0.8) +
    geom_step(data = others[others$method == "NPMLE", ], aes(x = x * scale, y = cdf, color = "NPMLE"),
              linewidth = 0.5) +
    geom_line(data = others[others$method != "NPMLE", ],
              aes(x = x * scale, y = cdf, color = method), linewidth = 0.5) +
    scale_color_manual(values = METHOD_COL[c("bayesEfron", "NPMLE", "Normal", "Unimodal")],
                       breaks = c("bayesEfron", "NPMLE", "Normal", "Unimodal"), name = NULL) +
    labs(x = xlab, y = "Distribution function") + theme_paper()
  if (!is.null(truth)) {
    p <- p + stat_ecdf(data = data.frame(x = truth), aes(x = x), geom = "step", color = PAL$ink,
                       linewidth = 0.4, linetype = "dotted")
  }
  p
}
pS1a <- cdf_plot(krw$ebnm$cdf, scale = 100, xlab = "Callback gap (percentage points)") +
  theme(legend.position = c(0.8, 0.3))
pS1b <- cdf_plot(twin$ebnm$cdf, xlab = "Effect (dotted: the true effects)", truth = sim$theta_true) +
  theme(legend.position = "none")
figS1 <- tags(pS1a | pS1b)
save_fig(figS1, "figS1_cdf_comparison", MM_FULL, 62, dir = OSMDIR)

# S2: estimates and standard errors of the two data sets
pS2a <- ggplot(krw$data, aes(x = yi * 100, y = sei * 100)) +
  geom_vline(xintercept = 0, color = PAL$ink3, linewidth = 0.3) +
  geom_point(color = PAL$blue, size = 1.5, alpha = 0.8) +
  labs(x = "Estimate, callback gap (percentage points)", y = "Standard error (percentage points)") +
  theme_paper()
pS2b <- ggplot(sim, aes(x = yi, y = sei, color = ifelse(theta_true < 0, "lower group", "upper group"))) +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_color_manual(values = c("lower group" = PAL$vermillion, "upper group" = PAL$blue), name = NULL) +
  labs(x = "Estimate", y = "Standard error") + theme_paper() +
  theme(legend.position = c(0.85, 0.2)) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
figS2 <- tags(pS2a | pS2b)
save_fig(figS2, "figS2_data", MM_FULL, 60, dir = OSMDIR)

message("all figures written")
