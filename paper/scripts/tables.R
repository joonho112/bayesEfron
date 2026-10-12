# ==============================================================================
# tables.R: the LaTeX tables whose numbers come from data/ (so that no number
# is typed by hand). Run from the paper folder root:  Rscript scripts/tables.R
#
#   tables/tab-ebnm.tex       Table 3 of the article: bayesEfron beside
#                             empirical Bayes fits of the same data
#   tables/tab-krw-ebnm-full.tex, tables/tab-twin-ebnm-full.tex
#                             the shipped comparison with all ebnm families (OSM)
#   tables/tab-coverage.tex   the simulation designs and coverage (OSM)
#   tables/tab-convergence.tex  the three settings on the 20 large panels (OSM)
#   tables/tab-few.tex        the 19-study simulation (OSM)
# ==============================================================================
.libPaths(c(file.path(normalizePath("."), "r-library"), .libPaths()))
print(packageVersion("bayesEfron"))
stopifnot(packageVersion("bayesEfron") == "0.4.0")
dir.create("tables", showWarnings = FALSE)
krw <- readRDS("data/krw.rds")
twin <- readRDS("data/twin.rds")
ev <- readRDS("data/evidence.rds")
registry <- read.csv("data/numbers.csv", stringsAsFactors = FALSE)
regnum <- function(name) as.numeric(registry$value[registry$name == name])
f <- function(x, d = 3) formatC(x, format = "f", digits = d)
pct <- function(x, d = 0) paste0(formatC(100 * x, format = "f", digits = d), "\\%")
tex_escape <- function(x) gsub("%", "\\\\%", x, fixed = FALSE)

# ---- Table 3 (main text) -----------------------------------------------------
k <- krw$ebnm$table; rownames(k) <- k$method
t <- twin$ebnm$table; rownames(t) <- t$method
order_methods <- c("Normal", "Unimodal", "NPMLE", "Plug-in g", "bayesEfron")
label <- c("Normal" = "Normal (ebnm)", "Unimodal" = "Unimodal (ebnm)", "NPMLE" = "NPMLE (ebnm)",
           "Plug-in g" = "bayesEfron $g$, plugged in", "bayesEfron" = "bayesEfron")
rows <- vapply(order_methods, function(m) {
  paste0(label[[m]], " & ", f(k[m, "correlation"], 3), " & ", f(k[m, "length_ratio"], 3),
         " & ", f(t[m, "rmse"], 3), " & ", pct(t[m, "inclusion"], 1), " & ", f(t[m, "mean_length"], 2),
         " \\\\")
}, character(1))
tab3 <- c(
  "\\begin{table}[!t]", "\\centering", "\\footnotesize",
  "\\begin{tabular}{@{}l rr rrr@{}}", "\\toprule",
  " & \\multicolumn{2}{c}{97 firms} & \\multicolumn{3}{c}{1,500 simulated sites} \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-6}",
  "Method & Correlation & Length ratio & RMSE & Inclusion & Length \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}",
  paste0("\\caption{Comparison of bayesEfron with empirical Bayes fits to the same data using ebnm ",
         "(version ", ev$about$versions[["ebnm"]], "; locations estimated where applicable; 90\\% ",
         "equal-tailed intervals). For the 97 firms, the correlation of a method's posterior means ",
         "with those of bayesEfron and the length of its intervals relative to bayesEfron; for the ",
         "1,500 simulated sites, the root mean squared error of the posterior means, the inclusion ",
         "proportion and mean interval length. ",
         "The plug-in row treats the posterior mean of $g$ from bayesEfron as known. The last two ",
         "rows therefore differ only in accounting for uncertainty in $g$. The REML normal model ",
         "of Section 3.2 gives similar results to the ebnm normal prior (inclusion ", f(regnum("twin_normal_inside_90_pct"), 1), "\\%, RMSE ",
         f(regnum("twin_rmse_normal"), 3), ").}"),
  "\\label{tab:ebnm}", "\\end{table}"
)
writeLines(tab3, "tables/tab-ebnm.tex")

# ---- OSM: the shipped comparison with all ebnm families ----------------------
sk <- ev$stored_krw
rows <- apply(sk, 1, function(r) {
  paste0(r[["method"]], " & ", f(as.numeric(r[["correlation"]]), 3), " & ",
         f(as.numeric(r[["largest_difference_in_se"]]), 2), " & ",
         f(as.numeric(r[["sd_ratio"]]), 3), " & ", f(as.numeric(r[["length_ratio"]]), 3), " & ",
         ifelse(is.na(r[["loglik"]]), "--", f(as.numeric(r[["loglik"]]), 1)), " \\\\")
})
writeLines(c("\\begin{tabular}{@{}l rrrrr@{}}", "\\toprule",
             "Method & Correlation & Largest difference (SE) & SD ratio & Length ratio & Log likelihood \\\\",
             "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), "tables/tab-krw-ebnm-full.tex")
st <- ev$stored_twin
st$method[st$method == "Known g"] <- "Empirical reference $g$"
rows <- apply(st, 1, function(r) {
  paste0(r[["method"]], " & ", f(as.numeric(r[["rmse"]]), 3), " & ",
         pct(as.numeric(r[["inclusion"]]), 1), " (", f(100 * as.numeric(r[["inclusion_se"]]), 1), ") & ",
         f(as.numeric(r[["mean_length"]]), 2), " & ",
         ifelse(is.na(r[["w1"]]), "--", f(as.numeric(r[["w1"]]), 3)), " \\\\")
})
writeLines(c("\\begin{tabular}{@{}l rrrr@{}}", "\\toprule",
             "Method & RMSE & Inclusion (ref. SE) & Mean length & $W_1$ to the truth \\\\",
             "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), "tables/tab-twin-ebnm-full.tex")

# ---- OSM: coverage over simulated data sets ----------------------------------
ev4 <- read.csv("data/figure4-table.csv", stringsAsFactors = FALSE)
ev4$design <- gsub("\n", " ", ev4$design, fixed = TRUE)
design_order <- c("Two groups, K = 100 (100 data sets)", "Two groups, K = 1,500 (20 data sets)",
                  "19 studies, normal SD 0.14 (200 data sets)", "19 studies, normal SD 0.31 (200 data sets)",
                  "19 studies, fitted g (200 data sets)", "Two groups, equal SE, K = 1,500 (20 data sets)")
method_order <- c("bayesEfron", "bayesEfron g, plugged in", "Normal (plug-in)", "NPMLE",
                  "Efron's estimator", "Known g")
ev4 <- ev4[order(match(ev4$design, design_order), match(ev4$method, method_order)), ]
rows <- character()
for (d in design_order) {
  rows <- c(rows, paste0("\\multicolumn{4}{@{}l}{\\emph{", sub("fitted g", "fitted $g$", d, fixed = TRUE), "}} \\\\"))
  sub <- ev4[ev4$design == d, ]
  for (i in seq_len(nrow(sub))) {
    m <- sub$method[i]
    m <- sub("bayesEfron g, plugged in", "bayesEfron $g$, plugged in", m, fixed = TRUE)
    m <- sub("Known g", "Known $g$", m, fixed = TRUE)
    rows <- c(rows, paste0("\\quad ", m, " & ", f(sub$inclusion[i], 3), " (", f(sub$inclusion_se[i], 3),
                           ") & ", f(sub$length[i], 3), " (", f(sub$length_se[i], 3), ") & ",
                           f(sub$rmse[i], 3), " (", f(sub$rmse_se[i], 3), ") \\\\"))
  }
}
writeLines(c("\\begin{tabular}{@{}l rrr@{}}", "\\toprule",
             "Design and method & Coverage (MC SE) & Mean length (MC SE) & RMSE (MC SE) \\\\",
             "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), "tables/tab-coverage.tex")

# ---- OSM: convergence of the three settings ---------------------------------
cv <- ev$panels_convergence
rows <- apply(cv, 1, function(r) {
  paste0(gsub('"', "", r[["setting"]]), " & ", r[["valid"]], " & ", f(as.numeric(r[["rhat_max"]]), 3),
         " & ", round(as.numeric(r[["ess_bulk_min"]])), " & ", f(as.numeric(r[["treedepth_share_max"]]), 3),
         " & ", round(as.numeric(r[["sampling_seconds_median"]])), " & ",
         round(as.numeric(r[["pilot_seconds_median"]])), " & ", f(as.numeric(r[["inclusion"]]), 3), " \\\\")
})
writeLines(c("\\begin{tabular}{@{}l rrrrrrr@{}}", "\\toprule",
             "Setting & \\shortstack{Valid\\\\of 20} & \\shortstack{Largest\\\\R-hat} & \\shortstack{Smallest\\\\bulk ESS} & \\shortstack{Tree depth\\\\share, max} & \\shortstack{Sampling\\\\s} & \\shortstack{Pilot\\\\s} & Coverage \\\\",
             "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), "tables/tab-convergence.tex")

# ---- OSM: the 19-study simulation --------------------------------------------
few <- ev$few
few <- few[few$method != "bayesEfron, adapt_delta = 0.99", ]
method_lab <- c("bayesEfron" = "bayesEfron", "normal" = "Normal, REML plug-in", "known g" = "Known $g$")
rows <- character()
for (sc in unique(few$scenario)) {
  sub <- few[few$scenario == sc, ]
  rows <- c(rows, paste0("\\multicolumn{8}{@{}l}{\\emph{", sub("fitted g", "fitted $g$", sc), ", true SD ",
                         f(sub$true_sd[1], 2), "}} \\\\"))
  for (i in seq_len(nrow(sub))) {
    rows <- c(rows, paste0("\\quad ", method_lab[[sub$method[i]]], " & ", f(sub$sd_mean[i], 2), " & ",
                           ifelse(is.na(sub$sd_inclusion[i]), "--", f(sub$sd_inclusion[i], 2)), " & ",
                           f(sub$inclusion[i], 3), " (", f(sub$inclusion_se[i], 3), ") & ",
                           f(sub$length[i], 2), " & ", f(sub$rmse[i], 3), " & ",
                           f(sub$largest_inclusion[i], 3), " & ", f(sub$largest_error[i], 2), " \\\\"))
  }
}
writeLines(c("\\begin{tabular}{@{}l rrrrrrr@{}}", "\\toprule",
             "\\shortstack{Distribution\\\\and method} & \\shortstack{Mean\\\\SD} & \\shortstack{SD interval\\\\covers} & \\shortstack{Coverage\\\\(MC SE)} & Length & RMSE & \\shortstack{Largest:\\\\coverage} & \\shortstack{Largest:\\\\bias} \\\\",
             "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), "tables/tab-few.tex")
message("tables written")

# ---- OSM: the three ranges for the 19 studies --------------------------------
# Put the three specifications in columns to avoid shrinking the type.
# Keep the previous formatting and precision before transposing the fields.
g <- readRDS("data/grid-comparison.rds")
tg <- g$table[order(g$table$expansion), ]
formatted <- apply(tg, 1, function(r) {
  paste0(f(as.numeric(r[["expansion"]]), 2), " & ", f(as.numeric(r[["lower"]]), 2), " to ",
         f(as.numeric(r[["upper"]]), 2), " & ", f(as.numeric(r[["elpd_loo"]]), 2), " (",
         f(as.numeric(r[["se_elpd_loo"]]), 2), ") & ", f(as.numeric(r[["elpd_diff"]]), 2), " (",
         f(as.numeric(r[["se_diff"]]), 2), ") & ", f(as.numeric(r[["mean_g"]]), 3), " & ",
         f(as.numeric(r[["sd_g"]]), 3), " & ", f(as.numeric(r[["site_4"]]), 3), " & ",
         f(as.numeric(r[["site_8"]]), 3), " & ", f(as.numeric(r[["max_rhat"]]), 3), " & ",
         as.integer(r[["divergences"]]), " & ", f(as.numeric(r[["end_mass"]]), 3), " & ",
         ifelse(as.logical(r[["adequate"]]), "yes", "no"), " \\\\")
})

cells <- vapply(strsplit(formatted, " & ", fixed = TRUE), function(r) {
  r[length(r)] <- trimws(substr(r[length(r)], 1, nchar(r[length(r)]) - 2L))
  r
}, character(12))
grid_row <- function(label, values) {
  paste0(paste(c(label, values), collapse = " & "), " \\\\")
}
rows <- c(
  "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lrrr@{}}",
  "\\toprule",
  " & \\multicolumn{3}{c}{Extension} \\\\",
  "\\cmidrule(l){2-4}",
  grid_row("", cells[1, ]),
  "\\midrule",
  grid_row("Range", cells[2, ]),
  "\\addlinespace",
  grid_row("elpd (SE)", cells[3, ]),
  grid_row("Difference (SE)", cells[4, ]),
  "\\addlinespace",
  grid_row("Mean of $g$", cells[5, ]),
  grid_row("SD of $g$", cells[6, ]),
  grid_row("Study 4", cells[7, ]),
  grid_row("Study 8", cells[8, ]),
  "\\addlinespace",
  grid_row("Largest R-hat", cells[9, ]),
  grid_row("Divergences", cells[10, ]),
  grid_row("End mass", cells[11, ]),
  grid_row("Resolution adequate", cells[12, ]),
  "\\bottomrule", "\\end{tabular*}"
)
writeLines(rows, "tables/tab-grid.tex")
message("grid table written")


# ---- Table E.2: coverage of the intervals for subsets of the sites ------------
# Stack the two K values and use full-width subset headers. The four body
# columns remain legible without scaling the table down to fit the page.
tc <- read.csv("data/tail-coverage.csv", stringsAsFactors = FALSE)
subset_label <- c("All sites" = "All sites",
                  "Estimate between the groups (-0.7, 0.7)" = "Estimate between the groups",
                  "Estimates outside the sample 2.5th and 97.5th percentiles" =
                    "Estimates outside the sample 2.5th and 97.5th percentiles",
                  "True effect in the outer tenth at either end of its group" =
                    "True effect in the outer tenth at either end of its group")
stopifnot(all(unique(tc$subset) %in% names(subset_label)))
method_label <- c("bayesEfron" = "bayesEfron", "bayesEfron g, plugged in" = "bayesEfron $g$, plugged in",
                  "True g, plugged in" = "True $g$, plugged in", "Normal (ebnm)" = "Normal (ebnm)",
                  "NPMLE (ebnm)" = "NPMLE (ebnm)")
stopifnot(all(unique(tc$method) %in% names(method_label)))
rows <- c("\\begin{tabular}{@{}l rrr@{}}", "\\toprule",
          "Method & Sites per set & Coverage & MC SE \\\\", "\\midrule")
for (K in c(100L, 1500L)) {
  nsets <- unique(tc$panels[tc$K == K])
  stopifnot(length(nsets) == 1L)
  if (K == 1500L) rows <- c(rows, "\\midrule")
  rows <- c(rows, paste0("\\multicolumn{4}{@{}l}{\\textbf{$K = ",
                         if (K == 1500L) "1{,}500" else "100", "$; ",
                         nsets, " data sets}} \\\\"))
  for (s in names(subset_label)) {
    rows <- c(rows, paste0("\\multicolumn{4}{@{}l}{\\emph{", subset_label[[s]], "}} \\\\"))
    for (m in names(method_label)) {
      r <- tc[tc$K == K & tc$method == m & tc$subset == s, ]
      stopifnot(nrow(r) == 1L)
      rows <- c(rows, paste0("\\quad ", method_label[[m]], " & ", f(r$sites_per_panel, 1), " & ",
                            f(r$inclusion, 3), " & ", f(r$inclusion_se, 3), " \\\\"))
    }
  }
}
rows <- c(rows, "\\bottomrule", "\\end{tabular}")
writeLines(rows, "tables/tab-tail.tex")
cat("subset table written\n")

# Format the supplementary tables without changing their numerical values.
source("scripts/table-layout.R")
