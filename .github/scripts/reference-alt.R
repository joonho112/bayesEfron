#!/usr/bin/env Rscript
# Add descriptions to the plots of the reference pages and repair the script
# of the 404 page, after pkgdown has built the site.
# Usage, from the package root: Rscript .github/scripts/reference-alt.R [docs]
#
# The list below has one description for each plot that the examples of a
# help page draw. Change it when those examples change. A plot without a
# description, or a description without a plot, gives a warning and does not
# stop the build: the site is complete without the descriptions.
args <- commandArgs(trailingOnly = TRUE)
site <- if (length(args)) args[[1L]] else "docs"

problems <- 0L
problem <- function(...) {
  text <- paste0(...)
  problems <<- problems + 1L
  if (identical(Sys.getenv("GITHUB_ACTIONS"), "true")) {
    # Shown as an annotation of the workflow run.
    cat("::warning::", text, "\n", sep = "")
  } else {
    message("Warning: ", text)
  }
}

if (!requireNamespace("xml2", quietly = TRUE)) {
  problem("The package xml2 is not installed; the site is left as pkgdown built it.")
  quit(save = "no", status = 0L)
}

# Texts that several pages share.
extremes <- paste(
  "Two panels with the ten firms with the highest and the ten with the",
  "lowest posterior mean callback gaps, in percentage points: posterior",
  "means, 50% and 90% intervals, the estimates as open circles, and a",
  "dashed line at the mean of the distribution of effects."
)
five_panels <- paste(
  "Five panels for the fit of the callback data: R-hat and the bulk and",
  "tail effective sample sizes of the monitored quantities in blocks,",
  "against their thresholds; the largest mass that a site has at each end",
  "of the grid; and the standard deviation of the distribution of effects",
  "under the prior and under the posterior."
)
callback_grid <- paste(
  "The six spline basis functions and the knots of the grid for the",
  "callback data, over the 97 estimates in increasing order, each with two",
  "standard errors on each side, in percentage points."
)
three_ranges <- paste(
  "Three rows, one for each range of the grid for the 19 studies. Each row",
  "has the estimated density with its band, the mean and the standard",
  "deviation of the distribution of effects, and the posterior means and",
  "intervals of studies 4 and 8. The mean of study 4 falls as the range",
  "widens."
)

plots <- list(
  "bayes_efron_fit.html" = c(
    "bayes_efron_fit-1.png" = paste(
      "The estimated density of the callback gaps of the 97 firms, with a",
      "pointwise 90% band."
    )
  ),
  "make_efron_grid.html" = c("make_efron_grid-1.png" = callback_grid),
  "plot.bef_fit_re.html" = c(
    "plot.bef_fit_re-1.png" = extremes,
    "plot.bef_fit_re-2.png" = paste(
      "The estimated density of the callback gaps in percentage points, with",
      "a pointwise 90% band."
    ),
    "plot.bef_fit_re-3.png" = five_panels,
    "plot.bef_fit_re-4.png" = paste(
      "The posterior mean of the gap of each firm against its estimate, in",
      "percentage points, colored by the standard error, with a dashed",
      "diagonal. The horizontal axis is labeled Estimated gap (percentage",
      "points)."
    )
  ),
  "plot_g.html" = c(
    "plot_g-1.png" = paste(
      "The estimated density of the callback gaps in percentage points with",
      "a pointwise 90% band, over a histogram of the 97 estimates, and the",
      "normal distribution of the random-effects model as a dashed curve."
    ),
    "plot_g-2.png" = paste(
      "The estimated distribution function of the callback gaps, with a 50%",
      "band inside a 90% band."
    ),
    "plot_g-3.png" = paste(
      "The estimated density for the 1,500 simulated sites with its band,",
      "over a histogram of the true effects, which has two groups, and a",
      "dashed normal curve, which has one."
    ),
    "plot_g-4.png" = paste(
      "The estimated density for the 19 teacher expectancy studies with its",
      "band, and forty thin curves, each a posterior draw of the",
      "distribution."
    )
  ),
  "plot_sites.html" = c(
    "plot_sites-1.png" = extremes,
    "plot_sites-2.png" = paste(
      "The 19 studies in the order of their standard errors, the most",
      "precise at the top: posterior means, 50% and 90% intervals, the",
      "estimates as open circles, and a dashed line at zero."
    ),
    "plot_sites-3.png" = paste(
      "The 90% posterior intervals of the 1,500 simulated sites as thin",
      "lines without labels, in the order of the posterior means, with a",
      "dashed line at the mean of the distribution of effects."
    )
  ),
  "plot_shrinkage.html" = c(
    "plot_shrinkage-1.png" = paste(
      "The posterior mean of the callback gap of each firm against its",
      "estimate, in percentage points, colored by the standard error. The",
      "points lie closer to the center than the dashed diagonal, those of",
      "the least precise estimates the most."
    ),
    "plot_shrinkage-2.png" = paste(
      "For each of the 19 studies an arrow from the estimate to the",
      "posterior mean, at the height of its standard error. The arrows of",
      "the imprecise studies are the longest."
    ),
    "plot_shrinkage-3.png" = paste(
      "Two panels for the 1,500 simulated sites, this fit and a normal",
      "model: posterior means against estimates, colored by the true effect.",
      "The fit moves the estimates toward two groups and the normal model",
      "toward one mean."
    )
  ),
  "plot_exceedance.html" = c(
    "plot_exceedance-1.png" = paste(
      "The probability that the callback gap of each of the 97 firms is",
      "positive, the largest at the top, with a dashed line at 0.95. The",
      "firms at or above it are filled points."
    ),
    "plot_exceedance-2.png" = paste(
      "The probability that the gap of each firm exceeds two percentage",
      "points, the largest at the top."
    ),
    "plot_exceedance-3.png" = paste(
      "The probabilities that the effects of the 1,500 simulated sites are",
      "positive, in increasing order, as a curve from 0 to 1, with a dashed",
      "line at 0.9."
    )
  ),
  "plot_site_posterior.html" = c(
    "plot_site_posterior-1.png" = paste(
      "Three panels, each with the posterior distribution of the effect of a",
      "simulated site, the likelihood of its estimate and the distribution",
      "of effects. For two imprecise estimates between the groups the",
      "posterior distribution has two modes; for the most precise estimate",
      "it is a narrow curve at the estimate."
    ),
    "plot_site_posterior-2.png" = paste(
      "Two panels for the sites that lean most on the ends of the grid: the",
      "posterior distribution of the effect, the likelihood of the estimate",
      "and the distribution of effects."
    )
  ),
  "plot_diagnostics.html" = c(
    "plot_diagnostics-1.png" = five_panels,
    "plot_diagnostics-2.png" = paste(
      "One panel with the R-hat of all monitored quantities of the fit of",
      "the simulated sites as a box, against a dashed line at 1.01 and a",
      "dotted line at 1.05, with the worst quantity named."
    )
  ),
  "plot_loo.html" = c(
    "plot_loo-1.png" = paste(
      "The Pareto k value of each of the 97 firms, all below the dashed line",
      "at 0.656."
    ),
    "plot_loo-2.png" = paste(
      "The Pareto k value of each of the 19 studies, all below the dashed",
      "line at 0.656."
    )
  ),
  "plot_ppc.html" = c(
    "plot_ppc-1.png" = paste(
      "The empirical distribution function of the PIT values of the 97",
      "estimates less the uniform one, as a step line inside a gray band",
      "around zero."
    ),
    "plot_ppc-2.png" = paste(
      "A histogram of the PIT values of the 1,500 estimates, with a dashed",
      "line at the height that a uniform distribution gives every bar and a",
      "gray band around it."
    ),
    "plot_ppc-3.png" = paste(
      "A histogram of the 97 estimates in percentage points, with the",
      "density that the fitted model gives to an estimate as a curve."
    )
  ),
  "plot_mcmc.html" = c(
    "plot_mcmc-1.png" = paste(
      "Trace plots of lambda and of the mean and the standard deviation of",
      "the distribution of effects, four chains each."
    ),
    "plot_mcmc-2.png" = paste(
      "Rank plots of the six spline coefficients, one line for each chain."
    ),
    "plot_mcmc-3.png" = paste(
      "Densities of the mean and of the standard deviation of the",
      "distribution of effects, one line for each chain."
    )
  ),
  "plot.bef_grid.html" = c(
    "plot.bef_grid-1.png" = callback_grid,
    "plot.bef_grid-2.png" = paste(
      "The eight spline basis functions and the knots of the grid of the",
      "fit of the simulated sites, over a histogram of the 1,500 estimates,",
      "some of which lie outside the grid."
    )
  ),
  "plot.bef_grid_comparison.html" = c("plot.bef_grid_comparison-1.png" = three_ranges),
  "compare_efron_grids.html" = c("compare_efron_grids-1.png" = three_ranges),
  "plot_recovery.html" = c(
    "plot_recovery-1.png" = paste(
      "The 90% posterior intervals of 150 of the 1,500 simulated sites in",
      "the order of their posterior means, with the true effects as crosses,",
      "and as orange triangles where they lie outside the interval."
    ),
    "plot_recovery-2.png" = paste(
      "The share of the true effects inside the posterior interval against",
      "the level of the interval, from 5% to 95%, close to the dashed",
      "diagonal."
    ),
    "plot_recovery-3.png" = paste(
      "The posterior means of the 1,500 simulated sites against their true",
      "effects, around the dashed diagonal."
    )
  )
)

verified <- changed <- 0L
for (page in names(plots)) {
  file <- file.path(site, "reference", page)
  if (!file.exists(file)) {
    problem("The reference page ", file, " does not exist.")
    next
  }
  text <- paste(readLines(file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  positions <- gregexpr("<img\\b[^>]*>", text, perl = TRUE)
  tags <- regmatches(text, positions)[[1L]]
  seen <- character()
  revised <- vapply(tags, function(tag) {
    node <- xml2::xml_find_first(xml2::read_html(tag), "//img")
    src <- basename(xml2::xml_attr(node, "src"))
    descriptions <- plots[[page]]
    if (!src %in% names(descriptions)) return(tag)
    seen <<- c(seen, src)
    description <- unname(descriptions[[src]])
    verified <<- verified + 1L
    if (identical(xml2::xml_attr(node, "alt"), description)) return(tag)
    changed <<- changed + 1L
    attribute <- paste0('alt="', description, '"')
    pattern <- "\\balt\\s*=\\s*(\"[^\"]*\"|'[^']*')"
    if (grepl(pattern, tag, perl = TRUE)) {
      sub(pattern, attribute, tag, perl = TRUE)
    } else {
      sub("<img\\b", paste("<img", attribute), tag, perl = TRUE)
    }
  }, character(1L), USE.NAMES = FALSE)
  if (anyDuplicated(seen) || !setequal(seen, names(plots[[page]]))) {
    problem(
      "The plots of ", page, " are not those listed in ",
      ".github/scripts/reference-alt.R; update the descriptions."
    )
  }
  regmatches(text, positions) <- list(revised)
  if (!identical(tags, revised)) writeLines(text, file, useBytes = TRUE)
}

# A plot that an example has gained since the list was written.
for (file in list.files(file.path(site, "reference"), "[.]html$", full.names = TRUE)) {
  nodes <- xml2::xml_find_all(xml2::read_html(file), "//img")
  src <- xml2::xml_attr(nodes, "src")
  alt <- xml2::xml_attr(nodes, "alt")
  informative <- !basename(src) %in% c("logo.png", "logo.svg")
  if (any(informative & (is.na(alt) | !nzchar(trimws(alt))))) {
    problem("An image in ", file, " has no description.")
  }
}
cat(sprintf("Descriptions of reference plots: %d checked, %d written.\n", verified, changed))

# When pkgdown makes the links of the 404 page absolute, some versions give
# the inline script that configures MathJax a source ending in /NA. A script
# with a source does not run its own content, so the configuration is lost.
# Only that attribute is removed.
file <- file.path(site, "404.html")
if (!file.exists(file)) {
  problem("The page ", file, " does not exist.")
} else {
  text <- paste(readLines(file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  positions <- gregexpr("(?s)<script\\b[^>]*>.*?</script>", text, perl = TRUE)
  tags <- regmatches(text, positions)[[1L]]
  repaired <- 0L
  revised <- vapply(tags, function(tag) {
    node <- xml2::xml_find_first(xml2::read_html(tag), "//script")
    src <- xml2::xml_attr(node, "src")
    if (is.na(src) || !grepl("/NA$", src) ||
        !grepl("window.MathJax", xml2::xml_text(node), fixed = TRUE)) return(tag)
    repaired <<- repaired + 1L
    sub("\\s+src\\s*=\\s*(\"[^\"]*\"|'[^']*')", "", tag, perl = TRUE)
  }, character(1L), USE.NAMES = FALSE)
  regmatches(text, positions) <- list(revised)
  if (!identical(tags, revised)) writeLines(text, file, useBytes = TRUE)
  nodes <- xml2::xml_find_all(xml2::read_html(file), "//script[@src]")
  if (any(grepl("/NA$", xml2::xml_attr(nodes, "src")))) {
    problem("A script of ", file, " has a source ending in /NA that was not repaired.")
  }
  cat(sprintf("404 page: %d script source removed.\n", repaired))
}

if (problems > 0L) {
  cat(sprintf("%d warning(s); the site was not changed where they apply.\n", problems))
}
