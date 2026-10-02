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

plots <- list(
  "bayes_efron_fit.html" = c(
    "bayes_efron_fit-1.png" = paste(
      "The posterior mean of the probability that the distribution of",
      "effects gives to each grid point, with a pointwise 90% band, for the",
      "teacher expectancy studies."
    )
  ),
  "make_efron_grid.html" = c(
    "make_efron_grid-1.png" = paste(
      "The six natural cubic spline functions on the default grid for the",
      "teacher expectancy studies."
    )
  ),
  "plot.bef_fit_re.html" = c(
    "plot.bef_fit_re-1.png" = paste(
      "The posterior mean of the effect of each study with its 90% interval",
      "and its central 50% interval, in the order of the posterior means. A",
      "dashed line marks the mean of the distribution of effects."
    ),
    "plot.bef_fit_re-2.png" = paste(
      "The same plot with 95% intervals and the studies in the order of",
      "their standard errors."
    ),
    "plot.bef_fit_re-3.png" = paste(
      "The posterior mean of the probability that the distribution of",
      "effects gives to each grid point, with a pointwise 90% band."
    ),
    "plot.bef_fit_re-4.png" = paste(
      "Four panels with the largest R-hat, the smallest bulk and tail",
      "effective sample sizes, the number of divergent transitions and the",
      "number of iterations at the maximum tree depth."
    ),
    "plot.bef_fit_re-5.png" = paste(
      "The plot of the distribution of effects with the title Teacher",
      "expectancy effects."
    ),
    "plot.bef_fit_re-6.png" = paste(
      "The Pareto k value of each study from leave-one-out",
      "cross-validation, with reference lines at 0.5 and 0.7."
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
