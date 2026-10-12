# Layout adjustments for the OSM tables.
# Preserve the already formatted values and precision.
table_lines <- function(name) readLines(file.path("tables", paste0(name, ".tex")))
table_body <- function(x) {
  x[(match("\\midrule", x) + 1L):(match("\\bottomrule", x) - 1L)]
}
table_cells <- function(x) {
  x <- x[grepl(" & ", x, fixed = TRUE)]
  x <- sub("^[[:space:]]*\\\\quad[[:space:]]*", "", x)
  x <- vapply(x, function(z) trimws(substr(z, 1L, nchar(z) - 2L)), character(1))
  do.call(rbind, strsplit(x, " & ", fixed = TRUE))
}
tex_row <- function(x, keep = FALSE) {
  paste0(paste(x, collapse = " & "), " \\\\", if (keep) "*" else "")
}
write_table <- function(name, x) writeLines(x, file.path("tables", paste0(name, ".tex")))

# D.1 and D.2: wrap long method names and column headings.
for (name in c("tab-krw-ebnm-full", "tab-twin-ebnm-full")) {
  x <- table_lines(name)
  if (name == "tab-krw-ebnm-full") {
    x[1] <- "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}L{0.28\\linewidth}rrrrr@{}}"
    x[3] <- tex_row(c("Method", "Correlation", "\\shortstack{Largest\\\\difference\\\\(SE)}",
                      "\\shortstack{SD\\\\ratio}", "\\shortstack{Length\\\\ratio}",
                      "\\shortstack{Log\\\\likelihood}"))
  } else {
    x[1] <- "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}L{0.37\\linewidth}rrrr@{}}"
    x[3] <- tex_row(c("Method", "RMSE", "\\shortstack{Inclusion\\\\(ref. SE)}",
                      "\\shortstack{Mean\\\\length}", "\\shortstack{$W_1$ to\\\\the truth}"))
  }
  x[length(x)] <- "\\end{tabular*}"
  write_table(name, x)
}

# E.1: longtable wrappers and repeated headings live in appendix_E.tex.
# Keep each design together while allowing page breaks between designs.
x <- table_body(table_lines("tab-coverage"))
groups <- which(startsWith(x, "\\multicolumn"))
for (i in seq_along(x)) {
  if (startsWith(x[i], "\\multicolumn")) {
    x[i] <- sub("{@{}l}", "{@{}p{\\linewidth}@{}}", x[i], fixed = TRUE)
  }
  if (i < length(x) && !(i + 1L) %in% groups) x[i] <- paste0(x[i], "*")
}
for (i in rev(groups[-1])) x <- append(x, "\\addlinespace[5pt]", after = i - 1L)
write_table("tab-coverage", x)

# E.2: pair the sample sizes; retain every sites-per-set, coverage and SE value.
x <- table_body(table_lines("tab-tail"))
values <- table_cells(x)
stopifnot(nrow(values) == 40L, ncol(values) == 4L)
labels <- x[grepl("\\emph{", x, fixed = TRUE)][1:4]
rows <- character()
for (j in 1:4) {
  if (j > 1L) rows <- c(rows, "\\addlinespace[5pt]")
  heading <- sub("\\multicolumn{4}{@{}l}", "\\multicolumn{5}{@{}p{\\linewidth}@{}}",
                 labels[j], fixed = TRUE)
  rows <- c(rows, paste0(heading, "*"))
  for (i in ((j - 1L) * 5L + 1L):(j * 5L)) {
    a <- values[i, ]; b <- values[i + 20L, ]
    stopifnot(a[1] == b[1])
    rows <- c(rows, tex_row(c(a[1], a[2], paste0(a[3], " (", a[4], ")"),
                             b[2], paste0(b[3], " (", b[4], ")")), keep = i %% 5L != 0L))
  }
}
write_table("tab-tail", rows)

# E.3: methods are columns within each generating distribution.
x <- table_body(table_lines("tab-few"))
values <- table_cells(x)
labels <- x[startsWith(x, "\\multicolumn")]
stopifnot(nrow(values) == 9L, ncol(values) == 8L, length(labels) == 3L)
metrics <- c("Mean SD", "SD interval covers", "Coverage (MC SE)", "Length", "RMSE",
             "Largest: coverage", "Largest: bias")
rows <- character()
for (j in 1:3) {
  if (j > 1L) rows <- c(rows, "\\addlinespace[5pt]")
  heading <- sub("\\multicolumn{8}{@{}l}", "\\multicolumn{4}{@{}p{\\linewidth}@{}}",
                 labels[j], fixed = TRUE)
  rows <- c(rows, paste0(heading, "*"))
  block <- values[((j - 1L) * 3L + 1L):(j * 3L), , drop = FALSE]
  stopifnot(identical(unname(block[, 1]), c("bayesEfron", "Normal, REML plug-in", "Known $g$")))
  for (k in seq_along(metrics)) rows <- c(rows, tex_row(c(metrics[k], block[, k + 1L]), k < 7L))
}
write_table("tab-few", rows)

# E.4: settings are columns, with diagnostics and computing times in rows.
x <- table_cells(table_body(table_lines("tab-convergence")))
stopifnot(nrow(x) == 3L, ncol(x) == 8L)
metrics <- c("Valid of 20", "Largest R-hat", "Smallest bulk ESS", "Tree depth share, max",
             "Sampling (s)", "Pilot (s)", "Coverage")
rows <- c("\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lrrr@{}}", "\\toprule",
          tex_row(c("", "defaults", "M = 8", "\\shortstack{settings =\\\\0.3.0}")), "\\midrule")
for (k in seq_along(metrics)) rows <- c(rows, tex_row(c(metrics[k], x[, k + 1L])))
write_table("tab-convergence", c(rows, "\\bottomrule", "\\end{tabular*}"))
message("OSM table layouts updated")
