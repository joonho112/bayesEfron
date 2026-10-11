# Builds `raudenbush1985`, the teacher expectancy studies, in data/. Run from
# the package root:
#
#   Rscript data-raw/raudenbush.R
#
# The model fitted to these data, `raudenbush_fit`, is built by
# data-raw/example-fits.R together with the other two stored fits.

devtools::load_all(".", quiet = TRUE)

# ---- the data ----------------------------------------------------------------
# Raudenbush and Bryk (1985, Table 1); `setting` and `tester` from Raudenbush
# (1984). The same data are distributed as `dat.raudenbush1985` in the metadat
# package; if it is installed, the script compares the values with it.
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

# ---- save --------------------------------------------------------------------
dir.create("data", showWarnings = FALSE)
save(raudenbush1985, file = "data/raudenbush1985.rda", compress = "xz", version = 3)
message(sprintf("data/raudenbush1985.rda %d bytes", file.size("data/raudenbush1985.rda")))
