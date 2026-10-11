# Builds `star_schools`, the kindergarten schools of Tennessee's Project
# STAR, in data/. Run from the package root:
#
#   Rscript data-raw/star.R <STAR_Students.tab> [audit-dir]
#
# The student file belongs to "Tennessee's Student Teacher Achievement Ratio
# (STAR) project" (Achilles et al., 2008) on Harvard Dataverse,
# https://doi.org/10.7910/DVN/SIWH9F, version 1.0, released under CC0 1.0.
# The file that `star_schools` was built from was retrieved on 2026-10-11
# (00:17 UTC) and has 13,094,524 bytes. Needs lme4 for the models and
# digest for the hashes. With `audit-dir` the script also writes the
# estimates and the counts of the sample there as CSV files.
#
# Each school is one site. Its estimate is the difference between small
# classes and regular classes (with or without a teacher's aide) in a
# composite of the reading and mathematics scores of kindergarten, in
# standard deviations of the composite. The model that gives the estimates
# and their standard errors has a random intercept for each class, because
# classes, and not single students, are what a school has few of.
#
# If the directory vignettes/articles/ exists, the script also writes
# class-size-first-stage.rds there, which the article "Small classes in 79
# schools" reads: the same contrasts under other assumptions about the
# variance of the classes.
#
# If data/star_schools.rda exists, the new estimates are compared with it.
# The script stops if they differ by more than 1e-6, or if another column
# differs. If the estimates differ by less, as they can with another version
# of lme4, it writes neither file and says so: the fit and the tables of the
# article were made from the stored numbers. Remove the file to build it
# anew.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || !file.exists(args[1L])) {
  stop("Supply the path to STAR_Students.tab from Harvard Dataverse.")
}
for (package in c("lme4", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Install ", package, ".")
}

# ---- the file ------------------------------------------------------------------
# Dataverse generates the header line of the tabular file, and four of its
# column names are placeholders (var6 to var9). The data lines are therefore
# hashed apart from the header. The estimates below were built from data
# lines with this hash.
body_md5 <- "c918b8f2233a2688396f0e5a0f4b0ce6"
bytes <- readBin(args[1L], "raw", n = file.size(args[1L]))
first_newline <- which(bytes == as.raw(10L))[1L]
file_hashes <- list(
  md5 = digest::digest(bytes, algo = "md5", serialize = FALSE),
  sha256 = digest::digest(bytes, algo = "sha256", serialize = FALSE),
  body_md5 = digest::digest(bytes[-seq_len(first_newline)], algo = "md5", serialize = FALSE),
  bytes = length(bytes)
)
rm(bytes)
if (!identical(file_hashes$body_md5, body_md5)) {
  stop("The data lines of this file differ from those that the estimates were built from.")
}

raw <- utils::read.delim(args[1L], check.names = FALSE, stringsAsFactors = FALSE,
                         na.strings = c("", "NA"))
required <- c("stdntid", "gkclasstype", "gkschid", "gksurban", "gktchid",
              "gktreadss", "gktmathss")
stopifnot(all(required %in% names(raw)), !anyDuplicated(raw$stdntid))

# ---- the kindergarten sample ---------------------------------------------------
# A student with a class type in kindergarten took part in that year. Class
# types: 1 small, 2 regular, 3 regular with a teacher's aide.
k <- raw[!is.na(raw$gkclasstype), required]
stopifnot(nrow(k) == 6325L, !anyNA(k[c("gkschid", "gktchid", "gksurban")]),
          all(k$gkclasstype %in% 1:3), length(unique(k$gkschid)) == 79L)
# A teacher ID is one class, in one school, of one type.
classes <- unique(k[c("gkschid", "gktchid", "gkclasstype")])
stopifnot(!anyDuplicated(classes$gktchid), nrow(classes) == 325L,
          identical(as.vector(table(classes$gkclasstype)), c(127L, 99L, 99L)))

complete <- k[!is.na(k$gktreadss) & !is.na(k$gktmathss), ]
stopifnot(nrow(complete) == 5786L, length(unique(complete$gkschid)) == 79L,
          length(unique(complete$gktchid)) == 325L)
standardize <- function(x) (x - mean(x)) / stats::sd(x)
complete$y <- standardize((standardize(complete$gktreadss) +
                             standardize(complete$gktmathss)) / 2)
complete$small <- as.integer(complete$gkclasstype == 1L)
complete$regular <- 1L - complete$small
complete$school <- factor(complete$gkschid)
complete$class <- factor(complete$gktchid)

# ---- estimates and standard errors ---------------------------------------------
# One intercept and one contrast for each school, and a random intercept for
# each class with one variance for all the classes, small and regular, of all
# the schools.
model <- lme4::lmer(y ~ 0 + school + school:small + (1 | class), data = complete,
                    REML = TRUE,
                    control = lme4::lmerControl(check.rankX = "stop.deficient"))
stopifnot(!lme4::isSingular(model), length(model@optinfo$conv$lme4$messages) == 0L)
coefficients <- lme4::fixef(model)
contrast <- grep(":small$", names(coefficients))
contrast_school <- sub("^school", "", sub(":small$", "", names(coefficients)[contrast]))
standard_errors <- sqrt(diag(as.matrix(stats::vcov(model))))[contrast]
components <- as.data.frame(lme4::VarCorr(model))
class_sd <- components$sdcor[components$grp == "class"]
residual_sd <- components$sdcor[components$grp == "Residual"]

# The codes of `gksurban` in the variable labels of the Dataverse file: 1
# inner city, 2 suburban, 3 rural, 4 urban.
location_labels <- c("inner city", "suburban", "rural", "urban")
schools <- lapply(levels(complete$school), function(id) {
  d <- complete[complete$school == id, ]
  stopifnot(length(unique(d$gksurban)) == 1L, any(d$small == 1L), any(d$small == 0L))
  # The same contrast with the students of a school taken as independent.
  ols <- stats::lm(y ~ small, data = d)
  data.frame(
    school_id = as.integer(id),
    location = location_labels[d$gksurban[1L]],
    yi_ols = unname(stats::coef(ols)[["small"]]),
    sei_ols = sqrt(diag(stats::vcov(ols)))[["small"]],
    n_students = nrow(d), n_students_small = sum(d$small),
    n_classes = length(unique(d$class)),
    n_classes_small = length(unique(d$class[d$small == 1L])),
    stringsAsFactors = FALSE
  )
})
schools <- do.call(rbind, schools)
schools <- schools[order(schools$school_id), ]
at <- match(as.character(schools$school_id), contrast_school)
stopifnot(!anyNA(at))
star_schools <- data.frame(
  site = seq_len(nrow(schools)), school_id = schools$school_id,
  location = schools$location,
  yi = unname(coefficients[contrast][at]), sei = unname(standard_errors[at]),
  stringsAsFactors = FALSE
)
star_schools$vi <- star_schools$sei^2
star_schools <- cbind(star_schools, schools[c("yi_ols", "sei_ols", "n_students",
                                              "n_students_small", "n_classes",
                                              "n_classes_small")])
rownames(star_schools) <- NULL
stopifnot(
  nrow(star_schools) == 79L, !anyNA(star_schools),
  all(is.finite(star_schools$yi)), all(star_schools$sei > 0),
  sum(star_schools$n_students) == 5786L, sum(star_schools$n_classes) == 325L,
  sum(star_schools$n_classes_small) == 127L,
  identical(as.vector(table(factor(star_schools$location, location_labels))),
            c(16L, 18L, 38L, 7L)),
  # Taking the classes into account widens every standard error.
  all(star_schools$sei > star_schools$sei_ols)
)

# ---- what the standard errors rest on ------------------------------------------
# The contrast of a school and its standard error in closed form, given the
# variance of the class intercepts in small and in regular classes and the
# residual variance: in each kind of class, a mean of the class means with
# weights 1 / (class variance + residual variance / class size).
by_class <- do.call(rbind, lapply(split(complete, complete$class, drop = TRUE), function(d) {
  data.frame(school_id = d$gkschid[1L], small = d$small[1L], n = nrow(d), mean = mean(d$y))
}))
contrasts_given <- function(variant, var_small, var_regular, var_residual) {
  w <- 1 / (ifelse(by_class$small == 1L, var_small, var_regular) + var_residual / by_class$n)
  kind <- function(keep) {
    id <- as.character(by_class$school_id[keep])
    total <- tapply(w[keep], id, sum)
    list(mean = tapply(w[keep] * by_class$mean[keep], id, sum) / total, variance = 1 / total)
  }
  small <- kind(by_class$small == 1L)
  regular <- kind(by_class$small == 0L)
  id <- as.character(star_schools$school_id)
  data.frame(variant = variant, site = star_schools$site,
             yi = as.vector(small$mean[id] - regular$mean[id]),
             sei = sqrt(as.vector(small$variance[id] + regular$variance[id])),
             stringsAsFactors = FALSE)
}
common <- contrasts_given("common", class_sd^2, class_sd^2, residual_sd^2)
stopifnot(max(abs(common$yi - star_schools$yi)) < 1e-10,
          max(abs(common$sei - star_schools$sei)) < 1e-10)

# How well the variance of the classes is known: the 95% interval of the
# ratio of the two standard deviations from the REML profile likelihood,
# with the residual variance set to its estimate at each ratio.
criterion <- lme4::lmer(y ~ 0 + school + school:small + (1 | class), data = complete,
                        REML = TRUE, devFunOnly = TRUE)
state <- environment(criterion)
residual_sd_at <- function(ratio) {
  criterion(ratio)
  sqrt((state$resp$wrss() + state$pp$sqrL(1)) / (nrow(complete) - length(coefficients)))
}
ratio <- lme4::getME(model, "theta")[[1L]]
stopifnot(abs(ratio - class_sd / residual_sd) < 1e-8, abs(residual_sd_at(ratio) - residual_sd) < 1e-8)
cut <- criterion(ratio) + stats::qchisq(0.95, 1)
ratio_limits <- c(
  lower = stats::uniroot(function(r) criterion(r) - cut, c(ratio / 3, ratio), tol = 1e-10)$root,
  upper = stats::uniroot(function(r) criterion(r) - cut, c(ratio, 3 * ratio), tol = 1e-10)$root
)
residual_sd_limits <- vapply(ratio_limits, residual_sd_at, numeric(1))
class_sd_limits <- ratio_limits * residual_sd_limits

# A variance of the class intercepts for each kind of class. Every class is
# of one kind, so the two random terms are the intercepts of the small and of
# the regular classes.
by_kind <- lme4::lmer(y ~ 0 + school + school:small + (0 + small | class) +
                        (0 + regular | class), data = complete, REML = TRUE)
stopifnot(!lme4::isSingular(by_kind), length(by_kind@optinfo$conv$lme4$messages) == 0L)
kinds <- as.data.frame(lme4::VarCorr(by_kind))
kind_sd <- c(small = kinds$sdcor[kinds$var1 %in% "small"],
             regular = kinds$sdcor[kinds$var1 %in% "regular"],
             residual = kinds$sdcor[kinds$grp == "Residual"])
kind_statistic <- as.numeric(2 * (stats::logLik(by_kind, REML = TRUE) -
                                    stats::logLik(model, REML = TRUE)))
variants <- rbind(
  common,
  contrasts_given("lower", class_sd_limits[["lower"]]^2, class_sd_limits[["lower"]]^2,
                  residual_sd_limits[["lower"]]^2),
  contrasts_given("upper", class_sd_limits[["upper"]]^2, class_sd_limits[["upper"]]^2,
                  residual_sd_limits[["upper"]]^2),
  contrasts_given("by_kind", kind_sd[["small"]]^2, kind_sd[["regular"]]^2, kind_sd[["residual"]]^2)
)
# The closed form gives what the model with two variances gives.
kind_coefficients <- lme4::fixef(by_kind)
kind_contrast <- grep(":small$", names(kind_coefficients))
kind_at <- match(as.character(star_schools$school_id),
                 sub("^school", "", sub(":small$", "", names(kind_coefficients)[kind_contrast])))
stopifnot(
  !anyNA(variants), all(variants$sei > 0),
  max(abs(variants$yi[variants$variant == "by_kind"] -
            kind_coefficients[kind_contrast][kind_at])) < 1e-10,
  max(abs(variants$sei[variants$variant == "by_kind"] -
            sqrt(diag(as.matrix(stats::vcov(by_kind))))[kind_contrast][kind_at])) < 1e-10
)

# The variance of the effect across schools estimated in one step from the
# students, with the same class intercepts.
one_step <- lme4::lmer(y ~ 0 + school + small + (0 + small | school) + (1 | class),
                       data = complete, REML = TRUE)
no_variation <- lme4::lmer(y ~ 0 + school + small + (1 | class), data = complete, REML = TRUE)
stopifnot(!lme4::isSingular(one_step))
one_step_terms <- as.data.frame(lme4::VarCorr(one_step))
first_stage <- list(
  variants = variants,
  class_sd = class_sd, residual_sd = residual_sd,
  class_sd_limits = class_sd_limits, residual_sd_limits = residual_sd_limits,
  # Classes less the cells of school by kind of class: the comparisons that
  # the variance of the classes is estimated from.
  between_class_df = nlevels(complete$class) - nrow(unique(complete[c("school", "small")])),
  by_kind = list(sd = kind_sd, statistic = kind_statistic, df = 1L,
                 p = stats::pchisq(kind_statistic, 1, lower.tail = FALSE)),
  one_step = list(
    mean = unname(lme4::fixef(one_step)[["small"]]),
    mean_se = sqrt(diag(as.matrix(stats::vcov(one_step))))[["small"]],
    sd = one_step_terms$sdcor[one_step_terms$grp == "school"],
    statistic = as.numeric(2 * (stats::logLik(one_step, REML = TRUE) -
                                  stats::logLik(no_variation, REML = TRUE)))
  ),
  lme4_version = as.character(utils::packageVersion("lme4")),
  r_version = paste(R.version$major, R.version$minor, sep = ".")
)

attr(star_schools, "source") <- list(
  study = "Tennessee's Student Teacher Achievement Ratio (STAR) project, kindergarten",
  data = "Achilles et al. (2008), Harvard Dataverse, doi:10.7910/DVN/SIWH9F, version 1.0",
  license = "CC0 1.0",
  source_url = "https://doi.org/10.7910/DVN/SIWH9F",
  raw_filename = "STAR_Students.tab",
  raw_bytes = file_hashes$bytes,
  raw_md5 = file_hashes$md5,
  raw_sha256 = file_hashes$sha256,
  raw_body_md5 = file_hashes$body_md5,
  sample = "6,325 students in kindergarten; 5,786 with a reading and a mathematics score; 325 classes; 79 schools",
  outcome = "mean of the standardized reading and mathematics scale scores (Stanford Achievement Test), standardized again",
  contrast = "small classes minus regular classes, with or without a teacher's aide",
  estimator = "lme4::lmer(y ~ 0 + school + school:small + (1 | class), REML = TRUE)",
  class_sd = class_sd, residual_sd = residual_sd,
  ols = "yi_ols and sei_ols: lm(y ~ small) in each school, students taken as independent",
  lme4_version = first_stage$lme4_version, r_version = first_stage$r_version
)

# ---- save ----------------------------------------------------------------------
data_file <- file.path("data", "star_schools.rda")
first_stage_file <- file.path("vignettes", "articles", "class-size-first-stage.rds")
write_files <- TRUE
if (file.exists(data_file)) {
  stored <- new.env()
  load(data_file, envir = stored)
  estimates <- c("yi", "sei", "vi", "yi_ols", "sei_ols")
  others <- setdiff(names(star_schools), estimates)
  if (!identical(names(stored$star_schools), names(star_schools)) ||
      !identical(unclass(stored$star_schools)[others], unclass(star_schools)[others])) {
    stop(sprintf("Columns other than the estimates differ from those in %s. Remove the file to replace it.",
                 data_file))
  }
  moved <- max(abs(as.matrix(stored$star_schools[estimates]) - as.matrix(star_schools[estimates])))
  if (moved > 1e-6) {
    stop(sprintf("The estimates differ from those in %s by up to %.3g. Remove the file to replace them.",
                 data_file, moved))
  }
  if (moved > 0) {
    write_files <- FALSE
    message(sprintf(paste0("The estimates differ from those in %s by up to %.3g. That file is left ",
                           "as it is, and %s is not written."), data_file, moved, first_stage_file))
  } else {
    message(sprintf("The estimates equal those in %s; the file is written again.", data_file))
  }
}
if (write_files) {
  dir.create("data", showWarnings = FALSE)
  save(star_schools, file = data_file, compress = "xz", version = 3)
  if (dir.exists(dirname(first_stage_file))) {
    saveRDS(first_stage, first_stage_file, compress = "xz")
  } else {
    message(sprintf("%s does not exist; %s is not written.", dirname(first_stage_file),
                    basename(first_stage_file)))
  }
}

if (length(args) >= 2L && nzchar(args[2L])) {
  dir.create(args[2L], recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(star_schools, file.path(args[2L], "star-schools.csv"), row.names = FALSE)
  utils::write.csv(
    data.frame(students = nrow(k), complete = nrow(complete), classes = nrow(classes),
               schools = nrow(star_schools), class_sd = class_sd,
               residual_sd = residual_sd, mean_yi = mean(star_schools$yi),
               mean_sei = mean(star_schools$sei), mean_sei_ols = mean(star_schools$sei_ols),
               raw_md5 = file_hashes$md5, raw_body_md5 = file_hashes$body_md5),
    file.path(args[2L], "star-sample.csv"), row.names = FALSE
  )
}
cat(sprintf(paste0("%d schools; %d students; %d classes; mean estimate %.4f; ",
                   "standard errors %.3f to %.3f; class SD %.4f; residual SD %.4f\n"),
            nrow(star_schools), sum(star_schools$n_students), sum(star_schools$n_classes),
            mean(star_schools$yi), min(star_schools$sei), max(star_schools$sei),
            class_sd, residual_sd))
cat(sprintf("class SD %.3f to %.3f (95%%); by kind of class %.3f small, %.3f regular (p = %.3f)\n",
            class_sd_limits[["lower"]], class_sd_limits[["upper"]], kind_sd[["small"]],
            kind_sd[["regular"]], first_stage$by_kind$p))
cat(sprintf("%s %d bytes\n", data_file, file.size(data_file)))
