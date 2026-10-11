# The table of site effects and the labels of the sites.

labelled_fit <- function(labels) {
  fit <- readRDS(test_path("_fixtures", "short_fit.rds"))
  fit$metadata$theta_summary$label <- labels
  validate_bef_fit_re(fit)
}

test_that("the table has the columns of this version and repeats three of them under the names of version 0.3.0", {
  fit <- readRDS(test_path("_fixtures", "short_fit.rds"))
  table <- as.data.frame(fit)
  expect_named(table, .bef_site_table_columns())
  expect_identical(table$label, as.character(1:5))
  expect_identical(table$estimate, fit$metadata$data_list$theta_hat)
  expect_identical(table$se, fit$metadata$data_list$sigma)
  expect_identical(table$hpdi_lower, table$lower)
  expect_identical(table$hpdi_upper, table$upper)
  expect_identical(table$map, table$mode)
  expect_identical(attr(table, "level"), 0.9)
  # The table is what the engine computes from log_g.
  engine <- .bef_engine(fit)
  expect_equal(table$mean, .bef_site_moments(engine)$mean)
  expect_equal(table$sd, sqrt(.bef_site_moments(engine)$var))
  expect_equal(unname(as.matrix(table[, c("lower", "upper")])),
               .bef_site_quantile(engine, c(0.05, 0.95)))
  expect_true(all(table$mode %in% engine$grid))
  expect_equal(unname(coef(fit, type = "mode")), table$mode)
})

test_that("labels name what coef(), vcov() and confint() return", {
  labels <- c("north", "south", "east", "west", "centre")
  fit <- labelled_fit(labels)
  expect_named(coef(fit), labels)
  expect_identical(dimnames(vcov(fit)), list(labels, labels))
  intervals <- confint(fit)
  expect_identical(rownames(intervals), labels)
  expect_identical(intervals$site, 1:5)
  expect_identical(as.data.frame(fit)$label, labels)
  expect_true(any(grepl("centre", format(summary(fit)), fixed = TRUE)))
  # Without labels the rows are numbered and the names are the site numbers.
  plain <- readRDS(test_path("_fixtures", "short_fit.rds"))
  expect_named(coef(plain), as.character(1:5))
  expect_identical(rownames(confint(plain)), as.character(1:5))
})

test_that("numbers in parm are positions and character strings are labels", {
  # Labels that are numbers in another order: "5" is the first site.
  fit <- labelled_fit(c("5", "4", "3", "2", "1"))
  by_position <- confint(fit, parm = 5)
  by_label <- confint(fit, parm = "5")
  expect_identical(by_position$site, 5L)
  expect_identical(by_label$site, 1L)
  expect_identical(rownames(by_label), "5")
  expect_equal(by_label$mean, unname(coef(fit)[1L]))
  expect_identical(confint(fit, parm = c("2", "5"))$site, c(4L, 1L))
  expect_error(confint(fit, parm = "6"), "does not have", class = "bef_invalid_args")
  expect_error(confint(fit, parm = 6), "between 1 and 5", class = "bef_invalid_args")
  expect_error(confint(fit, parm = TRUE), class = "bef_invalid_args")
})

test_that("as_bef_data() takes labels from names, a column and slab, and refuses repeats", {
  y <- c(a = 0.1, b = 0.2, c = -0.1, d = 0.4, e = 0.3)
  s <- rep(0.2, 5L)
  expect_identical(as_bef_data(list(theta_hat = y, sigma = s))$names, names(y))
  df <- data.frame(yi = unname(y), sei = s, school = c("P", "Q", "R", "S", "T"))
  data <- as_bef_data(df, site = "school")
  expect_identical(data$names, df$school)
  expect_identical(
    as.data.frame(data),
    data.frame(site = 1:5, label = df$school, estimate = unname(y), se = s)
  )
  expect_identical(as.data.frame(as_bef_data(df))$label, as.character(1:5))

  df$school[2L] <- "P"
  expect_error(as_bef_data(df, site = "school"), "must be distinct", class = "bef_invalid_args")
  expect_error(as_bef_data(list(theta_hat = c(a = 1, a = 2, b = 3, c = 4, d = 5), sigma = s)),
               "must be distinct", class = "bef_invalid_args")

  # metafor keeps study labels in the attribute `slab` of the estimates.
  es <- structure(
    data.frame(yi = unname(y), vi = s^2),
    class = c("escalc", "data.frame")
  )
  attr(es$yi, "slab") <- c("Smith", "Jones", "Lee", "Kim", "Park")
  expect_identical(as_bef_data(es)$names, c("Smith", "Jones", "Lee", "Kim", "Park"))
})

test_that("a long table shows its first and last sites with the omitted ones between", {
  lines <- format(summary(bayesEfron::krw_fit))
  from <- grep("Site effects (posterior summaries, 90% intervals):", lines, fixed = TRUE)
  expect_length(from, 1L)
  table <- lines[(from + 1L):length(lines)]
  expect_length(table, 14L)
  expect_match(table[1L], "site +mean +sd +lower +upper +mode")
  expect_match(table[8L], "... 85 sites omitted ...", fixed = TRUE)
  # The stored fit calls its sites by the names of the firms.
  expect_match(table[2L], "^  Rite Aid ")
  expect_match(table[14L], "^  UGI ")
  expect_identical(site_summary(bayesEfron::krw_fit)$label,
                   bayesEfron::krw_discrimination$firm_name)
  expect_true(any(grepl("each is a posterior mean", lines, fixed = TRUE)))
})

test_that("the note about the renamed columns is given once, in an interactive session only", {
  notes <- get(".bef_notes", asNamespace("bayesEfron"))
  rm(list = ls(notes), envir = notes)
  expect_silent(.bef_note_site_table_names(interactive_fun = function() FALSE))
  expect_message(.bef_note_site_table_names(interactive_fun = function() TRUE),
                 "`lower`, `upper` and `mode`")
  expect_silent(.bef_note_site_table_names(interactive_fun = function() TRUE))
  rm(list = ls(notes), envir = notes)
})

test_that("the fit and the grid take two vectors, a data frame or a bef_data object", {
  df <- data.frame(yi = c(-0.21, 0.04, 0.19, 0.38, 0.61, 0.12),
                   sei = c(0.18, 0.15, 0.22, 0.19, 0.24, 0.2),
                   school = c("P", "Q", "R", "S", "T", "U"))
  context <- function(...) {
    .bef_fit_prepare_context(quote(f()), ..., check_installed = FALSE)
  }
  vectors <- context(theta_hat = df$yi, sigma = df$sei)
  from_df <- context(theta_hat = df)
  from_data <- context(theta_hat = as_bef_data(df, site = "school"))
  expect_identical(from_df$stan_data, vectors$stan_data)
  expect_identical(from_data$stan_data, vectors$stan_data)
  expect_null(vectors$bef_data$names)
  expect_identical(from_data$bef_data$names, df$school)

  expect_identical(make_efron_grid(df, grid_method = "paper_realdata"), make_efron_grid(df$yi, df$sei, grid_method = "paper_realdata"))
  expect_identical(make_efron_grid(as_bef_data(df), grid_method = "paper_realdata"), make_efron_grid(df$yi, df$sei, grid_method = "paper_realdata"))
  expect_error(context(theta_hat = df, sigma = df$sei), "must be left out",
               class = "bef_invalid_args")
  expect_error(make_efron_grid(df, df$sei, grid_method = "paper_realdata"), "must be left out", class = "bef_invalid_args")
  # The smallest number of sites is the same everywhere.
  expect_error(context(theta_hat = df$yi[1:4], sigma = df$sei[1:4]), class = "bef_invalid_args")
  expect_error(as_bef_data(df[1:4, ]), "at least 5", class = "bef_invalid_args")
})

test_that("a variance is not taken for a standard error without notice", {
  df <- data.frame(yi = c(0.1, 0.2, -0.1, 0.4, 0.3), vi = c(0.04, 0.01, 0.09, 0.04, 0.0225))
  # Estimates and variances only: the square roots are used, and said so.
  expect_message(data <- as_bef_data(df), "square roots of the variances in `vi`")
  expect_equal(data$sigma, sqrt(df$vi))
  expect_silent(explicit <- as_bef_data(df, sigma = "vi", variance = TRUE))
  expect_equal(explicit$sigma, sqrt(df$vi))
  # A column named like a variance given as the standard errors.
  expect_warning(wrong <- as_bef_data(df, sigma = "vi"), "suggests variances")
  expect_equal(wrong$sigma, df$vi)
  # With a column of standard errors present nothing is said.
  df$sei <- sqrt(df$vi)
  expect_silent(as_bef_data(df))
})
