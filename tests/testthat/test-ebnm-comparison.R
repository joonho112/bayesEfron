# The stored comparison with empirical Bayes methods: its tables hold the
# methods and measures that the documentation reads, the rows of this
# package are those of the stored fits, and the criterion for the large
# panels is met. The file is made by data-raw/ebnm-comparison.R.

ebnm_stored <- function() {
  readRDS(system.file("examples", "ebnm_comparison.rds", package = "bayesEfron"))
}
ebnm_methods <- c("No shrinkage", "Normal", "Scale mixture of normals", "ash", "Unimodal", "NPMLE")

test_that("the file has its parts and says how it was made", {
  stored <- ebnm_stored()
  expect_named(stored, c("krw", "twin", "equal_se", "panels", "settings", "about"))
  expect_named(stored$krw, c("table", "firms", "cdf"))
  expect_named(stored$twin, c("table", "cdf"))
  expect_named(stored$equal_se, c("table", "g", "agreement", "repeated"))
  expect_named(stored$panels, c("table", "by_panel", "convergence"))
  expect_true(all(c("bayesEfron", "ebnm", "deconvolveR", "cmdstan") %in%
                    names(stored$about$versions)))
  expect_identical(stored$about$interval, "equal-tailed 90%")
  # Every method of the tables has its call.
  used <- unique(c(stored$krw$table$method, stored$twin$table$method,
                   stored$equal_se$table$method, stored$panels$table$method))
  expect_true(all(used %in% stored$settings$method))
})

test_that("the rows of the callback gaps for this package are those of the stored fit", {
  table <- ebnm_stored()$krw$table
  expect_identical(table$method, c(ebnm_methods, "g of bayesEfron, plugged in", "bayesEfron"))
  sites <- site_summary(bayesEfron::krw_fit)
  own <- table[table$method == "bayesEfron", ]
  expect_equal(own$mean_sd, mean(sites$sd))
  expect_equal(own$mean_length, mean(sites$upper - sites$lower))
  expect_identical(c(own$correlation, own$sd_ratio, own$length_ratio), c(1, 1, 1))
  expect_identical(own$largest_difference_in_se, 0)
  # Taking g as known gives nearly the same means and intervals that are
  # not longer.
  plugged <- table[table$method == "g of bayesEfron, plugged in", ]
  expect_gt(plugged$correlation, 0.99)
  expect_lte(plugged$length_ratio, 1)
  expect_gt(plugged$length_ratio, 0.95)
  expect_true(all(table$correlation > 0.9))
  # The nonparametric maximum is the largest of the log likelihoods.
  expect_identical(table$method[which.max(table$loglik)], "NPMLE")
  expect_true(is.na(table$loglik[table$method == "No shrinkage"]))
})

test_that("the firms and the distribution functions of the callback gaps are complete", {
  stored <- ebnm_stored()$krw
  firms <- stored$firms
  expect_identical(nrow(firms), 80L)
  expect_identical(as.integer(table(firms$method)), rep(10L, 8L))
  expect_true(all(firms$firm %in% bayesEfron::krw_discrimination$firm_name))
  own <- firms[firms$method == "bayesEfron", ]
  sites <- site_summary(bayesEfron::krw_fit)
  expect_identical(own$firm, sites$label[order(-sites$mean)[1:10]])
  expect_equal(own$mean, sort(sites$mean, decreasing = TRUE)[1:10])

  cdf <- stored$cdf
  expect_setequal(unique(cdf$method), c(setdiff(ebnm_methods, "No shrinkage"), "bayesEfron"))
  for (method in unique(cdf$method)) {
    values <- cdf$cdf[cdf$method == method]
    expect_true(all(diff(values) >= -1e-12), info = method)
    expect_true(all(values >= 0 & values <= 1 + 1e-12), info = method)
  }
  band <- cdf[cdf$method == "bayesEfron", ]
  expect_true(all(band$lower <= band$cdf + 1e-12 & band$cdf <= band$upper + 1e-12))
  # The curves can be laid over the plot of g.
  others <- lapply(split(cdf[cdf$method != "bayesEfron", c("x", "cdf")],
                         cdf$method[cdf$method != "bayesEfron"]), as.data.frame)
  plot <- plot_g(bayesEfron::krw_fit, type = "cdf", compare = others[c("NPMLE", "Normal")])
  expect_s3_class(plot, "ggplot")
})

test_that("the rows of the simulated sites for this package are those of the stored fit", {
  table <- ebnm_stored()$twin$table
  expect_identical(table$method, c(ebnm_methods, "Known g", "g of bayesEfron, plugged in",
                                   "bayesEfron", "bayesEfron, range from the true effects"))
  truth <- bayesEfron::twin_towers$theta_true
  sites <- site_summary(bayesEfron::twin_towers_fit)
  own <- table[table$method == "bayesEfron", ]
  inside <- truth >= sites$lower & truth <= sites$upper
  expect_equal(own$rmse, sqrt(mean((sites$mean - truth)^2)))
  expect_equal(own$inclusion, mean(inside))
  expect_equal(own$mean_length, mean(sites$upper - sites$lower))
  # Every method has every measure, apart from what a flat prior cannot have.
  measures <- c("rmse", "rmse_precise", "rmse_middle", "rmse_imprecise", "inclusion",
                "inclusion_se", "mean_length")
  expect_false(anyNA(table[measures]))
  expect_false(anyNA(table[table$method != "No shrinkage", c("w1", "loglik")]))
  expect_equal(table$inclusion_se, sqrt(table$inclusion * (1 - table$inclusion) / 1500))
  expect_true(all(table$rmse_precise < table$rmse_imprecise))
  # The true distribution is at distance zero from itself, and no estimate
  # of g is closer to it.
  expect_lt(table$w1[table$method == "Known g"], 1e-8)
  expect_true(all(table$w1[!table$method %in% c("Known g", "No shrinkage")] > 0.01))
})

test_that("the comparison with equal standard errors has its three rows and its repetitions", {
  stored <- ebnm_stored()$equal_se
  expect_identical(stored$table$method,
                   c("Efron's estimator", "g of bayesEfron, plugged in", "bayesEfron"))
  expect_false(anyNA(stored$table[setdiff(names(stored$table), "seconds")]))
  expect_identical(nrow(stored$g), as.integer(stored$agreement$L))
  expect_equal(sum(stored$g$bayesEfron) * diff(stored$g$x[1:2]), 1, tolerance = 1e-6)
  expect_equal(sum(stored$g$efron) * diff(stored$g$x[1:2]), 1, tolerance = 1e-6)
  repeated <- stored$repeated
  expect_identical(nrow(repeated), 6L)
  expect_identical(repeated$M, rep(c(6L, 8L), each = 3L))
  expect_identical(unique(repeated$panels), 20L)
  expect_true(all(repeated$inclusion_se > 0 & repeated$inclusion_se < 0.02))
})

test_that("the panels have every method with its Monte Carlo error", {
  stored <- ebnm_stored()$panels
  table <- stored$table
  expect_identical(table$K, rep(c(100L, 1500L), each = 4L))
  expect_identical(table$method, rep(c("bayesEfron", "g of bayesEfron, plugged in", "Normal",
                                       "NPMLE"), 2L))
  expect_identical(table$panels, rep(c(100L, 20L), each = 4L))
  errors <- c("rmse_se", "inclusion_se", "mean_length_se")
  expect_false(anyNA(table))
  expect_true(all(unlist(table[errors]) > 0))
  # The table is the mean of the values of the panels.
  by_panel <- stored$by_panel
  expect_identical(nrow(by_panel), 4L * (100L + 20L))
  one <- by_panel[by_panel$K == 100L & by_panel$method == "NPMLE", ]
  row <- table[table$K == 100L & table$method == "NPMLE", ]
  expect_equal(row$inclusion, mean(one$inclusion))
  expect_equal(row$inclusion_se, stats::sd(one$inclusion) / sqrt(100))
})

test_that("no fit of the large panels with the defaults has an R-hat above 1.05", {
  convergence <- ebnm_stored()$panels$convergence
  expect_identical(convergence$setting, c("defaults", "M = 8", "settings = \"0.3.0\""))
  expect_identical(convergence$panels, rep(20L, 3L))
  defaults <- convergence[convergence$setting == "defaults", ]
  expect_identical(defaults$rhat_above_1.05, 0L)
  expect_lt(defaults$rhat_max, 1.05)
  expect_equal(defaults$divergences, 0)
})

test_that("the tables of the panels in inst/verification are those of the stored file", {
  stored <- ebnm_stored()$panels
  read <- function(name) {
    utils::read.csv(system.file("verification", name, package = "bayesEfron"),
                    stringsAsFactors = FALSE, check.names = FALSE)
  }
  convergence <- read("panels-convergence.csv")
  methods <- read("panels-methods.csv")
  shared <- names(stored$convergence)
  expect_equal(convergence[shared], stored$convergence, ignore_attr = TRUE)
  expect_equal(methods[names(stored$table)], stored$table, ignore_attr = TRUE)
  expect_true(all(c("generated_on", "cmdstan_version") %in% names(convergence)))
})

test_that("the stored examples stay within their size", {
  files <- list.files(system.file("examples", package = "bayesEfron"), full.names = TRUE)
  expect_lte(sum(file.size(files)), 0.7 * 2^20)
})
