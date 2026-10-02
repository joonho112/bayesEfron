# bayes_efron_compile(): argument checks and what it hands to the model
# loader. Compilation itself is covered in test-stan-model.R and test-live.R.

test_that("bayes_efron_compile() rejects invalid arguments", {
  testthat::local_mocked_bindings(
    .bef_model = function(...) stop("the model loader should not be reached"),
    .package = "bayesEfron"
  )

  cases <- list(
    model_family = quote(bayes_efron_compile("HE")),
    quiet = quote(bayes_efron_compile(quiet = NA)),
    force_recompile = quote(bayes_efron_compile(force_recompile = c(TRUE, FALSE)))
  )
  for (arg in names(cases)) {
    err <- tryCatch(eval(cases[[arg]]), error = identity)
    expect_s3_class(err, "bef_invalid_args")
    expect_s3_class(err, "bef_error")
    expect_identical(err$arg, arg)
    expect_match(conditionMessage(err), arg, fixed = TRUE)
  }
})

test_that("bayes_efron_compile() passes its arguments on and returns the model invisibly", {
  received <- NULL
  model <- structure(list(), class = "fake_cmdstan_model")
  testthat::local_mocked_bindings(
    .bef_model = function(...) {
      received <<- list(...)
      model
    },
    .package = "bayesEfron"
  )

  expect_invisible(out <- bayes_efron_compile(quiet = FALSE, force_recompile = TRUE))
  expect_identical(out, model)
  expect_identical(
    received,
    list(model_family = "RE", force_recompile = TRUE, quiet = FALSE)
  )

  bayes_efron_compile()
  expect_identical(
    received,
    list(model_family = "RE", force_recompile = FALSE, quiet = TRUE)
  )
})
