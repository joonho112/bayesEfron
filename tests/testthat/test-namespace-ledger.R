test_that("EXPORTS_LEDGER matches generated NAMESPACE directives", {
  root <- bef_test_source_root(required = c("NAMESPACE", file.path("inst", "EXPORTS_LEDGER.md")))
  namespace <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
  namespace <- namespace[grepl("^(export|S3method)\\(", namespace)]
  ledger <- readLines(file.path(root, "inst", "EXPORTS_LEDGER.md"), warn = FALSE)

  expect_true(all(vapply(
    namespace,
    function(directive) any(grepl(directive, ledger, fixed = TRUE)),
    logical(1)
  )))

  ledger_directives <- regmatches(
    ledger,
    gregexpr("`((export|S3method)[^`]+)`", ledger)
  )
  ledger_directives <- gsub("`", "", unlist(ledger_directives, use.names = FALSE))
  ledger_directives <- ledger_directives[grepl("^(export|S3method)\\(", ledger_directives)]
  expect_setequal(ledger_directives, namespace)
  expect_false(anyDuplicated(namespace) > 0L)
  expect_false(anyDuplicated(ledger_directives) > 0L)

  exported <- namespace[grepl("^export\\(", namespace)]
  s3_methods <- namespace[grepl("^S3method\\(", namespace)]

  expect_length(exported, 6L)
  expect_length(s3_methods, 25L)
  expect_false(any(grepl("^export\\((new_|validate_)", exported)))

  input_converters <- s3_methods[grepl("^S3method\\(as_bef_data,", s3_methods)]
  diagnostic_producers <- s3_methods[grepl("^S3method\\(diagnose,", s3_methods)]
  bef_fit_methods <- c(
    "S3method(format,bef_fit)",
    "S3method(format,summary.bef_fit)",
    "S3method(logLik,bef_fit)",
    "S3method(nobs,bef_fit)",
    "S3method(posterior::as_draws,bef_fit)",
    "S3method(print,bef_fit)",
    "S3method(print,summary.bef_fit)",
    "S3method(summary,bef_fit)"
  )
  bef_fit_re_methods <- c(
    "S3method(as.data.frame,bef_fit_re)",
    "S3method(coef,bef_fit_re)",
    "S3method(confint,bef_fit_re)",
    "S3method(plot,bef_fit_re)",
    "S3method(summary,bef_fit_re)",
    "S3method(vcov,bef_fit_re)"
  )
  bef_data_methods <- c(
    "S3method(format,bef_data)",
    "S3method(print,bef_data)",
    "S3method(summary,bef_data)"
  )
  bef_diagnostic_methods <- c(
    "S3method(format,bef_diagnostic)",
    "S3method(print,bef_diagnostic)",
    "S3method(summary,bef_diagnostic)"
  )

  expect_length(input_converters, 3L)
  expect_length(diagnostic_producers, 2L)
  expect_setequal(intersect(s3_methods, bef_fit_methods), bef_fit_methods)
  expect_setequal(intersect(s3_methods, bef_fit_re_methods), bef_fit_re_methods)
  expect_setequal(intersect(s3_methods, bef_data_methods), bef_data_methods)
  expect_setequal(intersect(s3_methods, bef_diagnostic_methods), bef_diagnostic_methods)

  ledger_text <- paste(ledger, collapse = "\n")
  expect_match(ledger_text, "\\| Active exported functions \\|\\s*6\\s*\\|")
  expect_match(ledger_text, "\\| Input-conversion S3 registrations \\|\\s*3\\s*\\|")
  expect_match(ledger_text, "\\| Structural method registrations \\|\\s*20\\s*\\|")
  expect_match(ledger_text, "\\| Diagnostic-producer S3 registrations \\|\\s*2\\s*\\|")
  expect_match(ledger_text, "\\| Total S3 registrations \\|\\s*25\\s*\\|")

  internal_names <- c(
    "new_bef_data",
    "new_bef_fit",
    "new_bef_fit_re",
    "new_bef_diagnostic",
    "validate_bef_data",
    "validate_bef_fit",
    "validate_bef_fit_re",
    "validate_bef_diagnostic"
  )
  exported_names <- sub("^export\\(([^\\)]+)\\)$", "\\1", exported)
  expect_false(any(internal_names %in% exported_names))
  expect_true(all(vapply(
    internal_names,
    exists,
    logical(1),
    envir = asNamespace("bayesEfron"),
    inherits = FALSE
  )))
})
