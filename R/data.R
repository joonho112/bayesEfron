#' Teacher expectancy studies
#'
#' Results of 19 experiments on whether telling teachers that some of their
#' pupils are likely to show intellectual growth raises those pupils' IQ
#' scores. Raudenbush (1984) related the size of the effect to the time
#' that teachers had known their pupils before the expectation was induced,
#' and Raudenbush and Bryk (1985) used the studies to illustrate empirical
#' Bayes meta-analysis. The studies differ widely in size, and so do their
#' standard errors.
#'
#' @format A data frame with 19 rows and 11 columns:
#' \describe{
#'   \item{study}{Study number.}
#'   \item{author, year}{Authors and year of publication.}
#'   \item{weeks}{Weeks of contact between the teacher and the pupils
#'     before the expectation was induced.}
#'   \item{setting}{Whether the test was given to a group (`"group"`) or
#'     individually (`"indiv"`).}
#'   \item{tester}{Whether the person who gave the test knew which pupils
#'     had been named to the teacher (`"aware"`) or not (`"blind"`).}
#'   \item{n1i, n2i}{Numbers of pupils in the experimental and the control
#'     group.}
#'   \item{yi}{Standardized mean difference in IQ score. Positive values
#'     mean that the pupils named to the teacher scored higher on average.}
#'   \item{vi}{Sampling variance of `yi`.}
#'   \item{sei}{Standard error of `yi`, the square root of `vi`.}
#' }
#'
#' @source Raudenbush, S. W. and Bryk, A. S. (1985). Empirical Bayes
#'   meta-analysis. *Journal of Educational Statistics*, 10(2), 75-98.
#'   \doi{10.3102/10769986010002075}. The columns `setting` and `tester`
#'   are from Raudenbush, S. W. (1984). Magnitude of teacher expectancy
#'   effects on pupil IQ as a function of the credibility of expectancy
#'   induction: A synthesis of findings from 18 experiments. *Journal of
#'   Educational Psychology*, 76(1), 85-97.
#'   \doi{10.1037/0022-0663.76.1.85}. The values were taken from
#'   `dat.raudenbush1985` in the metadat package, which distributes the
#'   same data.
#'
#' @seealso [raudenbush_fit], the model fitted to these data.
#'
#' @examples
#' head(raudenbush1985)
#' as_bef_data(raudenbush1985)
"raudenbush1985"

#' Model fitted to the teacher expectancy studies
#'
#' The object returned by [bayes_efron_fit()] for the [raudenbush1985] data
#' with the default settings. The examples and the vignettes use it, so
#' that they run without CmdStan.
#'
#' The model was fitted with
#' `bayes_efron_fit(raudenbush1985$yi, raudenbush1985$sei, seed = 1985)`,
#' that is, with four chains of 1,000 warmup and 3,000 sampling iterations
#' each. To keep the package small, one draw in fifteen is stored, so the
#' object holds 800 of the 12,000 draws, and its posterior summaries are
#' computed from those 800. The convergence diagnostics that [diagnose()]
#' returns are those of the full run. One of the 12,000 iterations reached
#' the maximum tree depth, which [diagnose()] reports as a warning;
#' `vignette("diagnostics")` explains how to read it.
#'
#' @format An object of class `bef_fit_re`; see [bef_fit] for its
#'   components.
#'
#' @source Fitted with CmdStan 2.38.0 on 2026-10-02. The script is
#'   `data-raw/raudenbush.R` in the source repository.
#'
#' @seealso [raudenbush1985], [bayes_efron_fit()]
#'
#' @examples
#' raudenbush_fit
#' summary(raudenbush_fit)
"raudenbush_fit"
