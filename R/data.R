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
#' @seealso [raudenbush_fit], `vignette("teacher-expectancy")`
#'
#' @examples
#' head(raudenbush1985)
#' as_bef_data(raudenbush1985)
"raudenbush1985"

#' Model fitted to the teacher expectancy studies
#'
#' @description
#' A fitted model for the 19 studies in [raudenbush1985], as returned by
#' [bayes_efron_fit()] with the default settings. It is thinned to 800
#' posterior draws, is used in `vignette("teacher-expectancy")`, and was
#' created with
#'
#' ```
#' bayes_efron_fit(raudenbush1985, seed = 1985, thin = 15)
#' ```
#'
#' @details
#' Each of the 4 chains ran 1,000 warmup and 3,000 sampling iterations,
#' and one draw out of every 15 was stored. Summaries are computed from
#' the 800 stored draws, whereas the sampler diagnostics reported by
#' [diagnose()] are computed from all 12,000. For these 19 studies the
#' `"pilot"` method does not narrow the range, so the grid is the same as
#' with `grid_method = "paper_realdata"`. The object contains the draws of
#' \eqn{g} but no site-level draws, and it can be used without CmdStan.
#'
#' @format An object of class `bef_fit_re`; see [bef_fit] for its
#'   components and methods.
#'
#' @source Fitted with CmdStan 2.38.0 on 2026-10-10. The script is
#'   `data-raw/example-fits.R` in the source repository.
#'
#' @seealso [raudenbush1985], [bayes_efron_fit()],
#'   `vignette("teacher-expectancy")`
#'
#' @examples
#' raudenbush_fit
#' summary(raudenbush_fit)
"raudenbush_fit"

#' Racial callback gaps among 97 large firms
#'
#' Firm-specific differences in callback rates for applications with
#' distinctively white and distinctively Black names. The original
#' correspondence experiment is from Kline, Rose and Walters (2022).
#' Following their 2024 analysis and the replication code of Walters
#' (2024), firms with at least 40 sampled vacancies and an overall callback
#' rate of at least 3 percent are included: 97 firms, 78,910 applications
#' and 10,453 vacancies. Each firm is one site of the model; the standard
#' errors are clustered by vacancy within firm.
#'
#' @format A data frame with 97 rows and 8 columns:
#' \describe{
#'   \item{site}{Sequential site number, in ascending order of `firm_id`.}
#'   \item{firm_id}{Firm identifier in the original data.}
#'   \item{firm_name}{Name of the firm as supplied in the replication data.}
#'   \item{yi}{Difference in callback probability, white minus Black,
#'     estimated by ordinary least squares. Positive values favor
#'     applications with distinctively white names. The values are
#'     differences of proportions: multiply by 100 for percentage points.}
#'   \item{sei}{Standard error of `yi`, clustered by vacancy, with the CR1
#'     finite-sample correction used by Stata's clustered OLS regression.}
#'   \item{vi}{Sampling variance, `sei^2`.}
#'   \item{n_jobs}{Number of distinct sampled vacancies in the firm.}
#'   \item{n_obs}{Number of applications in the firm.}
#' }
#'
#' @source Kline, P., Rose, E. K. and Walters, C. R. (2022). Systemic
#'   discrimination among large U.S. employers. *The Quarterly Journal of
#'   Economics*, 137(4), 1963-2036. \doi{10.1093/qje/qjac024}. The analysis
#'   sample follows Kline, P., Rose, E. K. and Walters, C. R. (2024). A
#'   discrimination report card. *American Economic Review*, 114(8),
#'   2472-2525. \doi{10.1257/aer.20230700}. Data and regression code were
#'   obtained from the replication package accompanying Walters, C. R.
#'   (2024). Empirical Bayes methods in labor economics. In *Handbook of
#'   Labor Economics*, 5, 183-260.
#'   \doi{10.1016/bs.heslab.2024.11.001}. The package is available through
#'   \url{https://sites.google.com/view/christopher-walters/research}.
#'   The Bayesian reanalysis is described in Lee and Sui (2025),
#'   \doi{10.3390/math13162639}. `data-raw/krw-discrimination.R` in the
#'   source repository reproduces the estimates from `krw_data.dta`, and
#'   the attribute `source` of the data frame records where the data came
#'   from.
#'
#' @seealso [as_bef_data()], `vignette("bayesEfron")`
#'
#' @examples
#' head(krw_discrimination)
#' as_bef_data(krw_discrimination)
"krw_discrimination"

#' Twin Towers simulated data
#'
#' A simulated data set of 1,500 sites from the public replication data of
#' Lee and Sui (2025). The true effects are the 1,500 values of
#' `disjointTheta` in the deconvolveR package, which lie in two separated
#' groups; normally distributed sampling errors obscure the two groups in
#' the estimates. The largest standard error is nine times the smallest
#' (`R = 9`), and the design parameter `I` is 0.7; both are defined below.
#'
#' The sampling variances are equally spaced on a log scale and were
#' assigned to the sites at random before the normal errors were added. The
#' data set is taken as it stands from the replication data of that study,
#' where it is called a panel.
#'
#' @format A data frame with 1,500 rows and 7 columns:
#' \describe{
#'   \item{site}{Site number, in the order of the replication data.}
#'   \item{theta_true}{True effect, a value of `disjointTheta`.}
#'   \item{yi}{Estimate of the effect.}
#'   \item{sei}{Standard error of the estimate.}
#'   \item{vi}{Sampling variance, `sei^2`.}
#'   \item{I}{Design parameter, \eqn{1 / (1 + \bar v)}{1 / (1 + v_bar)}
#'     with \eqn{\bar v}{v_bar} the geometric mean of `vi`. It would be the
#'     reliability of the estimates if the true effects had variance one;
#'     their variance in this data set is 2.11.}
#'   \item{R}{Ratio of the largest to the smallest standard error. The
#'     ratio of the variances is `R^2`, 81 in this data set.}
#' }
#'
#' @source Lee, J. and Sui, D. (2025). Fully Bayesian inference for
#'   meta-analytic deconvolution using Efron's log-spline prior.
#'   *Mathematics*, 13(16), 2639. \doi{10.3390/math13162639}.
#'   Public replication data:
#'   \url{https://github.com/joonho112/fully-bayes-efron-prior}.
#'   The true effects are from the Twin Towers example distributed in
#'   deconvolveR; see Narasimhan, B. and Efron, B. (2020). deconvolveR:
#'   A G-modeling program for deconvolution and empirical Bayes estimation.
#'   *Journal of Statistical Software*, 94(11), 1-20.
#'   \doi{10.18637/jss.v094.i11}. `data-raw/twin-towers.R` in the source
#'   repository selects the data set and checks it, and the attribute
#'   `source` of the data frame records where it came from.
#'
#' @seealso [as_bef_data()], `vignette("twin-towers")`
#'
#' @examples
#' head(twin_towers)
#' as_bef_data(twin_towers)
"twin_towers"

#' Model fitted to the racial callback gaps
#'
#' @description
#' A fitted model for the 97 firms in [krw_discrimination], as returned by
#' [bayes_efron_fit()] with the default settings and with the firm names
#' as site labels. It is thinned to 800 posterior draws and was created
#' with
#'
#' ```
#' bayes_efron_fit(as_bef_data(krw_discrimination, site = "firm_name"),
#'                 seed = 20250807, thin = 15)
#' ```
#'
#' @details
#' Each of the 4 chains ran 1,000 warmup and 3,000 sampling iterations,
#' and one draw out of every 15 was stored. Summaries are computed from
#' the 800 stored draws, whereas the sampler diagnostics reported by
#' [diagnose()] are computed from all 12,000. For these data the `"pilot"`
#' method does not narrow the range: the grid has 101 points and extends
#' beyond the range of the estimates by half the width of that range on
#' each side. The object contains the draws of \eqn{g} but no site-level
#' draws, and it can be used without CmdStan.
#'
#' Sites are in the same order as the rows of [krw_discrimination].
#' Tables, plots and the names of the [coef()] vector use the firm names
#' as site labels, and a site can be selected by position or by name.
#'
#' @format An object of class `bef_fit_re`; see [bef_fit] for its
#'   components and methods.
#' @source Fitted with CmdStan 2.38.0 on 2026-10-10. The script is
#'   `data-raw/example-fits.R` in the source repository. The sources of the
#'   data are in [krw_discrimination].
#' @seealso [krw_discrimination], `vignette("bayesEfron")`
#' @examples
#' krw_fit
#' confint(krw_fit, type = "g")
#' site_summary(krw_fit, sites = c("KFC", "Rite Aid"), scale = 100)
"krw_fit"

#' Model fitted to the Twin Towers simulated data
#'
#' @description
#' A fitted model for the 1,500 simulated sites in [twin_towers], as
#' returned by [bayes_efron_fit()] with 8 spline functions. The true
#' effects, which are also in the data frame, are not used in the fit.
#' The fit is thinned to 800 posterior draws and was created with
#'
#' ```
#' bayes_efron_fit(twin_towers, M = 8, seed = 20250627, thin = 15)
#' ```
#'
#' @details
#' In this data set, 8 spline functions (`M = 8`) give an estimate of
#' \eqn{g} that is closer to the true effects than the default of 6, and,
#' unlike with 6, no iteration reaches the maximum tree depth. See
#' `vignette("twin-towers")` for a comparison of the two settings,
#' including the lower coverage of the site intervals with 8 functions in
#' simulated data sets. For these data
#' the `"pilot"` method narrows the range of `"paper_realdata"`, and the
#' printed fit shows by how much.
#'
#' Each of the 4 chains ran 1,000 warmup and 3,000 sampling iterations,
#' and one draw out of every 15 was stored. Summaries are computed from
#' the 800 stored draws, whereas the sampler diagnostics reported by
#' [diagnose()] are computed from all 12,000. The object contains the
#' draws of \eqn{g} but no site-level draws, and it can be used without
#' CmdStan. Sites are numbered in the order of the data.
#'
#' A comparison of this fit with fits of the same data on other ranges is
#' stored in the file `twin_towers_comparison.rds` in the `examples`
#' directory of the installed package. The other ranges include the range
#' around the true effects (`grid_method = "paper_simulation"`), which is
#' used in the simulations of Lee and Sui (2025), and the wide range of
#' `"paper_realdata"`.
#'
#' @format An object of class `bef_fit_re`; see [bef_fit] for its
#'   components and methods.
#' @source Fitted with CmdStan 2.38.0 on 2026-10-10. The script is
#'   `data-raw/example-fits.R` in the source repository. The sources of the
#'   data are in [twin_towers].
#' @seealso [twin_towers], `vignette("twin-towers")`
#' @examples
#' twin_towers_fit
#' confint(twin_towers_fit, type = "g")
"twin_towers_fit"

#' Small-class effects in 79 schools of Project STAR
#'
#' School-specific estimates of the effect of small classes on achievement
#' in kindergarten, from Tennessee's Student Teacher Achievement Ratio
#' (STAR) project. Within each of 79 schools, students and teachers were
#' assigned at random to small classes or to regular classes, the regular
#' ones with or without a teacher's aide, so the experiment is a multisite
#' trial with the schools as its sites. The estimate of a school is the
#' difference between its small classes and its regular classes of either
#' kind on a composite of the kindergarten reading and mathematics scores.
#'
#' The outcome is the mean of the standardized reading and mathematics
#' scale scores of the Stanford Achievement Test, standardized again, so
#' the estimates are in standard deviations of that composite. Of the 6,325
#' students in kindergarten, the 5,786 with both scores are used; they sat
#' in 325 classes, 127 of them small.
#'
#' The estimates and their standard errors come from one linear mixed model
#' fitted to all the schools by restricted maximum likelihood. The model has
#' an intercept and a small-class contrast for each school and a random
#' intercept for each class. The variance of the class intercepts is taken
#' to be the same for small and for regular classes and in every school,
#' and so is the residual variance, so the standard error of a school
#' depends only on the numbers and the sizes of its small and its regular
#' classes. A school has from one to four small classes, and that
#' number, more than the number of its students, limits how precise its
#' estimate can be.
#'
#' `yi_ols` and `sei_ols` are the same contrast computed within each school
#' by ordinary least squares, with the students of a school taken as
#' independent. That ignores what the students of a class have in common,
#' and the standard error is smaller for every school. The two columns are
#' there to show how much a fit depends on its standard errors; the article
#' "Small classes in 79 schools" on the package website,
#' \url{https://joonho112.github.io/bayesEfron/}, fits both, and it also
#' varies the assumption about the variance of the classes.
#'
#' @format A data frame with 79 rows and 12 columns:
#' \describe{
#'   \item{site}{Sequential site number, in ascending order of `school_id`.}
#'   \item{school_id}{School identifier in the public-use data.}
#'   \item{location}{Kind of area that the school serves: `"inner city"`,
#'     `"suburban"`, `"rural"` or `"urban"`.}
#'   \item{yi}{Estimate of the effect of small classes: small minus
#'     regular, in standard deviations of the composite score.}
#'   \item{sei}{Standard error of `yi`, from the model with a random
#'     intercept for each class.}
#'   \item{vi}{Sampling variance, `sei^2`.}
#'   \item{yi_ols, sei_ols}{The same contrast and its standard error by
#'     ordinary least squares within the school.}
#'   \item{n_students, n_students_small}{Numbers of students with both
#'     scores: in the school, and in its small classes.}
#'   \item{n_classes, n_classes_small}{Numbers of classes: in the school,
#'     and small ones.}
#' }
#'
#' @source Achilles, C. M., Bain, H. P., Bellott, F., Boyd-Zaharias, J.,
#'   Finn, J., Folger, J., Johnston, J. and Word, E. (2008). Tennessee's
#'   Student Teacher Achievement Ratio (STAR) project. Harvard Dataverse,
#'   V1. \doi{10.7910/DVN/SIWH9F}. The data are released under CC0 1.0. The
#'   experiment is described in Finn, J. D. and Achilles, C. M. (1990).
#'   Answers and questions about class size: A statewide experiment.
#'   *American Educational Research Journal*, 27(3), 557-577.
#'   \doi{10.3102/00028312027003557}. `data-raw/star.R` in the source
#'   repository builds the estimates from the student file
#'   `STAR_Students.tab`, and the attribute `source` of the data frame
#'   records the file, the sample and the model.
#'
#' @seealso [as_bef_data()], [bayes_efron_fit()]
#'
#' @examples
#' head(star_schools)
#' as_bef_data(star_schools)
#'
#' # The contrast and standard errors that take the students of a school as
#' # independent.
#' as_bef_data(star_schools, theta_hat = "yi_ols", sigma = "sei_ols")
"star_schools"
