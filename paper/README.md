# Reproduction material for the SoftwareX article on bayesEfron 0.4.0

This folder holds the scripts, data and recorded outputs behind the article
"bayesEfron: An R package for fully Bayesian estimation of the distribution
of effects and of site-specific effects from estimates and standard errors"
(SoftwareX, submitted in 2026). The LaTeX sources of the article and of its
supplementary material go to the journal as supplementary material and are
not kept here. Everything in this folder is excluded from the R package
build (see `.Rbuildignore`).

## Contents

| Folder or file | Content |
|:--|:--|
| `scripts/install-package.sh` | Installs bayesEfron 0.4.0 from a checkout of the tagged commit into a library inside this folder |
| `scripts/analysis.R` | The fits of the two examples (97 firms, Twin Towers) and of the 19 teacher expectancy studies, the comparisons with REML and with ebnm, the numbers registry and the console outputs |
| `scripts/stored-examples.R` | Numbers quoted from the package's stored examples (`inst/examples/`) |
| `scripts/tail-coverage.R` | Coverage of the 90% intervals for subsets of the sites in the two-group simulation designs, from the per-site results in `data-raw/benchmark/results-v040/` |
| `scripts/figures.R`, `figures/fig1_architecture.tex` | The figures of the article and of the supplement; Figure 1 is TikZ |
| `scripts/tables.R` | The tables whose numbers come from `data/` |
| `scripts/check-numbers.R` | Compares every number quoted in the article's sources with the registries; it needs those sources (supplementary material) and is kept here as the record of what was checked |
| `data/numbers.csv` | The 134 numbers quoted in the article, with the value shown |
| `data/stored-examples-numbers.csv`, `data/tail-coverage-numbers.csv` | The two smaller registries |
| `data/tail-coverage.csv` | Coverage by subset of the sites, design and method |
| `data/krw.rds`, `data/twin.rds`, `data/grid-comparison.rds`, `data/evidence.rds` | The summaries, curves and tables behind the figures |
| `data/console/` | The console sessions printed in the supplement |
| `figures/`, `osm/figures/` | The figures as PDF. Working names differ from the article's numbering: `fig1_architecture` is Figure 1, `fig4_evidence` Figure 2, `fig2_firms` Figure 3, `fig3_twin_towers` Figure 4, `fig5_grid_and_uncertainty` Figure 5; `osm/figures/` holds the two figures of the supplement |
| `tables/` | LaTeX tables built from `data/` (Table 3 of the article and the tables of the supplement) |
| `qa/` | Session information, the installed commit, and the logs of the analysis and figure scripts |

The large fit objects (`data/fits.rds`, 53 MB) are not kept; `analysis.R`
writes them again.

## Rerunning

Requirements: R (at least 4.1) with CmdStan installed through cmdstanr, and
the packages listed in `qa/sessionInfo-analysis.txt` (bayesEfron 0.4.0,
cmdstanr, metafor, ebnm, loo, ggplot2, patchwork, ragg). From this folder:

```sh
bash scripts/install-package.sh path/to/bayesEfron   # a checkout of the tag v0.4.0
Rscript scripts/analysis.R          # the fits and comparisons; writes data/ (about 2.5 minutes)
Rscript scripts/stored-examples.R   # numbers read from the package's stored examples
Rscript scripts/tail-coverage.R     # reads data-raw/benchmark/ of the checkout (about 10 minutes)
Rscript scripts/figures.R           # figures/ and osm/figures/
Rscript scripts/tables.R            # tables/
```

With the same versions of CmdStan and of the packages on the same kind of
machine, the fits reproduce the draws exactly, so every number agrees; on
another machine or with another CmdStan version the draws differ and the
summaries agree to within Monte Carlo error, which is below the last digit
shown in the article for everything but the timings.
