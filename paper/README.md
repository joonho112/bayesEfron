# Reproduction material for the SoftwareX article on bayesEfron 0.4.0

This folder contains the scripts, data and recorded outputs behind
"bayesEfron: Fully Bayesian inference with Efron’s log-spline prior in R,"
prepared for SoftwareX. The article and its online supplementary material
(OSM) are submitted as PDFs. The statistical results can be reproduced from
this folder without their LaTeX sources. This folder is excluded from the
R package build by `.Rbuildignore`.

The article's Data availability statement and software metadata (S3) link to
the commit containing these reproduction files. Use that commit when reproducing
the submitted article. The package itself remains at tag `v0.4.0`; that release
tag does not identify the revision of this folder.

## Contents

| Folder or file | Content |
|:--|:--|
| `scripts/install-package.sh` | Installs from a supplied checkout of bayesEfron 0.4.0 into this folder's `r-library/` and records the installed commit |
| `scripts/analysis.R` | Fits for the two examples and the 19-study grid comparison; empirical Bayes comparisons, numerical registry and console output |
| `scripts/stored-examples.R` | Numerical quotations from the package's stored examples |
| `scripts/tail-coverage.R` | Subset coverage in the two-group designs, using the saved benchmark results and the continuous generating distribution |
| `scripts/ratio-mcse.R` | Paired delta-method Monte Carlo standard errors for the length and RMSE ratios in Figure 2; reads saved benchmark results without refitting |
| `scripts/figures.R` | Figures 2–5 and both supplementary figures from saved data; no model fitting |
| `figures/fig1_architecture.tex` | Standalone TikZ source for Figure 1, compiled separately |
| `scripts/tables.R`, `scripts/table-layout.R` | Table contents generated from saved data, with the final OSM layout adjustments |
| `scripts/check-numbers.R` | In reproduction mode, compares rerun registries with archived values without manuscript sources |
| `data/numbers.csv`, `data/stored-examples-numbers.csv`, `data/tail-coverage-numbers.csv` | Registries of selected numerical quotations; these do not enumerate every numeric literal in the article |
| `data/reference/` | Preserved article registries with the same three filenames, used as the reference for a rerun |
| `data/tail-coverage.csv`, `data/tail-coverage-panels.csv` | Aggregate and panel-level subset coverage; panel rows support the replication-based MC SE calculation |
| `data/ratio-mcse.csv` | Plotted ratios and paired MC SEs, so figures do not need a source checkout to read benchmark results |
| `data/krw.rds`, `data/twin.rds`, `data/grid-comparison.rds`, `data/evidence.rds` | Summaries, curves and simulation tables used by the figures and tables |
| `data/console/` | Recorded console sessions |
| `figures/`, `osm/figures/` | Figure PDFs; `fig1_architecture` is Figure 1, `fig4_evidence` Figure 2, `fig2_firms` Figure 3, `fig3_twin_towers` Figure 4 and `fig5_grid_and_uncertainty` Figure 5 |
| `tables/` | Generated LaTeX table fragments; document preambles and surrounding table layouts are maintained with the manuscript sources |
| `qa/` | Recorded session information, installed commit, and analysis/figure logs |

The large `data/fits.rds` is omitted. `analysis.R` regenerates it; figure and
table regeneration does not require it. Keep `data/reference/` unchanged when
rerunning; replacing the reference with new results would defeat the comparison.

## Requirements and package version

Use R 4.1 or later and the dependencies recorded in
`qa/sessionInfo-analysis.txt`. The analysis uses bayesEfron 0.4.0, cmdstanr,
metafor, ebnm and loo; plotting also uses ggplot2, patchwork, scales and ragg.
CmdStan is needed for `analysis.R`, but not to redraw saved figures or tables.
The recorded fitting toolchain was R 4.6.0, cmdstanr 0.8.0 and CmdStan 2.38.0.

Use a source checkout at tag `v0.4.0` (commit `89a0a71`) for installation and
benchmark inputs. The reproduction folder can be a later revision than that
release; use a separate tagged checkout so switching the package checkout to
the tag does not remove this folder. Run all commands below from this folder:

```sh
BAYESEFRON_SOURCE=/absolute/path/to/bayesEfron-v0.4.0
bash scripts/install-package.sh "$BAYESEFRON_SOURCE"
```

The scripts look first in the local `r-library/`. The install command concerns
bayesEfron; install its other requirements beforehand.

## Redraw saved figures and tables without fitting

```sh
Rscript scripts/figures.R
Rscript scripts/tables.R
```

The figure script chooses `cairo_pdf` when Cairo is loadable and falls back to
base R's `pdf` device otherwise; PNG proofs use ragg. Figure 1 is separate and
requires a LaTeX installation with `latexmk`, the standalone class and TikZ:

```sh
(cd figures && latexmk -pdf -interaction=nonstopmode -halt-on-error fig1_architecture.tex)
```

## Recompute and compare the statistical results

```sh
Rscript scripts/tail-coverage.R "$BAYESEFRON_SOURCE" --check-inputs
Rscript scripts/ratio-mcse.R "$BAYESEFRON_SOURCE" --check-inputs
Rscript scripts/analysis.R
Rscript scripts/stored-examples.R
Rscript scripts/tail-coverage.R "$BAYESEFRON_SOURCE"
Rscript scripts/ratio-mcse.R "$BAYESEFRON_SOURCE"
Rscript scripts/figures.R
Rscript scripts/tables.R
Rscript scripts/check-numbers.R --mode=reproduction --reference=data/reference
```

The tail and ratio scripts accept the package checkout as their first argument.
When this folder is directly inside a suitable checkout, their default first
checks the parent; explicit paths also work for a separately downloaded copy of this folder.
`--check-inputs` verifies required benchmark files without fitting or writing
results. The analysis refits the article examples and reads the package's
stored simulation results; it does not rerun the full benchmark campaign.

Reproduction mode compares the three regenerated registries with their preserved
counterparts in `data/reference/` and works without manuscript sources. Runtime
changes and registered software-version differences are reported separately from
statistical discrepancies. The criterion is
`abs(new - recorded) <= 1e-8 * max(1, abs(recorded))`, allowing floating-point
roundoff when replaying the recorded calculation; it is not a Monte Carlo error
bound or an accuracy guarantee for all printed digits. The article's fits
reproduced their stored draws on the recorded machine and toolchain. Other
platforms or versions can give different draws, requiring an assessment of
summary differences against Monte Carlo error.

The script also retains an archival mode for checking registered quotations
against manuscript sources during document preparation. That mode requires
`main.tex` and `osm/appendices/`, which are not included here. Use the reproduction
command above for this folder. Neither mode checks every numeric literal in
the manuscript.
