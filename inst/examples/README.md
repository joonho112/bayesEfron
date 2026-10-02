# Stored examples

`grid_comparison_example.rds` holds the fits that
`vignette("choosing-a-grid")` and `vignette("diagnostics")` discuss, as a
list:

- `raudenbush`, the result of `compare_efron_grids()` with its default
  arguments for the teacher expectancy studies (`raudenbush1985`), with
  `seed = 1985`;
- `raudenbush_sensitivity`, summaries of nine fits of the same studies, for
  three widths of the grid and six, eight and ten spline functions;
- `two_groups`, the result of `compare_efron_grids()` for simulated data
  with 200 sites whose true effects form two groups, with `seed = 2025`;
- `two_groups_runs`, diagnostics of three fits of those data: four spline
  functions with the default length of run and with a run four times as
  long, and six spline functions;
- `two_groups_sensitivity`, summaries of six fits of those data, for three
  widths of the grid and six and eight spline functions;
- `two_groups_data`, the simulated effects, estimates and standard errors.

The elements `cmdstan_version` and `created` say with which CmdStan and when
the file was made. The script that makes it is `data-raw/grid-comparison.R`
in the source repository.
