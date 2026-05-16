# bayesEfron Export Ledger

This ledger records the active v0.1 public namespace for `bayesEfron`.
It follows the DEC-026 posture: user-facing producers and generics are exported,
S3 methods are registered, and constructors/validators remain internal.

## Active Exported Functions

| NAMESPACE directive | Role | Primary return |
|---|---|---|
| `export(as_bef_data)` | Input converter generic | `bef_data` |
| `export(bayes_efron_clear_cache)` | Cache maintenance utility | invisible integer counts |
| `export(bayes_efron_compile)` | CmdStan cache pre-warming entry | invisible `cmdstanr::CmdStanModel` |
| `export(bayes_efron_fit)` | Main fitting entry | `c("bef_fit_re", "bef_fit")` |
| `export(diagnose)` | Diagnostic producer generic | `bef_diagnostic` |
| `export(make_efron_grid)` | Grid-construction entry | plain grid list |

## Active S3 Method Registrations

### Input Conversion

| NAMESPACE directive | Status |
|---|---|
| `S3method(as_bef_data,default)` | active |
| `S3method(as_bef_data,escalc)` | active |
| `S3method(as_bef_data,list)` | active |

### `bef_fit` Parent Stratum

Family-agnostic methods safe for all future `bef_fit` children.

| NAMESPACE directive | Status |
|---|---|
| `S3method(format,bef_fit)` | active |
| `S3method(format,summary.bef_fit)` | active |
| `S3method(logLik,bef_fit)` | active |
| `S3method(nobs,bef_fit)` | active |
| `S3method(posterior::as_draws,bef_fit)` | active |
| `S3method(print,bef_fit)` | active |
| `S3method(print,summary.bef_fit)` | active |
| `S3method(summary,bef_fit)` | active |

### `bef_fit_re` Concrete Child Stratum

RE-specific methods that consume RE-only fields such as `theta_summary`,
`theta_rep_draws`, or grid-density draws.

| NAMESPACE directive | Status |
|---|---|
| `S3method(as.data.frame,bef_fit_re)` | active |
| `S3method(coef,bef_fit_re)` | active |
| `S3method(confint,bef_fit_re)` | active |
| `S3method(plot,bef_fit_re)` | active |
| `S3method(summary,bef_fit_re)` | active |
| `S3method(vcov,bef_fit_re)` | active |

### `bef_data` Standalone Class

| NAMESPACE directive | Status |
|---|---|
| `S3method(format,bef_data)` | active |
| `S3method(print,bef_data)` | active |
| `S3method(summary,bef_data)` | active |

### `bef_diagnostic` Standalone Class

| NAMESPACE directive | Status |
|---|---|
| `S3method(format,bef_diagnostic)` | active |
| `S3method(print,bef_diagnostic)` | active |
| `S3method(summary,bef_diagnostic)` | active |

### Diagnostic Producer Dispatch

`diagnose()` is a user-facing producer generic rather than part of the
20-dispatch fit/data/diagnostic method partition.

| NAMESPACE directive | Status |
|---|---|
| `S3method(diagnose,bef_fit)` | active |
| `S3method(diagnose,default)` | active |

## Internal Constructors and Validators

These names are intentionally not exported in v0.1.

| Internal name | Status |
|---|---|
| `new_bef_data` | internal |
| `new_bef_fit` | internal |
| `new_bef_fit_re` | internal |
| `new_bef_diagnostic` | internal |
| `validate_bef_data` | internal |
| `validate_bef_fit` | internal |
| `validate_bef_fit_re` | internal |
| `validate_bef_diagnostic` | internal |

## Reserved v0.2 Promotion Candidates

These names are reserved for possible future export if v0.2 introduces a
documented public construction or subclassing use case.

| Reserved name | Current status |
|---|---|
| `new_bef_fit` | internal at v0.1 |
| `new_bef_fit_he` | not implemented at v0.1 |
| `new_bef_fit_che` | not implemented at v0.1 |
| `new_bef_data` | internal at v0.1 |
| `new_bef_diagnostic` | internal at v0.1 |
| `validate_bef_fit` | internal at v0.1 |
| `validate_bef_fit_re` | internal at v0.1 |
| `validate_bef_fit_he` | not implemented at v0.1 |
| `validate_bef_fit_che` | not implemented at v0.1 |
| `validate_bef_data` | internal at v0.1 |
| `validate_bef_diagnostic` | internal at v0.1 |

## Count Summary

| Category | Count |
|---|---:|
| Active exported functions | 6 |
| Input-conversion S3 registrations | 3 |
| Structural method registrations | 20 |
| Diagnostic-producer S3 registrations | 2 |
| Total S3 registrations | 25 |

The structural 20-dispatch partition is `8 + 6 + 3 + 3`: eight methods on
`bef_fit`, six on `bef_fit_re`, three on `bef_data`, and three on
`bef_diagnostic`.
