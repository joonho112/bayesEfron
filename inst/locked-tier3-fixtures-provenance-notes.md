# Locked Tier 3 Fixture Provenance Notes

This note explains the Tier 3 Lee-Sui fixture SHA transitions relevant to the
strict v0.1.0 release audit. The authoritative current fixture checksum
manifest remains `inst/locked-tier3-fixtures-checksums.txt`.

## Current Locked Fixture SHAs

The current locked Tier 3 fixtures are recorded in the four-column TSV manifest
`inst/locked-tier3-fixtures-checksums.txt`.

| K | Fixture path | Current SHA-256 |
|---:|---|---|
| 50 | `tests/testthat/_fixtures/lee_sui_K50.rds` | `5278df09a077fae3b7b7a0764bd4b45416ae2085ce2339b73b7d3ddc6b809d4f` |
| 100 | `tests/testthat/_fixtures/lee_sui_K100.rds` | `e538319a47ed14d412c032d9f34269f75b9ab3a676ef9aad8321859776521fdb` |
| 200 | `tests/testthat/_fixtures/lee_sui_K200.rds` | `951a0510d19e41655f6e6dff5bd9736c3683917f6b0e7eeac5b47d2bfb64277e` |
| 500 | `tests/testthat/_fixtures/lee_sui_K500.rds` | `c5cdae545e9ecdd9f0374bd1985ef61abbd1ba7ccdb7d695011fc27dcb943b25` |
| 1500 | `tests/testthat/_fixtures/lee_sui_K1500.rds` | `b276c3056ae81dce7b70cbf2d72f719e17842bca96c933c9338e54a7075e7843` |

## Regeneration Timeline

The original single-replication Tier 3 fixtures were regenerated into v2
all-20-replication fixtures while preserving manifest target IDs and
relative paths. The canonical generation date is `generated_on =
"2026-05-13"` so regeneration does not drift by calendar date.

Accepted coverage evidence against the pre-repair v2 fixture bytes:

| K | Coverage | Replications fit | Failed replications | In band |
|---:|---:|---:|---:|---|
| 50 | 0.8940 | 20 | 0 | TRUE |
| 100 | 0.8910 | 20 | 0 | TRUE |
| 200 | 0.8940 | 20 | 0 | TRUE |
| 500 | 0.8916 | 20 | 0 | TRUE |

A subsequent regeneration pass repaired the K1500 fixture identity problem
and regenerated all five Tier 3 fixtures. For K50, K100, K200, and K500
this was a metadata and policy hardening pass: the fixtures remained
canonical to the archived Lee-Sui appendix object, retained 20 unique
input payloads per K, and gained the explicit
`replication_identity_policy = "unique_input_payload_hashes"` metadata
used by the hardened validators. The file SHAs changed because the RDS
metadata changed; this does not indicate a change in the K50-K500
archived appendix data-vector payloads used for the accepted coverage
evidence.

The K1500 case is different. The earlier K1500 attempt is
diagnostic-only because the archived appendix K1500 rows were
label-distinct but data-identical. The repair regenerated K1500 using
the public Part 01 full-panel simulation rule for data vectors while
retaining the archived appendix dataset and seed controls. Accepted
K1500 evidence against the repaired fixture:

| K | Coverage | Replications fit | Failed replications | In band |
|---:|---:|---:|---:|---|
| 1500 | 0.876433333333333 | 20 | 0 | TRUE |

## Canonical Input Sources

Current fixture metadata records the canonical input source for each K:

| K | `metadata$source$canonical_input` | Notes |
|---:|---|---|
| 50 | `archived_appendix_small_k_result` | Archived Lee-Sui appendix K-specific panels. |
| 100 | `archived_appendix_small_k_result` | Archived Lee-Sui appendix K-specific panels. |
| 200 | `archived_appendix_small_k_result` | Archived Lee-Sui appendix K-specific panels. |
| 500 | `archived_appendix_small_k_result` | Archived Lee-Sui appendix K-specific panels. |
| 1500 | `part01_full_panel_simulation_rule` | Repaired full-panel data vectors with archived K1500 dataset/seed controls. |

This distinction is intentional. Do not describe K1500 as if the archived
appendix supplied 20 independent K1500 data panels; it did not. Do not describe
K50-K500 as public Part 06 replay fixtures; their canonical payload source is
the archived appendix object.

## Final Review Disposition

The final release review flagged an audit-trail concern: earlier
K50-K500 coverage runs cite fixture SHAs that no longer match the current
locked fixture file SHAs. This note resolves the interpretation of that
mismatch.

For K50, K100, K200, and K500, the SHA mismatch is explained by v2
metadata and replication-identity-policy hardening after the coverage
runs, not by an intended change to the archived appendix data-vector
payloads. Active test coverage for the current fixture
payload/provenance invariants in `tests/testthat/` checks this
explanation mechanically.

