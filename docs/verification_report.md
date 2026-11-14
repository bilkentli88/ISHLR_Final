# Reproduction verification record

This report summarizes numerical checks against the frozen inputs and archived result files. It is not a statistical validation of the maintained identification assumptions.

## Synthetic benchmarks

- Primary experiment: 16 leaves, 12 internal nodes; 170 explicit flow variables versus 16 compact leaf variables and 135 compact affine rows. Archived compact and explicit leaf-bound optimizations agree to numerical precision.
- General-tree experiment: 24 tree fixtures and 240 random support objectives, with maximum archived discrepancy of approximately `1.11e-15`.
- Extended validation: 9 tree structures, 36 population scenarios, 72 sensitivity specifications, and 288 paired objectives. Independent rerun with HiGHS feasibility tolerances `1e-10` gives maximum objective discrepancy `2.22e-15`.
- Computational comparison: primary 16-leaf final-versus-full-history mean widths `0.11484075` and `0.10966039`. Separate fixtures scale from 16 to 256 leaves. Measured wall-clock times vary with the execution environment.

## Frozen iNaturalist study

- Initial: 695 observations, 22 species, 8 genera; screened: 639 observations, 21 species, 8 genera.
- The screening table reproduces all 695 initial observation identifiers. The quality screen excludes 56 *Larus smithsonianus* observations.
- Frozen compatible identification histories reproduce 3,475 initial and 3,195 screened annotation states, with no time-anchor mismatches.
- The raw Γ=12 LP is feasible for both cohorts. Unfiltered hierarchy/final/all mean widths are `0.02668410726`, `0.01615851365`, `0.01459302971`; screened mean widths are `0.00871898055`, `0.00563620739`, `0.00404156205`.
- Recomputed Γ=6/8/12/20 sensitivity path uses the same minimum all-horizon relaxation for final-horizon and all-horizon width comparisons. Reported widths and relaxations agree at the paper's stated precision.
- Reference frequencies within intervals at Γ=12 are 21/21 hierarchy-only, 10/21 final-horizon, and 7/21 all-horizon. Direct substitution of the reference vector violates 80 ordered sibling-comparison inequalities.

## Bootstrap stability

An independent implementation resampled 639 complete trajectories with replacement using NumPy's `default_rng(20261007)`, 200 resamples, and pseudo-final days 1, 7, and 90. The 600 replicate-level rows are archived under `results/paper/frozen_data/`.

| Pseudo-final day | Reference frequencies within intervals | Raw bootstrap feasibility | Median minimum relaxation | Median mean interval width |
|---|---:|---:|---:|---:|
| 1 | 19/21 | 100% | 0 | 0.02921149 |
| 7 | 17/21 | 90% | 0 | 0.01144046 |
| 90 | 7/21 | 37% | 0.00007285805 | 0.00410945 |

At day 90, the 95th percentile of minimum relaxation is `0.00198716206`, while the 2.5th and 97.5th percentiles of mean width are `0.00201140056` and `0.00761328391`. These agree with the published summary after rounding. The historical bootstrap notebook was not supplied, but these values were regenerated from the archived cohort with an independently written script using the documented protocol.

## Qualifications

The offline tests do not reproduce the original live-API candidate discovery from a full raw-response snapshot. Later reference species names are not independent population truth. Bootstrap widths are descriptive resampling diagnostics, not frequentist confidence intervals. None of the numerical checks proves the population identification assumptions.
