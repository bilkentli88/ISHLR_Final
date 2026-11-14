# Informatively Stopped Hierarchical Label Refinement (ISHLR)

Reproduction materials for **“Informatively Stopped Hierarchical Label Refinement: Sharp Identification and Sensitivity Analysis.”**

This repository contains numerical experiments, the associated result tables, and frozen processed iNaturalist records needed to reproduce the manuscript's selected-cohort analyses without querying the live iNaturalist API.

## Quick start

Use Python 3.10+ (the archived runs used Python 3.13.5). From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate          # Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -r requirements.txt
python scripts/verify_repository.py
```

`verify_repository.py` checks the archived paper results, runtime tables, extended synthetic output, and the historical iNaturalist cohort reconstruction. This entry point **does not access the internet** and does not alter archived benchmarks other than refreshing the offline verification JSON files. The longer 200-resample bootstrap and full synthetic benchmark can be rerun using their separate commands. Its verification scope is described in `docs/reproduction_map.md`.

## Reproduce the main numerical findings

| Manuscript result | Command or notebook | Reproduced evidence |
|---|---|---|
| Compact LP versus explicit latent-flow LP; primary 16-leaf synthetic fixture and 256-leaf scaling | `notebooks/paper/03_global_compatibility_gate.ipynb` | `results/paper/global_compatibility/` |
| 24-tree constructive coupling and historical nonredundancy | `notebooks/paper/04_general_tree_theorem_gate.ipynb` | `results/paper/general_tree/` |
| Sample-analogue endpoint convergence (50 runs per sample size) | `notebooks/paper/05_sample_analog_consistency_gate.ipynb` | `results/paper/sample_analog/` |
| Additional depth/horizon/class-imbalance validation (72 configurations, 288 paired LPs) | `python scripts/benchmark_extended_synthetic.py` | `results/paper/extended_synthetic/` |
| Final versus full-history timing, primary and 16–256-leaf scaling | `python scripts/benchmark_horizon_primary.py` and `python scripts/benchmark_horizon_scaling.py` | `results/paper/horizon_runtime/` |
| iNaturalist initial/screened cohorts, state reconstruction, and Γ=12 endpoints | `python scripts/reproduce_frozen_inaturalist.py` | `results/paper/frozen_data/` |
| iNaturalist full Γ sensitivity path with common relaxation | `python scripts/reproduce_sensitivity_path.py` | `results/paper/frozen_data/sensitivity_path_recomputed.json` |
| iNaturalist reference-frequency interval comparison and direct sensitivity violations | `python scripts/verify_reference_frequencies.py` | `results/paper/frozen_data/reference_frequency_audit.json` |
| iNaturalist 200-trajectory bootstrap at pseudo-final days 1, 7, and 90 | `python scripts/reproduce_bootstrap.py` | `results/paper/frozen_data/bootstrap_200_summary.json` and resample-level CSV |

For the synthetic notebooks, use Jupyter and execute all cells from the repository root. The three notebooks are the published experiment sources; outputs supplied under `results/paper` are the archived reference runs. Re-execution may overwrite local result files and create additional plots or archives. Copy the repository if you want to preserve the reference outputs.

Reproduce the extended synthetic benchmark without overwriting archived results:

```bash
python scripts/benchmark_extended_synthetic.py --out /tmp/ishlr_extended_check
```

For runtime scripts, set `ISHLR_OUTPUT_DIR` to a new directory if you wish to preserve the provided timing CSV files. Median times are machine-dependent; model dimensions and deterministic LP objective values are the comparison targets.

## iNaturalist data: which workflow is authoritative?

**The offline frozen-data workflow is authoritative for reproducing the published cohort.** It reconstructs the selection from a saved screening table (695 observations), applies the exploratory species-quality exclusion (639 observations), checks all 3,475 and 3,195 annotation states against the archived event audits, and recomputes the unrelaxed `Gamma=12` LP endpoints. All compatible identification events are included in the retrospective monotone envelope even when marked withdrawn in the archived export.

Historical acquisition notebooks are preserved under `notebooks/historical_provenance/` to document how records were originally requested and processed. They are not the default reproduction workflow; they may require the internet, original user interaction, and access to historical API responses no longer guaranteed to exist. The initial 4,000 candidate-observation acquisition process **cannot be independently replayed from the provided processed records alone**. The included event histories are not platform consensus trajectories, and reference species names are not independently measured ground-truth population labels.

For the archive contents, follow-up availability, temporal horizons, API-selection limits, and data-use qualifications, see [`docs/data_provenance.md`](docs/data_provenance.md).

## Contents

```text
notebooks/paper/                Primary synthetic and sample-analogue experiments
notebooks/historical_provenance/ Historical acquisition/validation notebooks (archival)
data/frozen/                   Frozen processed iNaturalist study files (ZIP)
scripts/                       Standalone reproduction and verification scripts
results/paper/                 Archived manuscript experiment results and new audits
docs/                          Protocols, source provenance, scope and limitations
requirements.txt               Dependencies
MANIFEST.sha256                Integrity checksums for repository files
```

## Reproducibility scope and limitations

- The compact-versus-explicit optimization results are **numerical experiments**, not a substitute for mathematical proofs.
- iNaturalist studies use a targeted, retrospectively screened cohort; neither the cohort nor the reference species frequencies are representative probability samples of all iNaturalist observations.
- Minimal LP sensitivity relaxation is a feasibility adjustment, not a confidence radius. Neither the reported bounds nor the archived bootstrap diagnostics are calibrated confidence intervals.
- The full original API-response snapshot and historical bootstrap notebook are not included in the supplied source archives. The bootstrap experiment has been independently reimplemented from the documented protocol and the frozen trajectories, and its numerical summaries have been verified against the published results. See [`docs/reproduction_map.md`](docs/reproduction_map.md).

## Licensing and data credit

The original source code is made available under the MIT License (`LICENSE`). iNaturalist observations and identification metadata originate from their respective platform records and may be subject to different licenses and reuse terms; the repository's code license does **not** relicense the underlying data. Please consult `docs/data_provenance.md` before further redistribution of the frozen processed data.
