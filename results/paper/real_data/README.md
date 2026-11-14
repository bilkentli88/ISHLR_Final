# iNaturalist results associated with the manuscript

The files in this directory are **archived numerical results**, not the primary reproduction inputs. The underlying processed records are included under `data/frozen/` and can be checked without internet access using `python scripts/reproduce_frozen_inaturalist.py`.

The analysis uses iNaturalist identification histories, not the iNaturalist computer-vision classification benchmark. Frozen events include compatible identifications marked withdrawn at export time; incompatible/off-lineage events are retained in audits but omitted from the monotone reconstructed states. A later live API query does not recreate the original historical selection.

- `analysis_configuration.json`: study parameters and reconstructed cohort criteria.
- `reported_species_resolution.csv`: reference-species resolution percentages.
- `reported_sensitivity_table.csv`: archived values for the sensitivity path, including numerical relaxations.
- `reported_real_data_summary.json`: cleaned/initial cohort metrics and feasibility thresholds.
- `reported_unfiltered_robustness.json`: pre-screening versus screened result.
- `reported_sensitivity_path.pdf` and `.png`: manuscript sensitivity-path graphic.

For the precise frozen-data acquisition and auditing scope, see `docs/data_provenance.md`. Some reported summary figures are rounded, and small solver-dependent variations in estimated feasibility thresholds are expected.
