# iNaturalist study data and provenance

## Source and sampling

The case study uses historical iNaturalist identification events for observations associated with Aves (birds). The original acquisition queried identification records at genus/family/order/class ranks for observations created during 2022-01-01 through 2024-12-31 and with a current species-level assignment. It paginated identification hits (200 records per page, up to 20 pages per rank), stopping upon reaching a candidate pool of at least 4,000 unique observations. The archived query-hit file includes 4,400 identification hits over 4,141 observations; only the genus and family searches were reached before stopping. Candidate observation identifiers were sorted and capped at 4,000 before detail retrieval.

The retained selected cohort required a research-grade current species assignment and a first **compatible coarse** identification. Within the frozen screening table, eligible species required at least 15 observations, and eligible genera at least two species. The eight most-supported genera were retained with at most five species each. This selection yields 695 observations, 22 species, and eight genera.

The post-hoc quality screen excluded a species if the compatible-event rate was below 85% or its withdrawn-identification rate exceeded 20%. It excluded 56 observations assigned to *Larus smithsonianus*, leaving 639 observations and 21 species in eight genera. The initial 695-observation results are retained as a sensitivity-to-screening comparison; neither cohort is a representative sample of the overall platform.

## Reconstructed time series

The fixed inferred tree is `ROOT -> genus -> species` (an induced label hierarchy, not the full biological taxonomic graph). A current reference species lineage is used to judge compatibility of past identification events. Compatible identifications, including those subsequently marked withdrawn in the frozen export, enter the retrospective reconstruction. Off-lineage events remain in the audit files but do not alter the monotone envelope.

Time zero is the earliest compatible coarse identification. States at 0, 1, 7, 30, and 90 elapsed days are the deepest compatible taxonomic node reached by each cutoff, collapsed into `ROOT`, a genus, or a species in the induced hierarchy. Thus the trajectories are **retrospective compatible envelopes**, not historical platform consensus states. Including compatible withdrawn events matters: 68 such event records occur in the screened cohort; excluding withdrawn records at the original time anchors changes 14 archived state cells.

The archived starting times are in early January 2022, and an event in the export is timestamped 2026-09-30. Thus the nominal 90-day follow-up is available to all selected observations. The precise API export timestamp is not independently recorded; the observed event timestamp only establishes a lower bound.

## Frozen inputs

`data/frozen/coarse_start_records.zip` contains:

- `coarse_identification_candidate_hits.csv`: saved identification candidate hits;
- `coarse_start_screening.csv`: screening table used for original observation selection;
- `identification_event_audit_initial.csv`, `identification_event_audit.csv`: initial and screened event histories with compatibility and current/withdrawn flags;
- `reconstructed_coarse_start_trajectories_initial.csv`, `reconstructed_coarse_start_trajectories.csv`: 695- and 639-observation reconstructed trajectories;
- recorded crossing masses, sensitivity results, and study summaries.

`data/frozen/cleaned_analysis_results.zip` and `data/frozen/cohort_comparison_results.zip` contain reference outputs for independent comparison of calculated intervals and study parameters.

## What is reproducible and what is not

The offline scripts reproduce the observation IDs selected **from the frozen screening table**, the quality exclusion, time anchors, all saved horizon states from frozen event audits, and the selected-cohort LP estimates. The scripts do not recreate the upstream historical API responses or the original candidate-pool acquisition from an independently saved raw API snapshot. Live retrieval can differ because identifiers, taxonomy, identifications, or withdrawal status may change. Results are conditional on the selected cohort and maintained sensitivity model.

## Use and citation

Data originate from the iNaturalist platform and contributors. Reference site: https://www.inaturalist.org/ ; API: https://api.inaturalist.org/v1/ . The publication's reference list should be used for full scholarly citation. Observation and identification identifiers link to public records, but this package does not establish license compatibility for every record's further redistribution. The source code's MIT License applies to code only. Verify current platform conditions, contributor licenses, and personal-data policies before publishing these processed files beyond this reproduction archive.
