#!/usr/bin/env python3
"""Check deterministic saved summaries for the extended synthetic and bootstrap runs."""
from pathlib import Path
import json
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
PAPER=ROOT/'results/paper'

def load(path): return json.loads(path.read_text())

ext=load(PAPER/'extended_synthetic/extended_synthetic_summary.json')
assert ext['base_scenarios']==36
assert ext['distinct_tree_structures']==9
assert ext['configurations']==72
assert ext['support_comparisons']==288
assert ext['depths']==[2,3,4]
assert ext['horizon_counts']==[2,5]
assert ext['gamma_values']==[3.0,6.0]
assert ext['leaf_count_max']==43
assert ext['max_flow_variables']==943
assert abs(ext['max_objective_difference']-2.220446049250313e-15)<1e-18
source=pd.read_csv(PAPER/'extended_synthetic/extended_synthetic_stress.csv')
assert len(source)==72 and source.support_objectives_checked.sum()==288
assert source.max_compact_flow_objective_gap.max()<1e-12

boot=load(PAPER/'frozen_data/bootstrap_200_summary.json')
assert boot['resamples']==200 and boot['seed']==20261007
counts=[boot['pseudo_final_results'][str(i)]['reference_frequencies_within_intervals'] for i in (1,7,90)]
feasible=[boot['pseudo_final_results'][str(i)]['raw_feasibility_fraction'] for i in (1,7,90)]
assert counts==[19,17,7] and feasible==[1.0,0.9,0.37]
assert abs(boot['pseudo_final_results']['90']['median_minimum_relaxation']-7.286e-5)<1e-8
assert abs(boot['pseudo_final_results']['90']['p95_minimum_relaxation']-.00198716)<1e-7
assert abs(boot['pseudo_final_results']['90']['median_mean_width']-.0041094)<1e-7
samples=pd.read_csv(PAPER/'frozen_data/bootstrap_200_replicates.csv')
assert len(samples)==600 and samples.replicate.nunique()==200
assert set(samples.pseudo_final_day)=={1,7,90}
print('Extended synthetic and 200-resample bootstrap records verified.')
