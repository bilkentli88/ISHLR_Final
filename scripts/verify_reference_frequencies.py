#!/usr/bin/env python3
"""Audit current-reference species frequencies against selected-cohort sensitivity LPs.

Read frozen cohort trajectories, not the live iNaturalist API.
"""
from pathlib import Path
from zipfile import ZipFile
import json
import numpy as np
from reproduce_frozen_inaturalist import read_archive_csv, make_problem, lp_matrices, bound, HORIZONS

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/paper/frozen_data/reference_frequency_audit.json'


def main():
    with ZipFile(ROOT/'data/frozen/coarse_start_records.zip') as z:
        trajectory = read_archive_csv(z,'reconstructed_coarse_start_trajectories.csv')
    problem = make_problem(trajectory)
    leaves = problem[0]
    p_ref = np.array([(trajectory.species_id==int(l.split(':')[1])).mean() for l in leaves])
    assert len(leaves) == 21
    runs = [('hierarchy_only',(90,),None),('final_horizon',(90,),12.0),('all_horizons',HORIZONS,12.0)]
    comparisons={}
    for name,hs,gamma in runs:
        lower,upper=bound(problem,hs,gamma)
        within=(p_ref>=lower-1e-9)&(p_ref<=upper+1e-9)
        max_gap=float(np.maximum.reduce([lower-p_ref,p_ref-upper,np.zeros(len(leaves))]).max())
        comparisons[name]={'within':int(within.sum()),'outside':int((~within).sum()),'max_endpoint_gap':max_gap}
    A,b=lp_matrices(problem,HORIZONS,12.0)
    residual=np.asarray(A@p_ref-b)
    count=int(np.count_nonzero(residual>1e-12))
    result={'cohort_observations':len(trajectory),'reference_species':len(leaves),
            'gamma':12,'tolerance':1e-9,'interval_comparisons':comparisons,
            'ordered_sibling_constraints_violated':count,
            'max_positive_lp_residual':float(np.maximum(residual,0).max()),
            'note':'Reference assignments are within-cohort counts, not independent population probabilities or confidence coverage.'}
    assert [comparisons[x]['within'] for x in ['hierarchy_only','final_horizon','all_horizons']] == [21,10,7]
    assert abs(comparisons['all_horizons']['max_endpoint_gap']-0.0010671656)<1e-8
    assert count==80
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
