#!/usr/bin/env python3
"""Recompute the iNaturalist full sensitivity path from the frozen selected cohort.

Uses one common minimum all-horizon relaxation at each sensitivity level for
both final-horizon and all-horizon interval calculations, as in the manuscript.
"""
from __future__ import annotations
from pathlib import Path
from zipfile import ZipFile
import json
import numpy as np
from scipy.optimize import linprog
from reproduce_frozen_inaturalist import read_archive_csv,make_problem,lp_matrices,HORIZONS

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/paper/frozen_data/sensitivity_path_recomputed.json'


def matrices(problem,horizons,gamma,final_day=90):
    # Existing LP builder fixes final-horizon compatibility to day 90.
    assert final_day == 90
    return lp_matrices(problem,horizons,gamma)


def endpoint_width(problem,horizons,gamma,eta=0):
    A,b=matrices(problem,horizons,gamma)
    n=len(problem[0]); n_compat=sum(map(len,problem[1].values()))
    if gamma is not None: b[n_compat:]+=eta
    kw=dict(A_ub=A,b_ub=b,A_eq=np.ones((1,n)),b_eq=[1.],bounds=[(0,None)]*n,method='highs')
    widths=[]
    for i in range(n):
        c=np.zeros(n); c[i]=1
        lo=linprog(c,**kw);hi=linprog(-c,**kw)
        if not(lo.success and hi.success):raise RuntimeError(f'Gamma {gamma}: failed leaf LP {i}: {lo.message} / {hi.message}')
        widths.append(-hi.fun-lo.fun)
    return float(np.mean(widths))


def min_relaxation(problem,gamma):
    A,b=matrices(problem,HORIZONS,gamma)
    n=len(problem[0]);ncompat=sum(map(len,problem[1].values()))
    R=np.c_[A,np.r_[np.zeros(ncompat),-np.ones(len(b)-ncompat)]]
    eq=np.r_[np.ones(n),0.][None,:]
    res=linprog(np.r_[np.zeros(n),1.],A_ub=R,b_ub=b,A_eq=eq,b_eq=[1.],
                bounds=[(0,None)]*(n+1),method='highs')
    if not res.success:raise RuntimeError(res.message)
    return max(0.,float(res.x[-1]))


def main():
    with ZipFile(ROOT/'data/frozen/coarse_start_records.zip') as z:
        clean=read_archive_csv(z,'reconstructed_coarse_start_trajectories.csv')
    problem=make_problem(clean)
    hw=endpoint_width(problem,(90,),None)
    results={}
    for gamma in (6.,8.,12.,20.):
        eta=min_relaxation(problem,gamma)
        guard=1e-12
        fw=endpoint_width(problem,(90,),gamma,eta+guard)
        aw=endpoint_width(problem,HORIZONS,gamma,eta+guard)
        results[str(int(gamma))]={
           'minimum_relaxation':eta,'raw_feasible':eta<=1e-10,
           'hierarchy_only_width':hw,'final_horizon_width_common_relaxation':fw,
           'all_horizon_width':aw,'tightening_vs_hierarchy':1-aw/hw,
           'tightening_vs_final':1-aw/fw}
    expected={6:(.000149707,.00566301,.00349597),
              8:(.0000366995,.005148,.003076),
              12:(0.,.00563620739,.00404156206),
              20:(0.,.006549,.005568)}
    checks={}
    for g,(eta,final,allv) in expected.items():
        row=results[str(g)]
        # Entries with rounded values in manuscript receive a wider tolerance.
        max_dev=max(abs(row['minimum_relaxation']-eta),abs(row['final_horizon_width_common_relaxation']-final),abs(row['all_horizon_width']-allv))
        checks[str(g)]={'max_absolute_difference_vs_published_precision':max_dev}
        assert max_dev<(2e-6 if g in (8,20) else 5e-8), (g,row,max_dev)
    output={'cohort_observations':len(clean),'horizons_days':list(HORIZONS),
            'method':'LP minimum relaxation of sibling sensitivity rows, common-eta endpoint comparison',
            'sensitivity_path':results,'comparisons_to_published_precision':checks}
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))

if __name__=='__main__':main()
