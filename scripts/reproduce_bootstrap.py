#!/usr/bin/env python3
"""Reproduce trajectory-level bootstrap stability statistics from frozen iNaturalist records.

Default: 200 trajectory resamples, NumPy default_rng seed 20261007,
Γ=12 and pseudo-final days 1, 7, 90. This script uses no API requests.
"""
from __future__ import annotations
import argparse
from pathlib import Path
from zipfile import ZipFile
import json
import numpy as np
import pandas as pd
from scipy.optimize import linprog
from reproduce_frozen_inaturalist import read_archive_csv, make_problem

ROOT = Path(__file__).resolve().parents[1]
DAYS=(0,1,7,30,90)
FINAL_DAYS=(1,7,90)
GAMMA=12.0
GUARD=1e-12


def prepared_data(df):
    leaves,children,vectors,_=make_problem(df)
    events={}
    for child in [ch for group in children.values() for ch in group]:
        possible = {child} if child.startswith('S:') else {child}|set(children[child])
        events[child]={h:df[f'state_{h}d'].isin(possible).to_numpy().astype(float) for h in DAYS}
    pairs=[(c,d) for group in children.values() for c in group for d in group if c!=d]
    return leaves,vectors,events,pairs


def lp_rows(prepared,indices,final_day):
    leaves,vectors,events,pairs=prepared
    included=[h for h in DAYS if h<=final_day]
    crossing={child:{h:events[child][h][indices].mean() for h in included} for child in events}
    compat_A=np.asarray([-vectors[child] for child in events]);compat_b=np.asarray([-crossing[child][final_day] for child in events])
    A=[];b=[]
    for h in included:
        for c,d in pairs:
            mc,md=crossing[c][h],crossing[d][h]
            A.append(mc*vectors[d]-GAMMA*md*vectors[c])
            b.append((1-GAMMA)*mc*md)
    return compat_A,compat_b,np.asarray(A),np.asarray(b)


def solve_endpoints(prepared,indices,final_day,with_endpoints=False):
    leaves=prepared[0];n=len(leaves)
    A,b,S,t=lp_rows(prepared,indices,final_day)
    relax_A=np.concatenate([np.column_stack([A,np.zeros(len(A))]),np.column_stack([S,-np.ones(len(S))])])
    relax_b=np.r_[b,t]
    c=np.r_[np.zeros(n),1.]
    one=np.r_[np.ones(n),0.][None,:]
    relax=linprog(c,A_ub=relax_A,b_ub=relax_b,A_eq=one,b_eq=[1.],bounds=[(0,None)]*(n+1),method='highs')
    if not relax.success:raise RuntimeError(f'Minimum relaxation LP failed: {relax.message}')
    eta=max(float(relax.x[-1]),0.0)
    opt_A=np.concatenate([A,S]);opt_b=np.r_[b,t+eta+GUARD]
    kwargs=dict(A_ub=opt_A,b_ub=opt_b,A_eq=np.ones((1,n)),b_eq=[1.],bounds=[(0,None)]*n,method='highs')
    lower=[];upper=[]
    for j in range(n):
        objective=np.zeros(n);objective[j]=1
        lo=linprog(objective,**kwargs);hi=linprog(-objective,**kwargs)
        if not (lo.success and hi.success):
            raise RuntimeError(f'Endpoint LP failed day={final_day}, leaf={j}: {lo.message} / {hi.message}')
        lower.append(lo.fun);upper.append(-hi.fun)
    mean_width=float(np.mean(np.asarray(upper)-np.asarray(lower)))
    if with_endpoints:return eta,mean_width,np.asarray(lower),np.asarray(upper)
    return eta,mean_width


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--resamples',type=int,default=200)
    ap.add_argument('--seed',type=int,default=20261007)
    ap.add_argument('--out',type=Path,default=ROOT/'results/paper/frozen_data')
    args=ap.parse_args()
    with ZipFile(ROOT/'data/frozen/coarse_start_records.zip') as z:
        df=read_archive_csv(z,'reconstructed_coarse_start_trajectories.csv')
    prepared=prepared_data(df)
    original=np.arange(len(df))
    ref=np.asarray([(df.species_id==int(l.split(':')[1])).mean() for l in prepared[0]])
    reference_counts={}
    for h in FINAL_DAYS:
        _,_,lo,hi=solve_endpoints(prepared,original,h,with_endpoints=True)
        reference_counts[str(h)]=int(((ref>=lo-1e-9)&(ref<=hi+1e-9)).sum())
    rng=np.random.default_rng(args.seed)
    records=[]
    for i in range(args.resamples):
        sample=rng.integers(0,len(df),size=len(df))
        for h in FINAL_DAYS:
            eta,width=solve_endpoints(prepared,sample,h)
            records.append({'replicate':i,'pseudo_final_day':h,'minimum_relaxation':eta,
                            'mean_width':width,'raw_feasible':eta<=1e-10})
        if (i+1)%50==0:
            print(f'Completed {i+1} of {args.resamples} complete-trajectory resamples',flush=True)
    frame=pd.DataFrame(records)
    summary={}
    for h in FINAL_DAYS:
        subset=frame.loc[frame.pseudo_final_day==h]
        summary[str(h)]={
            'reference_frequencies_within_intervals':reference_counts[str(h)],
            'raw_feasibility_fraction':float(subset.raw_feasible.mean()),
            'median_minimum_relaxation':float(subset.minimum_relaxation.median()),
            'p95_minimum_relaxation':float(subset.minimum_relaxation.quantile(.95)),
            'median_mean_width':float(subset.mean_width.median()),
            'p025_mean_width':float(subset.mean_width.quantile(.025)),
            'p975_mean_width':float(subset.mean_width.quantile(.975)),
        }
    result={'seed':args.seed,'resamples':args.resamples,'cohort_observations':len(df),
            'resampling_unit':'complete retrospective annotation trajectory',
            'gamma':GAMMA,'numerical_guard':GUARD,
            'raw_feasible_maximum_relaxation':1e-10,
            'pseudo_final_results':summary}
    if args.seed==20261007 and args.resamples==200:
        assert [reference_counts[str(h)] for h in FINAL_DAYS]==[19,17,7]
        assert [summary[str(h)]['raw_feasibility_fraction'] for h in FINAL_DAYS]==[1.0,.90,.37]
        assert abs(summary['90']['p95_minimum_relaxation']-.00198716)<1e-7
        assert abs(summary['90']['median_minimum_relaxation']-7.286e-5)<1e-8
        assert abs(summary['90']['median_mean_width']-.0041094)<1e-7
        assert abs(summary['90']['p025_mean_width']-.00201140)<1e-7
        assert abs(summary['90']['p975_mean_width']-.00761328)<1e-7
    args.out.mkdir(parents=True,exist_ok=True)
    frame.to_csv(args.out/'bootstrap_200_replicates.csv',index=False)
    (args.out/'bootstrap_200_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
