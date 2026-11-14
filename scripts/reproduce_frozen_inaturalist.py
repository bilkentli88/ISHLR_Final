#!/usr/bin/env python3
"""Offline reproduction of iNaturalist cohort selection, trajectory states, and gamma=12 LP bounds.

Requires only frozen archived CSV files included in this package; no API calls.
Run: python scripts/reproduce_frozen_inaturalist.py
"""
from __future__ import annotations

import io
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
from scipy.optimize import linprog

BASE = Path(__file__).resolve().parents[1]
HORIZONS = (0, 1, 7, 30, 90)
GAMMA = 12.0


def read_archive_csv(z: ZipFile, member: str) -> pd.DataFrame:
    return pd.read_csv(io.BytesIO(z.read(member)))


def selection_from_screen(screen: pd.DataFrame) -> set[int]:
    eligible = screen.loc[screen.coarse_start].copy()
    counts = eligible.groupby(['genus_id', 'genus_name', 'species_id', 'species_name']).size().reset_index(name='n')
    species = counts.loc[counts.n >= 15]
    genus = species.groupby(['genus_id','genus_name']).agg(nspecies=('species_id', 'nunique'), nobs=('n', 'sum')).reset_index()
    chosen_genus = genus.loc[genus.nspecies >= 2].sort_values('nobs', ascending=False).head(8)
    chosen_species = set()
    for gid in chosen_genus.genus_id:
        chosen_species.update(int(s) for s in species.loc[species.genus_id == gid].sort_values('n',ascending=False).head(5).species_id)
    return set(eligible.loc[eligible.species_id.isin(chosen_species), 'observation_id'].astype(int))


def verify_states(traj: pd.DataFrame, events: pd.DataFrame) -> tuple[int,int]:
    events = events.copy()
    events['when'] = pd.to_datetime(events.created_at, utc=True, format='mixed')
    groups = {int(oid):grp.loc[grp.compatible].sort_values('when') for oid,grp in events.groupby('observation_id')}
    bad_states = 0
    bad_starts = 0
    for row in traj.itertuples(index=False):
        compatible = groups.get(int(row.observation_id))
        if compatible is None or compatible.empty:
            raise AssertionError('No compatible events for observation '+str(row.observation_id))
        t0 = compatible.when.iloc[0]
        bad_starts += int(t0 != pd.to_datetime(row.first_compatible_time, utc=True))
        for day in HORIZONS:
            seen = compatible.loc[compatible.when <= t0+pd.Timedelta(days=day)]
            state = ('S:'+str(row.species_id) if (seen.taxon_id == row.species_id).any()
                     else ('G:'+str(row.genus_id) if (seen.taxon_id == row.genus_id).any() else 'ROOT'))
            bad_states += int(state != getattr(row, f'state_{day}d'))
    return bad_states, bad_starts


def make_problem(traj: pd.DataFrame):
    species = traj[['species_id','genus_id']].drop_duplicates().sort_values(['genus_id','species_id'])
    leaves = [f'S:{int(x)}' for x in species.species_id]
    genera = sorted(species.genus_id.unique())
    children = {'ROOT':[f'G:{int(x)}' for x in genera]}
    for gid, sub in species.groupby('genus_id'):
        children[f'G:{int(gid)}'] = [f'S:{int(x)}' for x in sub.species_id]
    vectors = {leaf:np.eye(len(leaves))[i] for i,leaf in enumerate(leaves)}
    for gid, sub in species.groupby('genus_id'):
        vectors[f'G:{int(gid)}'] = sum((vectors[f'S:{int(x)}'] for x in sub.species_id), np.zeros(len(leaves)))
    crossing = {}
    for parent, kids in children.items():
        for ch in kids:
            possible = {ch} if ch.startswith('S:') else {ch} | set(children[ch])
            crossing[ch] = {day:float(traj[f'state_{day}d'].isin(possible).mean()) for day in HORIZONS}
    return leaves, children, vectors, crossing


def lp_matrices(problem, horizons: tuple[int, ...] = HORIZONS, gamma: float | None = GAMMA):
    leaves, children, vectors, crossing = problem
    A=[]; b=[]
    for ch in [child for kids in children.values() for child in kids]:
        A.append(-vectors[ch]); b.append(-crossing[ch][90])
    if gamma is not None:
        for kids in children.values():
            for day in horizons:
                for c in kids:
                    for d in kids:
                        if c == d: continue
                        mc, md = crossing[c][day], crossing[d][day]
                        A.append(mc*vectors[d] - gamma*md*vectors[c]); b.append((1-gamma)*mc*md)
    return np.array(A, dtype=float), np.array(b, dtype=float)


def bound(problem, horizons: tuple[int, ...], gamma: float | None = GAMMA):
    leaves=problem[0]; n=len(leaves)
    A,b=lp_matrices(problem,horizons,gamma)
    lower,upper=np.zeros(n),np.zeros(n)
    for j in range(n):
        c=np.zeros(n);c[j]=1
        kw=dict(A_ub=A,b_ub=b,A_eq=np.ones((1,n)),b_eq=np.array([1.0]),bounds=[(0,None)]*n,method='highs')
        lo=linprog(c,**kw);hi=linprog(-c,**kw)
        if not lo.success or not hi.success:raise RuntimeError(f'LP failure on leaf {j}: {lo.message}, {hi.message}')
        lower[j],upper[j]=lo.fun,-hi.fun
    return lower,upper


def feasible(problem, horizons: tuple[int,...],gamma:float)->bool:
    n=len(problem[0]);A,b=lp_matrices(problem,horizons,gamma)
    return bool(linprog(np.zeros(n), A_ub=A,b_ub=b,A_eq=np.ones((1,n)),b_eq=[1.0],bounds=[(0,None)]*n,method='highs').success)


def gamma_threshold(problem, horizons:tuple[int,...]):
    if feasible(problem,horizons,1.0):return 1.0
    lo=1.0;hi=2.0
    while not feasible(problem,horizons,hi):lo,hi=hi,hi*2
    while hi-lo>1e-6*max(1.0,hi):
        mid=(hi+lo)/2
        if feasible(problem,horizons,mid):hi=mid
        else:lo=mid
    return hi


def main():
    with ZipFile(BASE/'data'/'frozen'/'coarse_start_records.zip') as z:
        initial=read_archive_csv(z,'reconstructed_coarse_start_trajectories_initial.csv')
        clean=read_archive_csv(z,'reconstructed_coarse_start_trajectories.csv')
        initial_events=read_archive_csv(z,'identification_event_audit_initial.csv')
        events=read_archive_csv(z,'identification_event_audit.csv')
        audit=read_archive_csv(z,'species_quality_audit.csv')
        screen=read_archive_csv(z,'coarse_start_screening.csv')
        candidates=read_archive_csv(z,'coarse_identification_candidate_hits.csv')
    with ZipFile(BASE/'data'/'frozen'/'cohort_comparison_results.zip') as z:
        archive_results=read_archive_csv(z,'unfiltered_vs_cleaned_gamma_profile.csv')

    assert len(initial)==695 and initial.species_id.nunique()==22 and initial.genus_id.nunique()==8
    assert len(clean)==639 and clean.species_id.nunique()==21 and clean.genus_id.nunique()==8
    assert selection_from_screen(screen)==set(initial.observation_id.astype(int)), 'Selection differs from saved initial cohort'
    excluded=audit.loc[(audit.compatible_event_rate < .85)|(audit.withdrawn_event_rate > .20)]
    assert len(excluded)==1 and excluded.species_id.iloc[0]==1578489
    assert set(clean.observation_id)==set(initial.loc[initial.species_id!=1578489,'observation_id'])
    assert len(events[(events.compatible)&(~events.current)])==68
    assert len(candidates)==4400 and candidates.observation_id.nunique()==4141
    verifications={}
    for name,df,ev in [('initial_695',initial,initial_events),('cleaned_639',clean,events)]:
        mismatches,start_mismatches=verify_states(df,ev)
        assert (mismatches,start_mismatches)==(0,0)
        verifications[name]={'observations':len(df),'species':int(df.species_id.nunique()),'genera':int(df.genus_id.nunique()),'checked_state_cells':len(df)*len(HORIZONS),'state_mismatches':mismatches,'start_mismatches':start_mismatches}
    outputs={}
    for key,df in [('unfiltered_695',initial),('cleaned_639',clean)]:
        problem=make_problem(df)
        h=bound(problem,(90,),None); f=bound(problem,(90,)); a=bound(problem,HORIZONS)
        hw=float(np.mean(h[1]-h[0]));fw=float(np.mean(f[1]-f[0]));aw=float(np.mean(a[1]-a[0]))
        metrics={'hierarchy_mean_width':hw,'final_mean_width':fw,'all_mean_width':aw,
                 'reduction_vs_hierarchy':float(1-aw/hw),'reduction_vs_final':float(1-aw/fw),
                 'leaf_intervals_tightened':int(np.sum((a[1]-a[0])<(f[1]-f[0])-1e-10)),
                 'raw_feasible_at_gamma12':feasible(problem,HORIZONS,12.0),
                 'min_full_history_gamma':float(gamma_threshold(problem,HORIZONS)),
                 'min_final_horizon_gamma':float(gamma_threshold(problem,(90,)))}
        row=archive_results.loc[(archive_results.cohort==key)&(archive_results.gamma==12)].iloc[0]
        pairs=[('hierarchy_mean_width','hierarchy_mean_width'),('final_mean_width','final_mean_width_common_eta'),('all_mean_width','all_mean_width_common_eta'),('reduction_vs_hierarchy','tightening_vs_hierarchy'),('reduction_vs_final','tightening_vs_final')]
        for ours,ref in pairs:
            assert abs(metrics[ours]-float(row[ref]))<2e-8, f'{key} {ours}: {metrics[ours]} vs {row[ref]}'
        outputs[key]=metrics
    summary={'verification':'PASS','method':'Recomputed from frozen screening, event, and trajectory CSV records, without live API access','selection_observation_ids_reproduced':len(initial), 'quality_excluded_species':'Larus smithsonianus','quality_excluded_observations':56,'candidate_identification_hits':len(candidates),'unique_candidate_observations':int(candidates.observation_id.nunique()),'compatible_withdrawn_events_cleaned':68, 'trajectory_checks':verifications,'gamma12':outputs}
    (BASE/'results'/'paper'/'frozen_data').mkdir(parents=True, exist_ok=True)
    (BASE/'results'/'paper'/'frozen_data'/'offline_verification_results.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
