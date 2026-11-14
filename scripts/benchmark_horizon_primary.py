"""Paired final-versus-all-horizon leaf-endpoint timing on the primary synthetic fixture.

Reuses the model-construction cells from the unmodified paper notebook.
Execute from any current directory: python scripts/benchmark_horizon_primary.py
The reported wall-clock times are environment-specific; numerical widths are not.
"""
import csv
import json
import os
import time
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.optimize import linprog

# Script is in scripts/; allow an explicit repository-root override.
repo_root = Path(os.environ.get('ISHLR_REPO', Path(__file__).resolve().parents[1]))
notebook = repo_root / 'notebooks/paper/03_global_compatibility_gate.ipynb'
if not notebook.exists():
    raise FileNotFoundError(f'Cannot find paper notebook at {notebook}. Set ISHLR_REPO.')
d = json.loads(notebook.read_text())
ns={'np':np,'defaultdict':defaultdict,'linprog':linprog,'SEED':20261004,'GAMMA':4.0,'H':3}
for cell_idx in [3,5,7,9]:
    exec(''.join(d['cells'][cell_idx]['source']),ns)
leaves=ns['leaves']; nL=ns['nL']
A_all,b_all=ns['A_all'],ns['b_all']
A_final,b_final=ns['build_compact_constraints'](horizons=[ns['H']-1])
Aeq,beq,bounds=ns['Aeq'],ns['beq'],ns['bounds']
assert len(b_all)==135 and len(b_final)==63


def endpoints(A,b):
    lower=[];upper=[]
    t0=time.perf_counter()
    for j in range(nL):
        c=np.zeros(nL);c[j]=1.0
        lo=linprog(c,A_ub=A,b_ub=b,A_eq=Aeq,b_eq=beq,bounds=bounds,method='highs')
        hi=linprog(-c,A_ub=A,b_ub=b,A_eq=Aeq,b_eq=beq,bounds=bounds,method='highs')
        if not lo.success or not hi.success:
            raise RuntimeError('infeasible at leaf '+str(j))
        lower.append(lo.fun);upper.append(-hi.fun)
    return time.perf_counter()-t0, np.asarray(upper)-np.asarray(lower)

for A,b in [(A_all,b_all),(A_final,b_final)]:
    endpoints(A,b)

timings={'final':[], 'all':[]}
width_by_mode={}
for i in range(7):
    # Alternate formulation order to avoid a systematic timing advantage.
    order = ['final','all'] if i%2==0 else ['all','final']
    for mode in order:
        A,b = (A_final,b_final) if mode=='final' else (A_all,b_all)
        seconds,widths=endpoints(A,b)
        timings[mode].append(seconds)
        width_by_mode[mode]=float(widths.mean())
results={mode: {'median_seconds':float(np.median(timings[mode])),
                'width':width_by_mode[mode],
                'individual_seconds':timings[mode]}
         for mode in ['final','all']}
assert nL == 16 and len(b_final) == 63 and len(b_all) == 135
assert abs(results['final']['width'] - 0.11484074798142302) < 1e-10
assert abs(results['all']['width'] - 0.10966039454569984) < 1e-10
print('main fixture:',nL,'leaves',len(b_final),len(b_all),'inequalities')
print('mean widths:',results['final']['width'],results['all']['width'])
print('temporal improvement:',1-results['all']['width']/results['final']['width'])
print('32 LP solves total elapsed median:',results['final']['median_seconds'],results['all']['median_seconds'])
print('ratio:',results['all']['median_seconds']/results['final']['median_seconds'])

outdir=Path(os.environ.get('ISHLR_OUTPUT_DIR', repo_root / 'results/paper/horizon_runtime'))
outdir.mkdir(parents=True, exist_ok=True)
out=outdir / 'primary_runtime.csv'
with out.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['leaves','horizons','gamma','formulation',
        'inequalities','endpoint_LPs','mean_leaf_interval_width','median_32_LP_time_seconds',
        'timing_repeats'])
    w.writeheader()
    for mode,b in [('final',b_final),('all',b_all)]:
        w.writerow(dict(leaves=nL,horizons=ns['H'],gamma=ns['GAMMA'],
            formulation=mode,inequalities=len(b),endpoint_LPs=2*nL,
            mean_leaf_interval_width=results[mode]['width'],
            median_32_LP_time_seconds=results[mode]['median_seconds'],
            timing_repeats=7))
print('Saved',out)
