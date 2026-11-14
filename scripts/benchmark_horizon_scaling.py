"""Paired final-versus-all-horizon timing on synthetic scaling fixtures.

Replicates the synthetic construction in the unmodified notebook's scaling audit,
using a common objective per horizon pair. After one warm-up, reports medians of
seven LP solver calls for each mode. Times are environment-specific.
Run: python scripts/benchmark_horizon_scaling.py
"""
import csv
import platform
import time
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import linprog

SEED = 20261004
GAMMA = 4.0
H = 3
REPEATS = 7
SIZES = [(4,4),(4,8),(8,8),(8,16),(16,16)]


def scaling_matrices(branches_n, leaves_per_branch, horizons, H=H, gamma=GAMMA):
    # Data-generating fixtures and sensitivity rows preserve the construction
    # from the repository's original scaling_problem function.
    rng = np.random.default_rng(SEED + branches_n*100 + leaves_per_branch)
    n = branches_n * leaves_per_branch
    p = rng.dirichlet(np.ones(n) * 2.0)
    branch_mass = p.reshape(branches_n, leaves_per_branch).sum(axis=1)
    base = np.array([0.28,0.52,0.76])
    qroot = np.zeros((branches_n,H))
    for i in range(branches_n):
        mult = 0.85 + 0.30 * (i / max(1, branches_n-1))
        qroot[i] = np.minimum(base * mult, 0.90)
    qleaf = np.zeros((branches_n,leaves_per_branch,H))
    for i in range(branches_n):
        for j in range(leaves_per_branch):
            factor = 0.52 + 0.26 * (j / max(1,leaves_per_branch-1))
            qleaf[i,j] = qroot[i] * factor

    A, b = [], []
    for i in range(branches_n):
        v = np.zeros(n)
        v[i*leaves_per_branch:(i+1)*leaves_per_branch] = -1
        A.append(v)
        b.append(-branch_mass[i] * qroot[i,-1])
        for j in range(leaves_per_branch):
            vv = np.zeros(n)
            vv[i*leaves_per_branch+j] = -1
            A.append(vv)
            b.append(-p[i*leaves_per_branch+j] * qleaf[i,j,-1])
    compatibility_rows = len(b)

    root_vec = []
    for i in range(branches_n):
        v = np.zeros(n)
        v[i*leaves_per_branch:(i+1)*leaves_per_branch] = 1
        root_vec.append(v)
    for h in horizons:
        rr = branch_mass*qroot[:,h]
        for i in range(branches_n):
            for j in range(branches_n):
                if i == j: continue
                A.append(rr[i]*root_vec[j] - gamma*rr[j]*root_vec[i])
                b.append((1-gamma)*rr[i]*rr[j])
    for i in range(branches_n):
        for h in horizons:
            rr = p[i*leaves_per_branch:(i+1)*leaves_per_branch] * qleaf[i,:,h]
            for j in range(leaves_per_branch):
                for k in range(leaves_per_branch):
                    if j == k: continue
                    v = np.zeros(n)
                    v[i*leaves_per_branch+k] += rr[j]
                    v[i*leaves_per_branch+j] -= gamma*rr[k]
                    A.append(v)
                    b.append((1-gamma)*rr[j]*rr[k])
    return np.asarray(A), np.asarray(b), compatibility_rows


def solve(matrices, n):
    A,b,_ = matrices
    obj = np.linspace(0.2, 2.0, n)
    t0 = time.perf_counter()
    res = linprog(-obj, A_ub=A, b_ub=b, A_eq=np.ones((1,n)),
                  b_eq=[1.0], bounds=[(0,None)]*n, method='highs')
    elapsed = time.perf_counter() - t0
    if not res.success:
        raise RuntimeError(f'LP failed (n={n}): {res.message}')
    return elapsed, -res.fun


def run():
    print('Python',platform.python_version(),'NumPy',np.__version__,
          'SciPy',scipy.__version__, 'Gamma',GAMMA, 'H',H,'repeats',REPEATS)
    rows=[]
    for b,k in SIZES:
        n=b*k
        builds={}
        m={}
        for name,horizons in [('final',[H-1]),('all',list(range(H)))]:
            t0=time.perf_counter()
            m[name]=scaling_matrices(b,k,horizons)
            builds[name]=time.perf_counter()-t0
        # Warm-up first, then interleave timings to limit order effects.
        for name in ['all','final']:
            solve(m[name],n)
        timings={'final':[],'all':[]}
        optima={'final':[],'all':[]}
        for i in range(REPEATS):
            order = ['final','all'] if i%2==0 else ['all','final']
            for name in order:
                sec, val=solve(m[name],n)
                timings[name].append(sec)
                optima[name].append(val)
        t_final=float(np.median(timings['final']))
        t_all=float(np.median(timings['all']))
        row=dict(branches=b,leaves_per_branch=k,leaves=n,
                 compatibility_rows=m['final'][2],
                 final_inequalities=len(m['final'][1]),
                 all_inequalities=len(m['all'][1]),
                 final_build_seconds=builds['final'],
                 all_build_seconds=builds['all'],
                 final_median_solve_seconds=t_final,
                 all_median_solve_seconds=t_all,
                 all_over_final_solve_ratio=t_all/t_final,
                 final_upper_objective=optima['final'][0],
                 all_upper_objective=optima['all'][0],
                 all_vs_final_upper_reduction=optima['final'][0]-optima['all'][0])
        assert row['all_inequalities'] >= row['final_inequalities']
        assert row['all_vs_final_upper_reduction'] >= -1e-9
        rows.append(row)
        print(f'{n:>3} leaves | rows {row["final_inequalities"]:>5}/{row["all_inequalities"]:>5} | '
              f'paired solve {t_final:.4f}/{t_all:.4f} s | ratio {t_all/t_final:.2f}x | '
              f'upper bound reduction {row["all_vs_final_upper_reduction"]:.8f}', flush=True)
    outdir=Path(__import__('os').environ.get('ISHLR_OUTPUT_DIR', Path(__file__).resolve().parents[1] / 'results/paper/horizon_runtime'))
    outdir.mkdir(parents=True, exist_ok=True)
    out=outdir / 'scaling_runtime.csv'
    with out.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0].keys())
        writer.writeheader();writer.writerows(rows)
    print('Saved',out)

if __name__=='__main__':
    run()
