#!/usr/bin/env python3
"""Extended synthetic compact-versus-explicit identification validation.

Self-contained test deriving utility functions from 04_general_tree_theorem_gate.ipynb.
Valid monotone trajectories use exact integration of the shared U mechanism.
The grid varies tree depth, horizon count, leaf-class imbalance, and sensitivity.
Only synthetic data; no iNaturalist records or external requests.

Run: python scripts/benchmark_extended_synthetic.py
"""
from __future__ import annotations
from collections import defaultdict
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from scipy.optimize import linprog


def generate_tree(seed, max_depth=3, min_branch=2, max_branch=3):
    rng = np.random.default_rng(seed)
    children = {}
    counter = [0]

    def rec(node, depth):
        if depth >= max_depth:
            return
        if depth > 0 and rng.random() < 0.18:
            return
        b = int(rng.integers(min_branch, max_branch + 1))
        ch = []
        for _ in range(b):
            counter[0] += 1
            ch.append(f"N{counter[0]}")
        children[node] = ch
        for c in ch:
            rec(c, depth + 1)

    rec("ROOT", 0)
    return children

def prep_tree(children):
    nodes = {"ROOT"}
    parent = {}
    for u, chs in children.items():
        nodes.add(u)
        for c in chs:
            nodes.add(c)
            parent[c] = u

    leaves = sorted(nodes - set(children))

    def path_to_root(node):
        p = [node]
        while p[-1] != "ROOT":
            p.append(parent[p[-1]])
        return list(reversed(p))

    def descendant_leaves(node):
        if node not in children:
            return [node]
        out = []
        for c in children[node]:
            out.extend(descendant_leaves(c))
        return out

    subleaves = {n: descendant_leaves(n) for n in nodes}

    def depth(node):
        d = 0
        while node != "ROOT":
            node = parent[node]
            d += 1
        return d

    return parent, leaves, path_to_root, subleaves, nodes, depth

def assign_valid_edge_curves(children, parent, seed, H=3):
    rng = np.random.default_rng(seed + 1234)
    q = {}
    base = np.linspace(0.25, 0.75, H)

    def rec(u, parent_q=None):
        if u not in children:
            return
        chs = children[u]
        sibling_factor = rng.uniform(0.90, 1.10, len(chs))

        if u == "ROOT":
            for c, f in zip(chs, sibling_factor):
                arr = np.clip(base * f, 0.03, 0.90)
                arr = np.maximum.accumulate(arr)
                q[c] = arr
                rec(c, arr)
        else:
            for c, f in zip(chs, sibling_factor):
                depth_factor = rng.uniform(0.55, 0.78)
                arr = np.minimum(parent_q * depth_factor * f, parent_q * 0.95)
                arr = np.maximum.accumulate(arr)
                q[c] = arr
                rec(c, arr)

    rec("ROOT")
    return q

def max_sibling_odds_ratio(children, q, H=3):
    ans = 1.0
    for u, chs in children.items():
        for h in range(H):
            odds = [q[c][h] / (1.0 - q[c][h]) for c in chs]
            for a in odds:
                for b in odds:
                    ans = max(ans, a / b)
    return ans



def exact_trajectory_distribution(leaves, path_to_root, q, p, H=3):
    group_mass = defaultdict(float)

    for leaf in leaves:
        path = path_to_root(leaf)
        thresholds = sorted(set(
            [0.0, 1.0] +
            [float(q[node][h]) for node in path[1:] for h in range(H)]
        ))

        for a, b in zip(thresholds[:-1], thresholds[1:]):
            if b - a <= 1e-14:
                continue
            U = 0.5 * (a + b)
            states = []
            for h in range(H):
                state = "ROOT"
                for node in path[1:]:
                    if U <= q[node][h] + 1e-15:
                        state = node
                    else:
                        break
                states.append(state)
            group_mass[tuple(states)] += p[leaf] * (b - a)

    assert abs(sum(group_mass.values()) - 1.0) < 1e-10
    return group_mass

def observed_quantities(children, path_to_root, group_mass, H=3):
    r = defaultdict(lambda: defaultdict(lambda: [0.0] * H))
    final_mass = defaultdict(float)

    for g, mass in group_mass.items():
        final_mass[g[-1]] += mass
        for h, state in enumerate(g):
            state_path = set(path_to_root(state))
            for u, chs in children.items():
                for c in chs:
                    if c in state_path:
                        r[u][c][h] += mass

    return r, final_mass



def build_compact(children, leaves, subleaves, r, H=3, gamma=3.0, horizons=None):
    if horizons is None:
        horizons = list(range(H))

    idx = {l:i for i,l in enumerate(leaves)}
    n = len(leaves)

    def vec(node):
        v = np.zeros(n)
        for l in subleaves[node]:
            v[idx[l]] = 1.0
        return v

    A, b = [], []

    # Final-node compatibility: s_c >= M_c = r_c,H.
    for u, chs in children.items():
        for c in chs:
            A.append(-vec(c))
            b.append(-r[u][c][H-1])

    # Sibling sensitivity half-spaces.
    for u, chs in children.items():
        vv = {c:vec(c) for c in chs}
        for h in horizons:
            for j in chs:
                for k in chs:
                    if j == k:
                        continue
                    rj, rk = r[u][j][h], r[u][k][h]
                    A.append(rj * vv[k] - gamma * rk * vv[j])
                    b.append((1.0 - gamma) * rj * rk)

    return np.asarray(A), np.asarray(b), idx

def build_full_flow(leaves, path_to_root, group_mass, A_compact, b_compact, idx):
    groups = list(group_mass)
    pairs = []

    for gi, g in enumerate(groups):
        final = g[-1]
        for leaf in leaves:
            if final in path_to_root(leaf):
                pairs.append((gi, leaf))

    P = np.zeros((len(leaves), len(pairs)))
    Aeq = np.zeros((len(groups), len(pairs)))
    beq = np.array([group_mass[g] for g in groups])

    for zi, (gi, leaf) in enumerate(pairs):
        P[idx[leaf], zi] = 1.0
        Aeq[gi, zi] = 1.0

    return groups, pairs, P, Aeq, beq, A_compact @ P, b_compact



def run_grid(output: Path, n_seeds=3, n_objectives=4):
    rows = []
    solver_options = {'primal_feasibility_tolerance': 1e-10,
                      'dual_feasibility_tolerance': 1e-10,
                      'ipm_optimality_tolerance': 1e-10}
    output.mkdir(parents=True, exist_ok=True)
    for depth in (2, 3, 4):
        for H in (2, 5):
            for concentration in (0.35, 2.0):
                for k in range(n_seeds):
                    seed = 100 + k
                    children = generate_tree(seed, max_depth=depth)
                    parent, leaves, path_to_root, subleaves, nodes, depth_fn = prep_tree(children)
                    q = assign_valid_edge_curves(children, parent, seed, H=H)
                    max_or = float(max_sibling_odds_ratio(children, q, H))
                    assert max_or < 3, (depth, H, seed, max_or)
                    rg = np.random.default_rng(4100 + seed + int(100 * concentration))
                    masses = rg.dirichlet(np.repeat(concentration, len(leaves)))
                    p_true = dict(zip(leaves, masses))
                    group_mass = exact_trajectory_distribution(leaves, path_to_root, q, p_true, H)
                    r, final_mass = observed_quantities(children, path_to_root, group_mass, H)
                    eq = np.ones((1,len(leaves)))
                    b_eq = np.array([1.0])
                    for gamma in (3.0, 6.0):
                        A,b,idx = build_compact(children,leaves,subleaves,r,H,gamma)
                        groups,pairs,P,AeqF,beqF,AF,bF = build_full_flow(
                            leaves,path_to_root,group_mass,A,b,idx)
                        true_error = float(max(np.maximum(A @ masses-b,0).max(),abs(masses.sum()-1)))
                        assert true_error < 1e-10, ('true not feasible',true_error)
                        rg_obj = np.random.default_rng(20000+100*seed+int(10*gamma)+H)
                        max_gap = 0.0
                        for oi in range(n_objectives):
                            direction=rg_obj.normal(size=len(leaves))
                            # Include leaf marginal objectives among arbitrary directions.
                            if oi==0:
                                direction[:]=0.0
                                direction[0]=1.0
                            if oi==1:
                                direction[:]=0.0
                                direction[-1]=-1.0
                            rc=linprog(-direction,A_ub=A,b_ub=b,A_eq=eq,b_eq=b_eq,
                                       bounds=[(0,None)]*len(leaves),method='highs',options=solver_options)
                            rf=linprog(-(direction@P),A_ub=AF,b_ub=bF,A_eq=AeqF,b_eq=beqF,
                                       bounds=[(0,None)]*len(pairs),method='highs',options=solver_options)
                            if not(rc.success and rf.success):
                                raise RuntimeError(f'LP error depth={depth} H={H} a={concentration} seed={seed} gamma={gamma} objective={oi}: compact={rc.message}; flow={rf.message}')
                            gap=abs(float(rc.fun-rf.fun))
                            max_gap=max(max_gap,gap)
                        rows.append({
                            'max_depth':depth,'horizons':H,'dirichlet_concentration':concentration,
                            'seed':seed,'gamma':gamma,'leaves':len(leaves),
                            'internal_nodes':len(children),'trajectory_groups':len(groups),
                            'explicit_flow_variables':len(pairs),'compact_variables':len(leaves),
                            'compact_inequalities':len(b),'max_true_odds_ratio':max_or,
                            'max_true_feasibility_violation':true_error,
                            'support_objectives_checked':n_objectives,
                            'max_compact_flow_objective_gap':max_gap,
                        })
    df=pd.DataFrame(rows)
    df.to_csv(output/'extended_synthetic_stress.csv',index=False,float_format='%.15g')
    summary={
       'base_scenarios':int(len(df)//2),
       'distinct_tree_structures':int(len(df[['max_depth','seed']].drop_duplicates())),
       'configurations':int(len(df)),
       'support_comparisons':int(df.support_objectives_checked.sum()),
       'depths':sorted(df.max_depth.unique().astype(int).tolist()),
       'horizon_counts':sorted(df.horizons.unique().astype(int).tolist()),
       'dirichlet_concentrations':sorted(df.dirichlet_concentration.unique().tolist()),
       'gamma_values':sorted(df.gamma.unique().tolist()),
       'leaf_count_min':int(df.leaves.min()),'leaf_count_max':int(df.leaves.max()),
       'max_flow_variables':int(df.explicit_flow_variables.max()),
       'max_compact_inequalities':int(df.compact_inequalities.max()),
       'max_true_sibling_odds':float(df.max_true_odds_ratio.max()),
       'max_truth_violation':float(df.max_true_feasibility_violation.max()),
       'max_objective_difference':float(df.max_compact_flow_objective_gap.max()),
    }
    (output/'extended_synthetic_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    assert summary['max_objective_difference'] < 1e-7
    return df

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=Path('results/paper/extended_synthetic'))
    parser.add_argument('--seeds',type=int,default=3)
    parser.add_argument('--objectives',type=int,default=4)
    args=parser.parse_args()
    run_grid(args.out,args.seeds,args.objectives)
