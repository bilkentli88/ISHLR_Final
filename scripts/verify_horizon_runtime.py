"""Check archived paired benchmark CSV structure and deterministic LP outputs.

Does not require time measurements to match another machine; does not rerun LPs.
To recompute use benchmark_horizon_primary.py and benchmark_horizon_scaling.py.
"""
from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]
DIR = ROOT / "results/paper/horizon_runtime"

with (DIR / "primary_runtime.csv").open(newline="") as f:
    primary = {r["formulation"]: r for r in csv.DictReader(f)}
assert set(primary) == {"final", "all"}
for mode, count, width in [("final", 63, 0.11484074798142302),
                           ("all", 135, 0.10966039454569984)]:
    r = primary[mode]
    assert int(r["leaves"]) == 16 and int(r["endpoint_LPs"]) == 32
    assert int(r["inequalities"]) == count
    assert abs(float(r["mean_leaf_interval_width"]) - width) < 1e-10
    assert int(r["timing_repeats"]) == 7
    assert float(r["median_32_LP_time_seconds"]) > 0

with (DIR / "scaling_runtime.csv").open(newline="") as f:
    scaling = list(csv.DictReader(f))
expected = [(16,80,200), (32,272,744), (64,576,1584),
            (128,2112,6064), (256,4352,12512)]
assert len(scaling) == len(expected)
for row, (n, final_rows, all_rows) in zip(scaling, expected):
    assert (int(row["leaves"]), int(row["final_inequalities"]),
            int(row["all_inequalities"])) == (n, final_rows, all_rows)
    assert float(row["final_median_solve_seconds"]) > 0
    assert float(row["all_median_solve_seconds"]) > 0
    assert abs(float(row["final_upper_objective"]) -
               float(row["all_upper_objective"])) < 1e-8

print("Paired horizon benchmark CSV checks passed (excluding variable timings).")
