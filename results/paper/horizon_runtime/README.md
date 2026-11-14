# Historical-horizon benchmark protocol

These benchmarks compare final-horizon and all-horizon sensitivity restrictions on identical synthetic data. They measure the additional cost of earlier-horizon constraints without changing the leaf-probability LP objectives.

Run from the repository root (requires NumPy, pandas, SciPy, and the original paper notebooks):

```bash
python scripts/benchmark_horizon_primary.py
python scripts/benchmark_horizon_scaling.py
python scripts/verify_horizon_runtime.py
```

To preserve the supplied CSV run, set the environment variable `ISHLR_OUTPUT_DIR` to a separate output directory. The scripts otherwise overwrite the saved runtime tables in this directory.

## Primary 16-leaf fixture

The primary experiment uses `H=3`, `Gamma=4`, and seed 20261004 from notebook `03_global_compatibility_gate.ipynb`. Each specification solves all 32 leaf-endpoint LPs. Final-horizon and all-horizon formulations have 63 and 135 affine inequality rows, respectively. Seven timed repetitions follow one warm-up, with alternating execution order. The archived median total times are 0.0415 s and 0.0483 s, and mean widths are 0.11484075 and 0.10966039.

## Separate scaling fixtures

The 16, 32, 64, 128, and 256-leaf scaling experiment uses one identical linear objective for each paired formulation. The 16-leaf scaling fixture is **different** from the primary fixture: it has 80/200 rows rather than 63/135. Reported times are seven-run medians following warm-up and exclude construction of LP coefficient matrices.

In the scaling fixtures, the tested objective does not tighten when historical constraints are added. This experiment measures the computational overhead of the additional inequalities, not a general identification gain or a solver-independent complexity bound. Exact wall-clock times depend on hardware, solver versions, process load, and numerical tolerances.
