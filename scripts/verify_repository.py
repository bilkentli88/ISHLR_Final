#!/usr/bin/env python3
"""Run deterministic, no-network reproduction checks for published numerical results."""
from pathlib import Path
import subprocess, sys
ROOT=Path(__file__).resolve().parents[1]
checks=['verify_packaged_results.py','verify_horizon_runtime.py',
        'verify_reported_derived_results.py',
        'reproduce_frozen_inaturalist.py','verify_reference_frequencies.py',
        'reproduce_sensitivity_path.py']
for script in checks:
    print(f'\n=== {script} ===',flush=True)
    subprocess.run([sys.executable,str(ROOT/'scripts'/script)],cwd=ROOT,check=True)
print('\nAll deterministic offline checks passed.',flush=True)
