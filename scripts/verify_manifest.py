#!/usr/bin/env python3
"""Verify the SHA-256 integrity manifest distributed with this repository."""
from pathlib import Path
import hashlib

ROOT=Path(__file__).resolve().parents[1]
MANIFEST=ROOT/'MANIFEST.sha256'
assert MANIFEST.is_file(),'Missing MANIFEST.sha256'
checks=0
for line in MANIFEST.read_text().splitlines():
    if not line.strip():continue
    digest, rel=line.split('  ',1)
    path=ROOT/rel
    if not path.is_file():raise FileNotFoundError(rel)
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    if actual!=digest:raise AssertionError(f'Checksum mismatch: {rel}')
    checks+=1
print(f'Integrity verified for {checks} repository files.')
