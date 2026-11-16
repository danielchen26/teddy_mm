#!/usr/bin/env python3
"""Run the SPEC files of a work folder through `anm_diag scope`, one after another.

  python run.py --work WORK --out ../results [--only ko_melanoma_pooled ko_thp1_pooled ...]
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from anm_diag.scope import main as scope_main  # noqa: E402

p = argparse.ArgumentParser()
p.add_argument("--work", type=Path, required=True)
p.add_argument("--out", type=Path, required=True)
p.add_argument("--only", nargs="*")
a = p.parse_args()
for spec in sorted(a.work.glob("*.yaml")):
    if a.only and spec.stem not in a.only:
        continue
    t0 = time.time()
    print(f"== {spec.stem}", flush=True)
    scope_main([str(spec), "--out", str(a.out / spec.stem)])          # the same entry as `python -m anm_diag scope`
    (a.out / spec.stem / "SPEC.yaml").write_text(spec.read_text())
    print(f"   done in {time.time() - t0:.0f} s", flush=True)
