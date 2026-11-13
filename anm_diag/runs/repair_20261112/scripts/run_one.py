"""Run one audit of this folder:  python run_one.py <audit name>   (for example pooled_frangieh2021_joint_s0)"""
import sys
import time

from _paths import OUT_FEATURES, RUN
from anm_diag.runner import run
from anm_diag.spec import load_declaration
from precomputed_adapter import PrecomputedAdapter

name = sys.argv[1]
t0 = time.time()
rec = run(load_declaration(str(RUN / name / "declaration.yaml")), PrecomputedAdapter(str(OUT_FEATURES / name)), str(RUN / name))
print("DONE", name, f"{time.time() - t0:.0f}s", {rid: {s: round(v["u"], 3) for s, v in r["states"].items() if v.get("role") == "chain"} for rid, r in rec["readouts"].items()}, flush=True)
