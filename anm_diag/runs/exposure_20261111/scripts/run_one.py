"""Run one declaration of this folder:  python run_one.py <arm>_<target>   (for example with_IRF1_IRF1)"""
import sys
import time

from _paths import ARMS, OUT_FEATURES, RUN
from anm_diag.runner import run
from anm_diag.spec import load_declaration
from precomputed_adapter import PrecomputedAdapter

name = sys.argv[1]
arm = next(a for a in ARMS if name.startswith(a + "_"))
decl = load_declaration(str(RUN / name / "declaration.yaml"))
t0 = time.time()
rec = run(decl, PrecomputedAdapter(str(OUT_FEATURES / arm)), str(RUN / name))
print("DONE", name, f"{time.time() - t0:.0f}s", {rid: {s: round(v["u"], 3) for s, v in r["states"].items() if v.get("role") == "chain"} for rid, r in rec["readouts"].items()})
