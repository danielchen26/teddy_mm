"""Run one declaration of this folder on exported features:  python run_one.py <run name>
The declaration is <run name>/declaration.yaml; the record is written next to it."""
import sys
import time

from _paths import DATASET_DIR, FEATURES, SPLIT
from anm_diag.runner import run
from anm_diag.spec import load_declaration
from precomputed_adapter import PrecomputedAdapter

name = sys.argv[1]
decl = load_declaration(str(SPLIT / name / "declaration.yaml"))
adapter = PrecomputedAdapter(str(FEATURES / DATASET_DIR[name.split("_")[0]]))
t0 = time.time()
rec = run(decl, adapter, str(SPLIT / name))
print("DONE", name, f"{time.time() - t0:.0f}s", {rid: {s: round(v["u"], 3) for s, v in r["states"].items() if v.get("role") == "chain"} for rid, r in rec["readouts"].items()})
