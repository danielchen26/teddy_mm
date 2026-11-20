#!/usr/bin/env python3
"""After the replication run: where is what Z leaves of MHC-I in co-culture? The same events read from RNA (the grid of
`anm_diag scope` on the SPEC mel_co_MHC1_rna of build_rna.py). This reads the head as trained only; no repaired output.

  python after.py --observe OBSERVE        (writes results/after_replication.json)
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import round as R  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--observe", type=Path, required=True)
a = ap.parse_args()
out = {"note": "computed after the replication run; the head as trained only"}
for g, r in (("mel_co", "MHC1"),):
    by, dz = R.observed(a.observe, g, r)
    rr = {x["event"]["unit"]: x for x in json.loads((a.observe / f"{g}_{r}_rna" / "scope.json").read_text())["results"]}
    dr = np.load(a.observe / f"{g}_{r}_rna" / "scope_draws.npz")
    view = {}
    for u in by:
        S = rr[u]["targets"]["shift"]["stages"]["RNA512"]
        Z = by[u]["targets"]["shift"]["stages"]["Z"]
        d = np.abs(dz[f"{u}|shift|Z|own|{Z['own']['class']}"].astype(np.float64)) - np.abs(dr[f"{u}|shift|RNA512|own|{S['own']['class']}"].astype(np.float64))
        iv = R.interval(d[1:])
        view[u] = {"at_Z": {k: {q: Z[k][q] for q in ("u", "ci", "verdict")} for k in ("units", "own")},
                   "from_RNA": {k: {q: S[k][q] for q in ("u", "ci", "verdict")} for k in ("units", "own")},
                   "added_from_RNA_to_Z_own_scope": {"added": float(d[0]), "ci": iv, "reading": "loss" if iv[0] > R.STEP else ("no added loss" if iv[1] < R.STEP else "unresolved")}}
    out[f"{g}_{r}"] = view
    p = view["pool"]
    print(g, r, "| at Z: units", round(p["at_Z"]["units"]["u"], 2), p["at_Z"]["units"]["verdict"], "own", round(p["at_Z"]["own"]["u"], 2), p["at_Z"]["own"]["verdict"],
          "| from RNA: units", round(p["from_RNA"]["units"]["u"], 2), p["from_RNA"]["units"]["verdict"], "own", round(p["from_RNA"]["own"]["u"], 2), p["from_RNA"]["own"]["verdict"],
          "| RNA -> Z", round(p["added_from_RNA_to_Z_own_scope"]["added"], 3), [round(x, 3) for x in p["added_from_RNA_to_Z_own_scope"]["ci"]], p["added_from_RNA_to_Z_own_scope"]["reading"])
(R.RES / "after_replication.json").write_text(json.dumps(out, indent=1))
