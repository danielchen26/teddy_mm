#!/usr/bin/env python3
"""The units of run 52, drawn by lot before anything is opened. Names only; nothing is read.

  python units.py

Half of the units never read, per stratum, by lot with seed 20261128. The other half stays closed.
-> results/units.json
"""
import json
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
RUNS = RUN.parent
SEED = 20261128


def main():
    census = json.loads((RUNS / "outer_step_20261122" / "results" / "census.json").read_text())
    read = set(census["batch_1"]["units"])
    for run in ("batch2_forecast_20261125", "batch3_frozen_20261126", "batch4_frozen_20261127", "final_audit_20261127"):
        read |= set(json.loads((RUNS / run / "results" / "roles.json").read_text())["units"])
    left = sorted(set(census["roles"]["unread_not_final"]) - read)
    assert not set(census["roles"]["final"]) & set(left), "a unit of the final audit is in the pool"
    by = {}
    for u in left:
        by.setdefault(u.split("@")[1], []).append(u)
    rng = np.random.default_rng(SEED)
    drawn, kept = [], []
    for k, v in sorted(by.items()):
        order = rng.permutation(len(v))
        half = len(v) // 2
        drawn += [v[i] for i in sorted(order[:half])]
        kept += [v[i] for i in sorted(order[half:])]
    out = {"rule": {"pool": "every unit of the census never read by any run", "draw": f"per stratum, a random half (seed {SEED}); the rest stays closed",
                    "control_cells": "the control cells of each stratum that no run has read; a unit is read against them"},
           "seed": SEED, "pool": len(left), "drawn": {"n": len(drawn), "by_stratum": {k: sum(1 for u in drawn if u.endswith(k)) for k in sorted(by)}, "units": drawn},
           "kept_closed": {"n": len(kept), "units": kept}}
    (RUN / "results" / "units.json").write_text(json.dumps(out, indent=1) + "\n")
    print("pool", len(left), "| drawn", len(drawn), out["drawn"]["by_stratum"], "| kept closed", len(kept))


if __name__ == "__main__":
    main()
