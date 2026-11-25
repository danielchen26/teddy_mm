#!/usr/bin/env python3
"""The units of batch 2, from the census of 22 November. results/roles.json. Record only.

  python roles.py

Targeted units are named below with the reason. Calibration units: per melanoma condition, the units that are unread and
not final, not of batch 1 and not targeted, in alphabetical order, permuted with numpy's default generator seeded with
[20261125, k] (k = 0 for the co-culture, 1 for the condition with no signal); the first four. THP-1: every unit that is
unread and not final. The script stops if a unit is final, was read, or is not a unit of the census.
"""
import json
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
CENSUS = RUN.parent / "outer_step_20261122" / "results" / "census.json"
SEED, DRAWN = 20261125, 4
CO, CT, TH = "frangieh2021 | Co-culture", "frangieh2021 | Control", "papalexi2021 | all"
TARGETED = {
    CO: {"CD274": "own gene of PD-L1", "HLA-B": "own gene of MHC-I", "HLA-C": "own gene of MHC-I", "HLA-E": "own gene of HLA-E",
         "STAT3": "flagged on MHC-I by the law of the library under interferon (run 41)", "PFN1": "named in the plan of the batch; no edge, no flag"},
    CT: {"CD274": "own gene of PD-L1", "HLA-A": "own gene of MHC-I", "HLA-B": "own gene of MHC-I", "HLA-E": "own gene of HLA-E",
         "SSR2": "flagged on MHC-I by the law of the library under interferon (run 41)", "CD44": "flagged on MHC-I by the law of the library under interferon (run 41)",
         "IFNGR1": "the pathway with no signal given", "IFNGR2": "the pathway with no signal given", "JAK2": "the pathway with no signal given"},
}


def main():
    c = json.load(open(CENSUS))
    roles, b1 = c["roles"], set(c["batch_1"]["units"])
    free = [u for u in roles["unread_not_final"] if u not in b1]
    out = {"source": "runs/outer_step_20261122/results/census.json", "seed": SEED, "units": {}, "strata": {}}
    for k, s in enumerate((CO, CT)):
        pool = sorted(u.split("@")[0] for u in free if u.endswith("@" + s) and u.split("@")[0] not in TARGETED[s])
        drawn = [pool[i] for i in np.random.default_rng([SEED, k]).permutation(len(pool))[:DRAWN]]
        out["strata"][s] = {"targeted": TARGETED[s], "calibration": sorted(drawn), "pool_of_the_draw": len(pool)}
    out["strata"][TH] = {"targeted": {}, "calibration": sorted(u.split("@")[0] for u in free if u.endswith("@" + TH)), "pool_of_the_draw": None}
    for s, e in out["strata"].items():
        for g in list(e["targeted"]) + e["calibration"]:
            u = f"{g}@{s}"
            assert u in roles["unread_not_final"] and u not in roles["final"] and u not in roles["development"] and u not in b1, u
            out["units"][u] = "targeted" if g in e["targeted"] else "calibration"
    out["counts"] = {"units": len(out["units"]), "targeted": sum(v == "targeted" for v in out["units"].values()), "calibration": sum(v == "calibration" for v in out["units"].values()),
                     "final_units_in_the_batch": 0, "unread_and_not_final_left_after_the_batch": len(free) - len(out["units"])}
    (RUN / "results" / "roles.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    for s, e in out["strata"].items():
        print(s, "| targeted", list(e["targeted"]), "| calibration", e["calibration"])
    print(out["counts"])


if __name__ == "__main__":
    main()
