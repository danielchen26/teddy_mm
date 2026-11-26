#!/usr/bin/env python3
"""The units of batch 3, from the census of 22 November and the roles of batch 2. results/roles.json. Record only.

  python roles.py

Seven units never read, each named below with the reason. The script stops if a unit is final, was read (development,
batch 1, batch 2), or is not a unit of the census.

The teachers of a type: for the readouts never read before, the recorded events of the interferon pathway of a stratum
are split by a seeded rule. In alphabetical order they are permuted with numpy's default generator seeded with
[20261126, k] (k = 0 for the co-culture, 1 for THP-1); the first `TEACH` are teachers, whose values are opened with the
control cells before any forecast, and the rest are forecast.
"""
import json
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
CENSUS = RUN.parent / "outer_step_20261122" / "results" / "census.json"
BATCH2 = RUN.parent / "batch2_forecast_20261125" / "results" / "roles.json"
SEED = 20261126
CO, CT, TH = "frangieh2021 | Co-culture", "frangieh2021 | Control", "papalexi2021 | all"
UNITS = {
    CO: {"CD47": "own gene of a readout never read (CD47)", "CD58": "own gene of a readout never read (CD58)", "CD59": "own gene of a readout never read (CD59)",
         "ITGA3": "no declared edge; the screen of 13 November lists it on CD49f: whether the law flags it under interferon, and whether that holds here, is asked by rule",
         "NONO": "an edge to PD-L1 proposed by the law of the library under interferon (run 43) that had no row elsewhere"},
    CT: {"CD47": "own gene of a readout never read (CD47)", "NONO": "an edge to PD-L1 proposed by the law of the library under interferon (run 43) that had no row elsewhere"},
}
PATHWAY = {CO: ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"], TH: ["IFNGR1", "IFNGR2", "JAK2", "STAT1"]}
TEACH = {CO: 3, TH: 2}


def main():
    c = json.load(open(CENSUS))
    roles, b1, b2 = c["roles"], set(c["batch_1"]["units"]), set(json.load(open(BATCH2))["units"])
    out = {"sources": ["runs/outer_step_20261122/results/census.json", "runs/batch2_forecast_20261125/results/roles.json"], "seed": SEED, "units": {}, "teachers": {}}
    for s, e in UNITS.items():
        for g, why in e.items():
            u = f"{g}@{s}"
            assert u in roles["unread_not_final"] and u not in roles["final"] and u not in roles["development"] and u not in b1 and u not in b2, u
            out["units"][u] = why
    for k, s in enumerate((CO, TH)):
        ev = sorted(PATHWAY[s])
        assert all(f"{g}@{s}" in roles["development"] for g in ev), s
        order = [ev[i] for i in np.random.default_rng([SEED, k]).permutation(len(ev))]
        out["teachers"][s] = {"teachers": sorted(order[:TEACH[s]]), "forecast": sorted(order[TEACH[s]:])}
    free = [u for u in roles["unread_not_final"] if u not in b1 and u not in b2]
    out["counts"] = {"units": len(out["units"]), "final_units_in_the_batch": 0, "unread_and_not_final_left_after_the_batch": len(free) - len(out["units"])}
    (RUN / "results" / "roles.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    for u, why in out["units"].items():
        print(u, "|", why)
    print(out["teachers"], out["counts"])


if __name__ == "__main__":
    main()
