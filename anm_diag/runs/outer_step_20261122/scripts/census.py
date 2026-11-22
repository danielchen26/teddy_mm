#!/usr/bin/env python3
"""The census of the pool from the record, the units the system has been read on, the final units, and a first batch.

  python census.py --record record_census.json --screen ../../loop_reading_20261112/results/data_census.json --out ../results

Environment: ANM_TEST_CELLS (library.npz, coculture.npz) and ANM_LOOP_WORK (thp1.npz): only the cell names and targets of
the cells that have states are read from them. No measurement and no state is read here.

record_census.json is written by record_census.py from the two protein files (cell names, guides, design fields only).

A unit is a target gene in a condition of a study, with at least 100 cells that carry guides for that gene only. A unit is
`read` if a state of the system exists for its cells (then some run has read the system on it), otherwise `unread`.
Final units: per stratum (study x condition), the unread units in alphabetical order are permuted with seed 20261122 and
the first third (rounded) is final. An unread control cell is final if sha256("final|<cell name>") is divisible by 3; the
rule is applied when states are computed, the count here is the expected one.

The batch is a proposal for which units get states next (the owner computes them). Chosen units: unread, not final, with
a clear response of one of the claim's readouts in the screen of 13 November, not the readout's own gene. Then half as
many again drawn at random from the other unread units that are not final, so that a third of the batch is a
calibration sample. The screen chooses what to compute; it admits nothing (admission is on a selection third).
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np

SEED, MIN_CELLS = 20261122, 100
READOUTS = {"CD274": "PD-L1", "PDL1": "PD-L1", "HLA_A": "MHC-I", "HLA_E": "HLA-E"}
OWN = {"PD-L1": {"CD274"}, "MHC-I": {"HLA-A", "HLA-B", "HLA-C"}, "HLA-E": {"HLA-E"}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--record", type=Path, required=True)
    ap.add_argument("--screen", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    rec, scr = json.loads(a.record.read_text()), json.loads(a.screen.read_text())
    cells = Path(os.environ["ANM_TEST_CELLS"])
    have = {}
    for key, f, tk in ((("frangieh2021", "IFNγ"), cells / "library.npz", "target"), (("frangieh2021", "Co-culture"), cells / "coculture.npz", "target"),
                       (("papalexi2021", "all"), Path(os.environ["ANM_LOOP_WORK"]) / "thp1.npz", "unit")):
        d = np.load(f)
        t, c = np.unique(d[tk].astype(str), return_counts=True)
        have[key] = dict(zip(t.tolist(), c.tolist()))
    strata = [("frangieh2021", s, v) for s, v in rec["frangieh2021"]["by_condition"].items()] + [("papalexi2021", "all", rec["papalexi2021"]["all"])]
    rng, rng2 = np.random.default_rng(SEED), np.random.default_rng(SEED + 1)      # the final draw (strata in sorted order), and the calibration sample
    out = {"rule": {"unit": f"a target gene in a condition of a study with at least {MIN_CELLS} cells that carry guides for that gene only",
                    "read": "a state of the system exists for cells of the unit", "final": "per stratum, the unread units in alphabetical order permuted with seed 20261122; the first third, rounded",
                    "final_control_cells": "an unread control cell is final if sha256('final|<cell name>') is divisible by 3"}, "strata": {}, "roles": {"final": [], "development": [], "unread_not_final": []}}
    clear = {}
    for study, conds in (("frangieh2021", scr["frangieh2021"]["conditions"]), ("papalexi2021", scr["papalexi2021"]["conditions"])):
        for cnd, v in conds.items():
            for x in v["clear_responses"]:
                r = READOUTS.get(x["antibody"])
                if r and x["target"] not in OWN[r]:
                    clear.setdefault((study, cnd), {}).setdefault(x["target"], []).append({"readout": r, "shift_sd": x["shift_sd"], "over_half_width": x["over_half_width"]})
    batch = []
    for study, cnd, v in sorted(strata, key=lambda x: (x[0], x[1])):
        h = have.get((study, cnd), {})
        units = sorted(u for u, n in v["units"].items() if u != "control" and n >= MIN_CELLS)
        read = [u for u in units if u in h]
        unread = [u for u in units if u not in h]
        perm = [unread[i] for i in rng.permutation(len(unread))]
        nf = int(round(len(unread) / 3))
        final, rest = sorted(perm[:nf]), sorted(perm[nf:])
        nc, nch = v["units"].get("control", 0), h.get("control", 0)
        name = f"{study} | {cnd}"
        cl = clear.get((study, cnd), {})
        chosen = [u for u in rest if u in cl]
        others = [u for u in rest if u not in cl]
        calib = sorted(others[i] for i in rng2.permutation(len(others))[: int(np.ceil(len(chosen) / 2))])
        out["strata"][name] = {"cells": v["cells"], "units": len(units), "read": len(read), "unread": len(unread), "final": len(final), "unread_not_final": len(rest),
                               "units_too_small": sorted(u for u, n in v["units"].items() if u != "control" and n < MIN_CELLS),
                               "control_cells": nc, "control_cells_with_states": nch, "control_cells_unread": nc - nch, "control_cells_final_expected": int(round((nc - nch) / 3)),
                               "screened_clear_and_final": sorted(u for u in final if u in cl), "screened_clear_and_read": sorted(u for u in read if u in cl),
                               "batch_chosen": {u: cl[u] for u in chosen}, "batch_calibration": calib,
                               "batch_cells": int(sum(v["units"][u] for u in chosen + calib))}
        out["roles"]["final"] += [f"{u}@{name}" for u in final]
        out["roles"]["development"] += [f"{u}@{name}" for u in read]
        out["roles"]["unread_not_final"] += [f"{u}@{name}" for u in rest]
        batch += [f"{u}@{name}" for u in chosen + calib]
        print(name, "| units", len(units), "read", len(read), "unread", len(unread), "final", len(final), "| controls", nc, "with states", nch,
              "| batch: chosen", chosen, "calibration", calib, "| clear but final:", out["strata"][name]["screened_clear_and_final"])
    out["batch_1"] = {"units": batch, "what_it_is_for": "units of the pool for the open items that need them (more events of a type, rows across conditions and cell types, the test of the split) "
                      "and a calibration sample; the states are to be computed by the owner, with a share of the unread control cells of each stratum that is not final"}
    out["totals"] = {k: len(v) for k, v in out["roles"].items()}
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "census.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    print(out["totals"], "| batch 1:", len(batch), "units")


if __name__ == "__main__":
    main()
