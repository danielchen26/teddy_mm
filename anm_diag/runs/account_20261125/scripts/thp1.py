#!/usr/bin/env python3
"""Section 7.2 (explored, not counted): THP-1, the five events on PD-L1. results/thp1[_EVENTS].json.

  python thp1.py --work WORK --states STATES --trained TRAINED [--only IRF1]

For every event: the account at Z, V, T, R, P (programs of the control cells) and ZV; for IRF1 also at the five last
layers trained on 24 November (the stage of a half is the state under the layer trained on that half, as in run 40).
"""
import argparse
import json
from pathlib import Path

import numpy as np

from _read import RUN, account, cells, typed, unit

PATHWAY = ["IFNGR1", "IFNGR2", "JAK2", "STAT1"]
ARMS = ("rows4", "rows5", "mixed", "cells5", "prop")
PAIRS = (("Z", "T"), ("R", "T"), ("Z", "V"), ("T", "V"), ("Z", "P"), ("Z", "ZV")) + tuple(("Z", f"trained:{a}") for a in ARMS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--states", type=Path, required=True)
    ap.add_argument("--trained", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=RUN / "results")
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    only = [t for t in a.only.split(",") if t] or ["IRF1"] + PATHWAY
    c = cells(a.work, "thp1", "PDL1")
    st, _worst = typed(a.states, "thp1", c)
    trained = {}
    for arm in ARMS:
        F = [np.load(a.trained / f"{arm}_{f}.npz") for f in (0, 1)]
        assert all(np.array_equal(d["ids"].astype(str), c["ids"]) for d in F), "the trained states are of other cells"
        trained[f"trained:{arm}"] = tuple(unit(d["ZT"]).astype(np.float64) for d in F)
    res = {"note": "THP-1, development cells; explored before and during the design: reported, not counted", "events": {}}
    for e in only:
        stages = {"Z": c["Z"], "V": c["V"], "T": st[f"T:A-{e}"] if f"T:A-{e}" in st else st["T"], "R": st["R"], "P": st["P"], "ZV": c["ZV"], **(trained if e == "IRF1" else {})}
        res["events"][e] = account(stages, c, e, [h for h in PATHWAY if h != e], pairs=PAIRS, say=lambda s_: print(s_, flush=True))
        (a.results / ("thp1.json" if not a.only else f"thp1_{'_'.join(only)}.json")).write_text(json.dumps(res, indent=1, default=float) + "\n")


if __name__ == "__main__":
    main()
