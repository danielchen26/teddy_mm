#!/usr/bin/env python3
"""Section 7.4 (declared): three melanoma targets, by event. results/melanoma.json. Made once.

  python melanoma.py --work WORK --states STATES [--results RESULTS]

WORK: the cells of the strata with the state Z and the input view V (those of runs/input_beside_z_20261123).
STATES: the typed states of runs/typed_repair_20261125 (<stratum>_typed.npz); their gene sets were chosen on the fitted half.
For every event of a target: the account at Z, V, T (typed by the genes of the other four events) and R (genes at random).
"""
import argparse
import json
from pathlib import Path

import numpy as np

from _read import RUN, account, cells, typed

TARGETS = {"MHC1_coculture": ("mel_co", "MHC1"), "PDL1_coculture": ("mel_co", "PDL1"), "PDL1_ifn": ("mel_ifn", "PDL1")}
EVENTS = ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"]
PAIRS = (("Z", "T"), ("R", "T"), ("Z", "V"), ("T", "V"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--states", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=RUN / "results")
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    only = [t for t in a.only.split(",") if t] or list(TARGETS)
    out_file = a.results / ("melanoma.json" if not a.only else f"melanoma_{'_'.join(only)}.json")
    if out_file.exists():
        raise SystemExit(f"{out_file.name} exists: this reading was made once and is not made again")
    res = {"note": "development cells; read by run 40 with its repair, not before with `carried`; nothing trained", "targets": {}}
    for t in only:
        g, readout = TARGETS[t]
        c = cells(a.work, g, readout)
        st, worst = typed(a.states, g, c)
        res["targets"][t] = {"stratum": g, "readout": readout, "the_mean_of_the_typed_states_against_the_state": worst, "events": {}}
        for e in EVENTS:
            stages = {"Z": c["Z"], "V": c["V"], "T": st[f"T:A-{e}"] if f"T:A-{e}" in st else st["T"], "R": st["R"]}
            res["targets"][t]["events"][e] = account(stages, c, e, [h for h in EVENTS if h != e], pairs=PAIRS, say=lambda s_: print(f"{t:15s} {s_}", flush=True))
            out_file.with_suffix(".partial").write_text(json.dumps(res, indent=1, default=float) + "\n")
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    out_file.with_suffix(".partial").unlink()


if __name__ == "__main__":
    main()
