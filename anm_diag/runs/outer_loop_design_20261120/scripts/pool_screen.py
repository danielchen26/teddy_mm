#!/usr/bin/env python3
"""Which knockouts of the library change a readout enough to be admitted: a census of events for the IFN-gamma condition.

  python pool_screen.py --work WORK --out ../results

Environment: ANM_TEST_CELLS (the folder written by runs/teg_repair_test_20261119/scripts/cells.py: library.npz).
WORK is the work folder of runs/teg_loop_20261120 (mel_ifn.npz: the control cells of the condition).

This reads responses of cells that earlier runs had read (runs/teg_repair_test_20261119 and the map of the whole
library). It is a planning number for the outer step, on all cells of an event, with an analytic interval (Welch). The
outer step itself screens on a selection split (OUTER_LOOP.md, Section 8.1).

A1: the change of the readout against the control cells is at least three times the half-width of its interval.
A2: the event does not act on the readout's own gene.
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np

OWN = {"PDL1": {"CD274"}, "MHC1": {"HLA-A", "HLA-B", "HLA-C"}, "HLAE": {"HLA-E"}}       # the readout's own genes (rule S3)
ROUND = {"IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"}                                   # the events of the round of 20 November
MIN_SNR = 3.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    lib = np.load(Path(os.environ["ANM_TEST_CELLS"]) / "library.npz", allow_pickle=True)
    ifn = np.load(a.work / "mel_ifn.npz", allow_pickle=True)
    C = ifn["A"] == 0
    tg = lib["target"].astype(str)
    out = {"condition": "melanoma, IFN-gamma", "control_cells": int(C.sum()), "targets": int(len(set(tg))), "min_snr": MIN_SNR,
           "note": "all cells of an event, analytic interval; cells read before; a planning number", "readouts": {}}
    for r in OWN:
        yc = ifn[f"y_{r}"][C].astype(float)
        rows = []
        for g in sorted(set(tg)):
            y = lib[f"readout__{r}"][tg == g].astype(float)
            d = y.mean() - yc.mean()
            half = 1.96 * np.sqrt(y.var(ddof=1) / len(y) + yc.var(ddof=1) / len(yc))
            rows.append({"event": g, "cells": int(len(y)), "change": float(d), "change_in_control_sd": float(d / yc.std(ddof=1)), "snr": float(abs(d) / half),
                         "A1": bool(abs(d) / half >= MIN_SNR), "A2_own_gene": g in OWN[r], "in_the_round": g in ROUND})
        rows.sort(key=lambda x: -x["snr"])
        ok = [x for x in rows if x["A1"] and not x["A2_own_gene"]]
        out["readouts"][r] = {"pass_A1": sum(x["A1"] for x in rows), "admissible": [x["event"] for x in ok], "admissible_not_in_the_round": [x["event"] for x in ok if not x["in_the_round"]],
                              "largest_snr_not_admissible_by_A1": max((x["snr"] for x in rows if not x["A1"]), default=None), "top": rows[:12]}
        print(r, "pass A1:", out["readouts"][r]["pass_A1"], "| admissible:", out["readouts"][r]["admissible"], "| new:", out["readouts"][r]["admissible_not_in_the_round"])
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "pool_screen.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
