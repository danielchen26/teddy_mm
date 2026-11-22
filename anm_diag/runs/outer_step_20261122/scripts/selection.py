#!/usr/bin/env python3
"""The selection split (protocol version 3, amendment 3d): events are admitted on a third of their cells, and those cells
enter no reading of u.

  python selection.py --work WORK --out ../results          (writes selection.json)

Environment: ANM_TEST_CELLS (library.npz); WORK is the work folder of runs/teg_loop_20261120 (mel_ifn.npz, mel_co.npz, thp1.npz).

A cell is in the selection third if the SHA-256 of "20261122|<cell name>" is divisible by 3. For a candidate event: its
cells in the third against the control cells in the third. A1: the change is resolved and at least three times the
half-width of its interval (Welch, cells as units). A2: the gene is not one of the readout's own.

1. Every knockout of the IFN-gamma library, for each readout. For MHC-I the events that pass are the events of the new
   target. For PD-L1 and HLA-E the count is reported (how many pass on a third against how many passed on all cells).
2. The events admitted in the round of 20 November keep their admission; how they fare on a third is reported.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np

OWN = {"PDL1": ["CD274"], "MHC1": ["HLA-A", "HLA-B", "HLA-C"], "HLAE": ["HLA-E"]}
ROUND = {"mel_ifn": ["PDL1"], "mel_co": ["PDL1", "MHC1"], "thp1": ["PDL1"]}
SEED, MIN_SNR = 20261122, 3.0
in_third = lambda cell: int(hashlib.sha256(f"{SEED}|{cell}".encode()).hexdigest(), 16) % 3 == 0


def screen(y, yc):
    d = y.mean() - yc.mean()
    half = 1.96 * np.sqrt(y.var(ddof=1) / len(y) + yc.var(ddof=1) / len(yc))
    return {"cells": int(len(y)), "change_in_control_sd": float(d / yc.std(ddof=1)), "snr": float(abs(d) / half), "A1": bool(abs(d) / half >= MIN_SNR)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    out = {"rule": f"a cell is in the selection third if sha256('{SEED}|<cell name>') is divisible by 3", "min_snr": MIN_SNR}
    lib = np.load(Path(os.environ["ANM_TEST_CELLS"]) / "library.npz")
    ids, tg = lib["ids"].astype(str), lib["target"].astype(str)
    third = np.array([in_third(i) for i in ids])
    C = (tg == "control") & third
    out["library"] = {"cells": int(len(ids)), "in_the_third": int(third.sum()), "control_cells_in_the_third": int(C.sum()), "candidates": int(len(set(tg)) - 1)}
    for r, own in OWN.items():
        y = lib[f"readout__{r}"].astype(float)
        rows = [{"event": g, **screen(y[(tg == g) & third], y[C]), "A2_own_gene": g in own} for g in sorted(set(tg) - {"control"})]
        rows.sort(key=lambda x: -x["snr"])
        out[r] = {"pass_A1": [x["event"] for x in rows if x["A1"]], "refused_A2": [x["event"] for x in rows if x["A1"] and x["A2_own_gene"]],
                  "admitted": [x["event"] for x in rows if x["A1"] and not x["A2_own_gene"]], "largest_not_passing": next((x for x in rows if not x["A1"]), None), "top": rows[:10]}
        print(r, "| pass A1 on a third:", out[r]["pass_A1"], "| admitted:", out[r]["admitted"], "| next:", out[r]["largest_not_passing"]["event"], round(out[r]["largest_not_passing"]["snr"], 2))
    out["events_of_the_round"] = {}
    for g, readouts in ROUND.items():
        K = np.load(a.work / f"{g}.npz")
        t3 = np.array([in_third(i) for i in K["ids"].astype(str)])
        u = K["unit"].astype(str)
        for r in readouts:
            y = K[f"y_{r}"].astype(float)
            rows = {s: screen(y[(u == s) & t3], y[(u == "control") & t3]) for s in sorted(set(u) - {"control"})}
            out["events_of_the_round"][f"{r} {g}"] = rows
            print(r, g, "| on a third:", {s: round(v["snr"], 1) for s, v in rows.items()})
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "selection.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
