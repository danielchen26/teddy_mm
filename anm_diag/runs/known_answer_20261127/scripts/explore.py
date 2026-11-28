#!/usr/bin/env python3
"""Explored after the reading; it describes and no line rests on it. Why do the five remove about half of the signal's
effect on a group's mean? On PD-L1, cell by cell: how many cells of a knockout look as if no signal had been given.

  python explore.py --cells final_cells.tsv --library library.npz --packs PACK...

-> results/explored.json. A cell is below the midpoint if its PD-L1 (centred log-ratio) is below the middle of the
medians of the control cells with and without the signal. The share of a knockout's cells that look unstimulated is
its share below the midpoint, placed between that of the control cells under the signal (0) and with no signal (1).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _ka as K  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--library", type=Path, required=True)
    ap.add_argument("--packs", type=Path, nargs="+", required=True)
    a = ap.parse_args()
    C = K.cells(a.cells, a.library)
    P = K.B3.Packs(a.packs)
    val = {k: P.readouts(K.STUDY, C[k]["ids"], ["PDL1"])[0]["PDL1"] for k in ("ifn", "ctl")}
    c0, c1 = val["ctl"][C["ctl"]["unit"] == "control"], val["ifn"][C["ifn"]["unit"] == "control"]
    mid = (np.median(c0) + np.median(c1)) / 2
    lo, hi = float((c1 < mid).mean()), float((c0 < mid).mean())
    others = np.isin(C["ifn"]["unit"], [g for g in K.knockouts(C, "ifn") if g not in K.FIVE and g != "CD274"])
    out = {"made_after_the_reading": True, "readout": "PDL1", "medians_of_the_control_cells": {"no signal": float(np.median(c0)), "IFN-γ": float(np.median(c1))}, "midpoint": float(mid),
           "below_the_midpoint": {"control cells under IFN-γ": lo, "control cells with no signal": hi, "the other knockouts under IFN-γ (CD274 left out)": float((val["ifn"][others] < mid).mean())}, "the_five_under_IFN-γ": {}}
    for g in K.FIVE:
        v = val["ifn"][C["ifn"]["unit"] == g]
        b = float((v < mid).mean())
        out["the_five_under_IFN-γ"][g] = {"cells": int(len(v)), "below_the_midpoint": b, "share_of_cells_that_look_unstimulated": (b - lo) / (hi - lo)}
        print(g, len(v), round(b, 3), round((b - lo) / (hi - lo), 2))
    (K.RUN / "results" / "explored.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    print(out["below_the_midpoint"])


if __name__ == "__main__":
    main()
