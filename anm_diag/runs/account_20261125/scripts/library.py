#!/usr/bin/env python3
"""Section 7.3 (explored, not counted): the law of the rows of a library, with the edges. results/library.json.

  python library.py --library LIBRARY.npz [--results RESULTS]

LIBRARY: the cells of the melanoma library under IFN-gamma: Z (the state), target (the event of a cell), A (0 for a control
cell), readout__<name> (the protein readouts). An edge of class "own": the event's site is a gene of the measured protein.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from _read import AC, RUN, unit

OWN = {"PDL1": ["CD274"], "MHC1": ["HLA-A", "HLA-B", "HLA-C", "B2M"], "HLAE": ["HLA-E", "B2M"]}     # the first three lists are events.FORWARD_READOUT; B2M is the second chain of both MHC molecules


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=RUN / "results")
    a = ap.parse_args()
    d = np.load(a.library)
    un = np.where(d["A"] == 0, "control", d["target"].astype(str))
    Y = {k: d["readout__" + k].astype(np.float64) for k in OWN}
    edges = {(e, r): "own" for r, genes in OWN.items() for e in genes}
    res = {"note": "development units; explored: reported, not counted", "cells": int(len(un)), "control_cells": int((un == "control").sum()), "edges": [list(k) for k in edges if k[0] in set(un)]}
    for name, ed in (("with_edges", edges), ("without_edges", {})):
        law = AC.Law(unit(d["Z"]).astype(np.float64), Y, un, edges=ed).fit()
        rows = law.forecast()
        res[name] = {"summary": law.summary(), "calibration": law.calibration(rows), "candidates": law.flags(rows),
                     "pairs_with_a_change_or_an_edge": [r for r in rows if abs(r["change"]) / r["se"] >= AC.MIN_SNR or (r["event"], r["readout"]) in edges]}
        s, cal = res[name]["summary"], res[name]["calibration"]
        print(f"{name}: {s['events']} events; widening {s['widening_of_the_control_prior']}; spread of a direct part without an edge {s['spread_of_a_direct_part_without_an_edge']:.3f}; classes {s['classes']}")
        print(f"   inside 95%: {cal['inside_95']:.2f} of {cal['pairs']}; mean z2 {cal['mean_z2']:.2f}; error: {cal['error']}; pairs with a change: {cal['pairs_with_a_change']}; by class: {cal.get('by_class')}")
        for f in res[name]["candidates"]:
            print(f"   candidate: {f['event']:8s} -> {f['readout']:5s} direct {f['direct']:+.2f} ({f['errors']:+.1f} errors; shift {f['shift_times_noise']:.1f}x noise): {f['candidate']}")
    for r in res["with_edges"]["pairs_with_a_change_or_an_edge"]:
        print(f"   {r['event']:8s} -> {r['readout']:5s} change {r['change']:+.2f} | read {r['read']:+.2f} | others {r['from_the_other_events']:+.2f} | edge {r['by_its_edge']:+.2f} | forecast {r['forecast']:+.2f} +- {r['sd']:.2f} (z {r['z']:+.1f}) {r['edge'] or ''}")
    (a.results / "library.json").write_text(json.dumps(res, indent=1, default=float) + "\n")


if __name__ == "__main__":
    main()
