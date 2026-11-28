#!/usr/bin/env python3
"""A made-up world in the formats the real scripts read, to find faults of the code and to show that C2 can be
answered either way, before anything real is opened. No real value is read.

  python rehearsal.py --out DIR [--seed 1] [--world control_is_enough|rows_teach]

The states are made up and the head is the real one. Every event and every drawn unit has its state pushed along one
drawn direction, so it has a real change and the state carries all of it. The readout of a cell is how far the
activation the frozen part of the head gives it lies along the direction those pushes take there, standardized by
the control cells, plus noise: a world where a row of the last layer can read it, so that what the judge gives can
be told from what it should give. The two worlds differ in whether the control
cells can pin the row down:

  control_is_enough   3,000 control cells and little noise: a fit on them alone finds the row, so the loop's rows
                      can add nothing. C2 should come out refuted or unresolved.
  rows_teach          200 control cells (fewer than the 512 numbers of a row) and four times the noise: the fit on
                      them is not determined, and the rows are the only usable information. C2 should hold.

The second is not the regime of the real data, where the teach control cells are thousands; it is there to show that
the difference C2 asks about is found when it is there. In both, the head as shipped knows nothing of the made-up
readout, so C1's H0 should be near zero.

-> DIR/{cells.tsv, units.json, states.npz, pack.npz, planted.json}
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _hc as K  # noqa: E402

PANEL = ["CD274", "HLA_A", "HLA_E", "CD117", "CD119", "CD140a", "CD140b", "CD172a", "CD184", "CD202b", "CD29", "CD309", "CD44", "CD47", "CD49f", "CD58", "CD59", "CD61", "CD9", "CD279",
         "Rat_IgG2a", "Mouse_IgG1", "Mouse_IgG2a", "Mouse_IgG2b"]
COLUMN = {"PDL1": "CD274", "MHC1": "HLA_A", "HLAE": "HLA_E"}
PUSH = 25.0                                                              # the size of the push along q, in the units of the state


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--world", choices=("control_is_enough", "rows_teach"), default="control_is_enough")
    ap.add_argument("--control", type=int)
    ap.add_argument("--noise", type=float)
    ap.add_argument("--recorded", type=int, default=20)
    ap.add_argument("--drawn", type=int, default=10)
    ap.add_argument("--per", type=int, default=600)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    n_ctl = a.control if a.control else (1500 if a.world == "control_is_enough" else 200)
    noise = a.noise if a.noise else (0.5 if a.world == "control_is_enough" else 2.0)
    rng = np.random.default_rng(a.seed)
    q = rng.normal(size=512)
    q /= np.linalg.norm(q)
    rows, ids, Z, kind, push = [], [], [], [], {}
    n = 0
    for c in K.CONDITIONS:
        for who, many, role in [("control", n_ctl, "control")] + [(f"REC{i:02d}", a.per, "recorded") for i in range(a.recorded)] + [(f"DRW{i:02d}", a.per, "drawn") for i in range(a.drawn)]:
            name = [f"C{n + i}" for i in range(many)]
            n += many
            z = rng.normal(size=(many, 512))
            if role != "control":
                push[(c, who)] = float(rng.normal(0, PUSH))
                z = z + push[(c, who)] * q                               # the state itself is pushed, as a real event pushes it
            rows += [("frangieh2021", c, who, s, role) for s in name]
            ids += name
            Z.append(z)
            kind += [role] * many
    Z = np.vstack(Z)
    kind, unit = np.array(kind), np.array([r[2] for r in rows])
    ctl = kind == "control"
    A, _ = K.head(Z, ids)
    M = np.array([A[unit == g].mean(0) for g in sorted(set(unit[~ctl].tolist()))])
    v = np.linalg.svd(M - M.mean(0), full_matrices=False)[2][0]          # the direction the pushes take in the activation
    along = A @ v                                                        # the readout: a row of the last layer can read it
    m, s = float(along[ctl].mean()), float(along[ctl].std())
    y = (along - m) / s + rng.normal(size=len(ids)) * noise
    log = rng.normal(0, 0.2, (len(ids), 24))
    for r, col in COLUMN.items():
        log[:, PANEL.index(col)] = y * (1.0 if r == "PDL1" else (0.6 if r == "MHC1" else 0.3)) + rng.normal(size=len(ids)) * 0.2
    with open(a.out / "cells.tsv", "w") as f:
        f.write("study\tcondition\ttarget\tcell\trole\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    drawn = sorted({f"{g}@frangieh2021 | {c}" for (c, g) in push if g.startswith("DRW")})
    (a.out / "units.json").write_text(json.dumps({"seed": a.seed, "drawn": {"n": len(drawn), "units": drawn}, "kept_closed": {"n": 0, "units": []}}, indent=1) + "\n")
    np.savez(a.out / "states.npz", ids=np.array(ids), z_last=Z.astype(np.float32), sets=np.array(["A|0"]))
    np.savez(a.out / "pack.npz", mel_panel=np.array(PANEL), mel_ids=np.array(ids), mel_log=log, thp_panel=np.array(["CD86", "PDL1", "PDL2", "CD366"]), thp_ids=np.array([], dtype="<U8"), thp_log=np.zeros((0, 4)))
    moved = float(np.mean([abs((along[unit == g] / 1).mean() - along[ctl].mean()) / s for g in sorted(set(unit[~ctl].tolist()))]))
    (a.out / "planted.json").write_text(json.dumps({"world": a.world, "noise": noise, "seed": a.seed, "cells": len(ids), "control_cells": n_ctl, "drawn_units": len(drawn),
                                                    "C2_should_come_out": "holds" if a.world == "rows_teach" else "refuted or unresolved",
                                                    "the_readout_among_control_cells": {"spread": 1.0, "against_noise": noise},
                                                    "a_push_moves_the_readout_by_spreads": moved}, indent=1) + "\n")
    print("cells", len(ids), "| drawn units", len(drawn), "|", a.world, "| control cells", n_ctl, "| noise", noise, "| a push moves the readout by", round(moved, 2), "spreads")


if __name__ == "__main__":
    main()
