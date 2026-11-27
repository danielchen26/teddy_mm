#!/usr/bin/env python3
"""A made-up file of the final cells, to find faults of read.py and score.py before the one reading.

  python rehearsal.py --cells final_cells.tsv --packs batch3_pack_A.npz batch3_pack_B.npz --out MADE_UP.npz

Nothing of a final cell is read but its line in the list. Every final cell gets the antibody columns of a control cell
of its condition read before, drawn at random, so no unit changes anything; then a few changes are planted.
"""
import argparse
from pathlib import Path

import numpy as np

PLANTED = {"mel": [("B2M", "HLA_A", -1.0), ("CD58", "CD58", -0.8), ("SSR2", "CD9", 0.4), ("ITGA3", "CD49f", -0.4), ("TIMP1", "CD44", -0.5)], "thp": [("CD86", "CD86", -1.0), ("STAT2", "PDL1", -0.4)]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--packs", type=Path, nargs="+", required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    rows = [l.rstrip("\n").split("\t") for l in a.cells.read_text().splitlines()[1:] if l.strip()]
    rng = np.random.default_rng(0)
    out = {}
    for k, study in (("mel", "frangieh2021"), ("thp", "papalexi2021")):
        log, panel = {}, None
        for f in a.packs:
            P = np.load(f)
            panel = [str(x) for x in P[f"{k}_panel"]]
            log.update(zip(P[f"{k}_ids"].astype(str).tolist(), P[f"{k}_log"]))
        pool = {}
        for s, cond, unit, cell, role in rows:
            if s == study and role == "control":
                pool.setdefault(cond, []).append(cell)
        keep, L = {x: [] for x in ("ids", "condition", "unit", "role")}, []
        for s, cond, unit, cell, role in rows:
            if s != study or role not in ("event", "reference"):
                continue
            assert cell not in log, "a final cell has values in the files exported before"
            v = np.array(log[pool[cond][rng.integers(len(pool[cond]))]], np.float32)
            for u, col, d in PLANTED[k]:
                if unit == u:
                    v[panel.index(col)] = max(v[panel.index(col)] + d, 0)
            L.append(v)
            for x, val in (("ids", cell), ("condition", cond), ("unit", unit), ("role", role)):
                keep[x].append(val)
        out.update({f"{k}_{x}": np.array(v) for x, v in keep.items() if not (k == "thp" and x == "condition")})
        out[f"{k}_log"], out[f"{k}_panel"] = np.array(L, np.float32), np.array(panel)
    np.savez(a.out, **out)
    print("made up:", out["mel_log"].shape, out["thp_log"].shape)


if __name__ == "__main__":
    main()
