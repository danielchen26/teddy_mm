#!/usr/bin/env python3
"""A made-up file of the batch, to find faults of read.py and score.py before the one reading.

  python rehearsal.py --cells batch4_cells.tsv --packs batch3_pack_A.npz batch3_pack_B.npz --out MADE_UP.npz

Nothing of the batch is read but the list of its cells. Every cell of the batch gets the antibody columns of a control
cell of its condition, drawn at random from the files exported before, so no unit changes anything; then a few changes
are planted so that the edge tests, the accounts and the step by rule have something to read.
"""
import argparse
from pathlib import Path

import numpy as np

PLANTED = [("CXCR4", "CD184", -0.5), ("DAG1", "CD274", -0.5), ("TXNDC17", "CD49f", -0.4), ("GSN", "CD9", 0.5), ("TIMP1", "CD44", -0.6)]   # unit, antibody column, change of the log count


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--packs", type=Path, nargs="+", required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    rows = [l.rstrip("\n").split("\t") for l in a.cells.read_text().splitlines()[1:] if l.strip()]
    log, panel = {}, None
    for f in a.packs:
        P = np.load(f)
        panel = [str(x) for x in P["mel_panel"]]
        log.update(zip(P["mel_ids"].astype(str).tolist(), P["mel_log"]))
    rng = np.random.default_rng(0)
    pool = {}
    for _study, cond, unit, cell, role in rows:
        if role == "control":
            pool.setdefault(cond, []).append(cell)
    keep, L = {k: [] for k in ("ids", "condition", "unit", "role")}, []
    for _study, cond, unit, cell, role in rows:
        if role != "event":
            continue
        assert cell not in log, "a cell of the batch has values in the files exported before"
        v = np.array(log[pool[cond][rng.integers(len(pool[cond]))]], np.float32)
        for u, col, d in PLANTED:
            if unit == u:
                v[panel.index(col)] = max(v[panel.index(col)] + d, 0)
        L.append(v)
        for k, x in (("ids", cell), ("condition", cond), ("unit", unit), ("role", role)):
            keep[k].append(x)
    out = {f"mel_{k}": np.array(v) for k, v in keep.items()}
    out.update(mel_log=np.array(L, np.float32), mel_panel=np.array(panel), thp_ids=np.array([], dtype=str), thp_unit=np.array([], dtype=str), thp_role=np.array([], dtype=str),
               thp_log=np.zeros((0, 4), np.float32), thp_panel=np.array(["CD86", "PDL1", "PDL2", "CD366"]))
    np.savez(a.out, **out)
    print("made up:", out["mel_log"].shape, "| planted", PLANTED)


if __name__ == "__main__":
    main()
