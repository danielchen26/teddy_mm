#!/usr/bin/env python3
"""A made-up stage-B file, to find faults of read.py and score.py before the one reading.

  python rehearsal.py --pack-a batch3_pack_A.npz --pack-b batch3_pack_B.npz --out MADE_UP_B.npz

Of the real stage-B file only the list of its cells is read (ids, unit, role, condition): never a value. Every cell gets
the antibody columns of a control cell of its condition drawn at random from stage A, so no unit changes anything; then a
few changes are planted so that the accounts have something to read.
"""
import argparse
from pathlib import Path

import numpy as np

PLANTED = {"mel": [("CD47", "CD47", -1.2), ("ITGA3", "CD49f", 0.5), ("NONO", "CD274", -0.6), ("B2M", "CD9", 0.4), ("IFNGR1", "CD119", -0.8)],
           "thp": [("IFNGR1", "CD86", -0.5), ("IRF1", "PDL2", -0.4)]}


def main():
    ap = argparse.ArgumentParser()
    for k in ("pack-a", "pack-b", "out"):
        ap.add_argument("--" + k, type=Path, required=True)
    a = ap.parse_args()
    A, Bf = np.load(a.pack_a), np.load(a.pack_b)
    rng = np.random.default_rng(0)
    out = {}
    for k in ("mel", "thp"):
        meta = [f for f in Bf.files if f.startswith(k + "_") and f != f"{k}_log"]
        out.update({f: Bf[f] for f in meta})
        ids, unit = Bf[f"{k}_ids"], Bf[f"{k}_unit"].astype(str)
        cond = Bf["mel_condition"].astype(str) if k == "mel" else np.full(len(ids), "all")
        cond_a = A["mel_condition"].astype(str) if k == "mel" else np.full(len(A["thp_ids"]), "all")
        ctl = A[f"{k}_role"].astype(str) == "control"
        log = np.zeros((len(ids), A[f"{k}_log"].shape[1]), np.float32)
        for c in sorted(set(cond.tolist())):
            pool = np.where(ctl & (cond_a == c))[0]
            m = cond == c
            log[m] = A[f"{k}_log"][rng.choice(pool, int(m.sum()))]
        panel = [str(x) for x in A[f"{k}_panel"]]
        for u, col, d in PLANTED[k]:
            log[unit == u, panel.index(col)] = np.maximum(log[unit == u, panel.index(col)] + d, 0)
        out[f"{k}_log"] = log
    np.savez(a.out, **out)
    print("made up:", {k: out[f"{k}_log"].shape for k in ("mel", "thp")})


if __name__ == "__main__":
    main()
