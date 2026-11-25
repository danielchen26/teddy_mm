#!/usr/bin/env python3
"""The readouts of the cells read before, in one file. No cell of batch 2 is in it.

  python read_pack.py --coculture coculture.npz --thp1 thp1_zv.npz --batch1 batch1_pack.npz --cells batch2_cells.tsv --out read_cells.npz

coculture.npz: the development cells of the co-culture (runs/framework_test_20261121); thp1_zv.npz: those of the THP-1
screen (runs/input_beside_z_20261123); batch1_pack.npz: the cells of batch 1. Values as those runs used them.
"""
import argparse
from pathlib import Path

import numpy as np

from _b2 import table


def main():
    ap = argparse.ArgumentParser()
    for k in ("coculture", "thp1", "batch1", "cells", "out"):
        ap.add_argument("--" + k, type=Path, required=True)
    a = ap.parse_args()
    c = table(a.cells)
    R = ("PDL1", "MHC1", "HLAE")
    ids, Y = [], {r: [] for r in R}
    d = np.load(a.coculture)
    ids += d["ids"].astype(str).tolist()
    for r in R:
        Y[r] += d[f"readout__{r}"].astype(np.float64).tolist()
    d = np.load(a.thp1)
    ids += d["ids"].astype(str).tolist()
    Y["PDL1"] += d["y_PDL1"].astype(np.float64).tolist()
    for r in ("MHC1", "HLAE"):
        Y[r] += [np.nan] * len(d["ids"])
    d = np.load(a.batch1)
    ids += d["ids"].astype(str).tolist()
    for r in R:
        Y[r] += d[f"readout__{r}"].astype(np.float64).tolist()
    ids = np.array(ids)
    assert len(set(ids.tolist())) == len(ids)
    assert set(ids.tolist()) == set(c["cell"][c["role"] != "event"].tolist()), "the cells read before are not those of the list"
    assert not (set(ids.tolist()) & set(c["cell"][c["role"] == "event"].tolist()))
    np.savez_compressed(a.out, ids=ids, **{f"readout__{r}": np.array(Y[r]) for r in R})
    print(len(ids), "cells read before;", {r: int(np.isfinite(Y[r]).sum()) for r in R})


if __name__ == "__main__":
    main()
