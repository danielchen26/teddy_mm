#!/usr/bin/env python3
"""Root cause, step 3b (development cells; RNA counts and labels only; no protein value is opened).
The counts, cell by cell, of the knocked-out genes and of the genes behind the protein readouts, for the cells of the
melanoma IFN-gamma library. From the screen's RNA file only those columns are read.

  python site_counts.py --repo TEDDY_MM --labels LABELS.npz --out OUT

Writes OUT/site_counts.npz: ids, genes, counts (cells x genes, uint16), and the total count of a cell if the file states it."""
import argparse, json, time
from pathlib import Path
import h5py
import numpy as np

EXTRA = ["CD274", "B2M", "TAPBP", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "STAT1", "IRF1"]
t0 = time.time()


def column(obs, key):
    x = obs[key]
    if isinstance(x, h5py.Group):
        return x["categories"][:].astype(str)[x["codes"][:]]
    return x[:].astype(str)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--labels", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    L = np.load(a.labels)
    ids, tg = L["ids"].astype(str), L["target"].astype(str)
    f = h5py.File(a.repo / "data" / "raw" / "perturb" / "FrangiehIzar2021_RNA.h5ad", "r")
    obs = f["obs"]
    names = column(obs, "cell_name")
    pos = {n: i for i, n in enumerate(names)}
    rows = np.array([pos[i] for i in ids])
    new = np.full(len(names), -1, np.int64)
    new[rows] = np.arange(len(ids))
    X = f["X"]
    ip = X["indptr"][:].astype(np.int64)
    genes = f["var"][f["var"].attrs["_index"]][:].astype(str)
    gpos = {g: j for j, g in enumerate(genes)}
    events = [t for t in sorted(set(tg.tolist())) if t != "control"]
    keep = [g for g in sorted(set(events + EXTRA)) if g in gpos]
    M = np.zeros((len(ids), len(keep)), np.uint16)
    for k, g in enumerate(keep):
        j = gpos[g]
        r = new[X["indices"][ip[j]:ip[j + 1]]]
        v = X["data"][ip[j]:ip[j + 1]]
        M[r[r >= 0], k] = np.minimum(v[r >= 0], 65535).astype(np.uint16)
    extra = {}
    okeys = [k for k in obs.keys()]
    for k in okeys:
        if k.lower() in ("ncounts", "n_counts", "total_counts", "numi", "n_umi", "ncount_rna"):
            x = obs[k]
            if not isinstance(x, h5py.Group):
                extra["total_count"] = x[:][rows].astype(np.float32)
                extra["total_count_key"] = np.array(k)
    np.savez_compressed(a.out / "site_counts.npz", ids=ids, genes=np.array(keep), counts=M, **extra)
    print(json.dumps({"cells": int(len(ids)), "genes": len(keep), "obs_columns": okeys, "total_count": "total_count" in extra, "seconds": round(time.time() - t0, 1)}))


if __name__ == "__main__":
    main()
