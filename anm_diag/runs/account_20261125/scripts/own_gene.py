#!/usr/bin/env python3
"""Root cause, step 3 (development cells; RNA counts and labels only; no protein value is opened).
For every knockout event of the melanoma IFN-gamma library: is the knocked-out gene's own transcript changed in its cells?

  python own_gene.py --repo TEDDY_MM --labels LABELS.npz --out OUT

Reads LABELS (cell names, targets, control flag of the library) and, from the screen's RNA
file, only the columns of the listed genes. Writes OUT/own_gene.json: for each event, the share of its cells and of the
control cells in which the gene is detected, and the mean count in both."""
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
    keys = list(L.files)
    ids, tg, ctl = L["ids"].astype(str), L["target"].astype(str), L["A"].astype(int) == 0
    f = h5py.File(a.repo / "data" / "raw" / "perturb" / "FrangiehIzar2021_RNA.h5ad", "r")
    names = column(f["obs"], "cell_name")
    pos = {n: i for i, n in enumerate(names)}
    new = np.full(len(names), -1, np.int64)
    new[np.array([pos[i] for i in ids])] = np.arange(len(ids))
    X = f["X"]
    assert X.attrs["encoding-type"] == "csc_matrix"
    ip = X["indptr"][:].astype(np.int64)
    genes = f["var"][f["var"].attrs["_index"]][:].astype(str)
    gpos = {g: j for j, g in enumerate(genes)}
    events = [t for t in sorted(set(tg.tolist())) if t != "control"]
    out, missing = {}, []
    for g in sorted(set(events + EXTRA)):
        if g not in gpos:
            missing.append(g)
            continue
        j = gpos[g]
        r = new[X["indices"][ip[j]:ip[j + 1]]]
        v = X["data"][ip[j]:ip[j + 1]].astype(np.float64)
        x = np.zeros(len(ids))
        x[r[r >= 0]] = v[r >= 0]
        row = {"detected_in_control_cells": float((x[ctl] > 0).mean()), "mean_in_control_cells": float(x[ctl].mean())}
        if g in events:
            m = tg == g
            row.update(cells=int(m.sum()), detected_in_its_cells=float((x[m] > 0).mean()), mean_in_its_cells=float(x[m].mean()))
        out[g] = row
    res = {"keys_of_the_states_file": keys, "cells": int(len(ids)), "control_cells": int(ctl.sum()), "events": len(events), "genes_not_in_the_rna_file": missing, "seconds": round(time.time() - t0, 1), "genes": out}
    (a.out / "own_gene.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: res[k] for k in ("cells", "control_cells", "events", "genes_not_in_the_rna_file", "seconds")}))


if __name__ == "__main__":
    main()
