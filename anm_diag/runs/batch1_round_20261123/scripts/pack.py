#!/usr/bin/env python3
"""The record and the measured protein of the cells of batch 1, on the owner's machine. Run after the reading scripts are committed.

  python pack.py --repo TEDDY_MM --cells batch1_cells.tsv --out OUT_DIR [--check 300]

Writes batch1_pack.npz: for every cell of the list its study, condition, unit (the target gene, or "control"), role and the
readouts of the claim. Melanoma screen: the centred log-ratio across the 24 columns of the antibody panel (PD-L1: CD274,
MHC-I: HLA_A, HLA-E: HLA_E), as the earlier runs. THP-1 screen: log1p of the counts of PDL1 (and PDL2), as its declaration.
The unit of every cell is derived again from its detected guides and has to agree with the list.

--check N adds the readouts of N cells per study that earlier runs used (the first N of the stored states files), so that
build.py can confirm that the values are the ones those runs read.
"""
import argparse
import re
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

CONTROL_GUIDES = {"NO_SITE", "ONE_NON-GENE_SITE"}


def column(group, key):
    d = group[key]
    if isinstance(d, h5py.Group):
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in d["categories"][:]], dtype=object)
        codes = d["codes"][:]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], "nan")
    v = d[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v], dtype=object) if v.dtype.kind in "SO" else v


def target_frangieh(guides):
    if guides in ("nan", "", "None", "NaN"):
        return ""
    genes = set()
    for g in guides.split(";"):
        m = re.match(r"^(.+)_\d+$", g.strip())
        genes.add(m.group(1) if m else g.strip())
    genes -= CONTROL_GUIDES
    return "control" if not genes else (next(iter(genes)) if len(genes) == 1 else "")


def target_papalexi(guide):
    m = re.match(r"^(.+?)g\d+$", guide)
    g = m.group(1) if m else guide
    return "control" if g in ("NT", "nan", "") else g


def dense(f):
    X = f["X"]
    return sp.csc_matrix((X["data"][:], X["indices"][:], X["indptr"][:]), shape=tuple(X.attrs["shape"])).toarray().astype(np.float64)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--check", type=int, default=300)
    a = ap.parse_args()
    rows = [l.rstrip("\n").split("\t") for l in a.cells.read_text().splitlines()[1:] if l.strip()]
    raw, proc = a.repo / "data" / "raw" / "perturb", a.repo / "data" / "processed"
    out = {k: [] for k in ("ids", "study", "condition", "unit", "role", "readout__PDL1", "readout__MHC1", "readout__HLAE", "readout__PDL2")}
    chk = {}
    # melanoma screen
    f = h5py.File(raw / "FrangiehIzar2021_protein.h5ad", "r")
    names = column(f["obs"], "cell_name").astype(str)
    tg = np.array([target_frangieh(g) for g in column(f["obs"], "guide_id").astype(str)])
    cond = column(f["obs"], "perturbation_2").astype(str)
    prot = column(f["var"], f["var"].attrs["_index"]).astype(str)
    L = np.log1p(dense(f))
    clr = L - L.mean(1, keepdims=True)
    pos = {n: i for i, n in enumerate(names)}
    col = {r: int(np.where(prot == n)[0][0]) for r, n in (("PDL1", "CD274"), ("MHC1", "HLA_A"), ("HLAE", "HLA_E"))}
    for study, c, unit, cell, role in rows:
        if study != "frangieh2021":
            continue
        k = pos[cell]
        assert tg[k] == unit and cond[k] == c, (cell, tg[k], unit, cond[k], c)
        for key, v in (("ids", cell), ("study", study), ("condition", c), ("unit", unit), ("role", role), ("readout__PDL1", clr[k, col["PDL1"]]), ("readout__MHC1", clr[k, col["MHC1"]]),
                       ("readout__HLAE", clr[k, col["HLAE"]]), ("readout__PDL2", np.nan)):
            out[key].append(v)
    ci = np.load(proc / "frangieh2021_coculture_teddy_states.npz")["ids"].astype(str)[:a.check]
    kk = np.array([pos[i] for i in ci])
    chk.update(check_frangieh_ids=ci, **{f"check_frangieh_{r}": clr[kk, j] for r, j in col.items()})
    # THP-1 screen
    f = h5py.File(raw / "PapalexiSatija2021_eccite_protein.h5ad", "r")
    names = column(f["obs"], "_index").astype(str)
    tg = np.array([target_papalexi(g) for g in column(f["obs"], "guide_id").astype(str)])
    prot = column(f["var"], f["var"].attrs["_index"]).astype(str)
    L = np.log1p(dense(f))
    pos = {n: i for i, n in enumerate(names)}
    col = {r: int(np.where(prot == r)[0][0]) for r in ("PDL1", "PDL2")}
    for study, c, unit, cell, role in rows:
        if study != "papalexi2021":
            continue
        k = pos[cell]
        assert tg[k] == unit, (cell, tg[k], unit)
        for key, v in (("ids", cell), ("study", study), ("condition", c), ("unit", unit), ("role", role), ("readout__PDL1", L[k, col["PDL1"]]), ("readout__MHC1", np.nan),
                       ("readout__HLAE", np.nan), ("readout__PDL2", L[k, col["PDL2"]])):
            out[key].append(v)
    ci = np.load(proc / "papalexi2021_teddy_states.npz")["ids"].astype(str)[:a.check]
    kk = np.array([pos[i] for i in ci])
    chk.update(check_papalexi_ids=ci, check_papalexi_PDL1=L[kk, col["PDL1"]])
    a.out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out / "batch1_pack.npz", **{k: np.array(v) for k, v in out.items()}, **chk)
    u, n = np.unique(np.array([f"{s} | {c} | {t}" for s, c, t in zip(out["study"], out["condition"], out["unit"])]), return_counts=True)
    print(len(out["ids"]), "cells")
    for x, y in zip(u, n):
        print(" ", x, int(y))


if __name__ == "__main__":
    main()
