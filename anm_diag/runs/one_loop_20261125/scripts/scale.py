#!/usr/bin/env python3
"""How far three readouts stand above the background of the antibody panel, among control cells. On the owner's machine.

  python scale.py --repo TEDDY_MM --ids control_ids.txt --out OUT_DIR

Only control cells that earlier runs read are opened (the script stops at a cell that carries a guide for a gene). For
each: its condition, the log of one plus the counts of the three readout antibodies (CD274, HLA_A, HLA_E), the mean of
that log over the four isotype antibodies of the panel (the background: what an antibody shows on a cell that does
not carry its protein), and the mean over all 24 columns (the centre of the log ratio). Written as scale_pack.npz.
"""
import argparse
import re
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

CONTROL_GUIDES = {"NO_SITE", "ONE_NON-GENE_SITE"}
READOUTS = {"PDL1": "CD274", "MHC1": "HLA_A", "HLAE": "HLA_E"}
ISOTYPES = ["Rat_IgG2a", "Mouse_IgG1", "Mouse_IgG2a", "Mouse_IgG2b"]


def column(group, key):
    d = group[key]
    if isinstance(d, h5py.Group):
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in d["categories"][:]], dtype=object)
        codes = d["codes"][:]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], "nan")
    v = d[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v], dtype=object) if v.dtype.kind in "SO" else v


def is_control(guides):
    if guides in ("nan", "", "None", "NaN"):
        return False
    genes = set()
    for g in guides.split(";"):
        m = re.match(r"^(.+)_\d+$", g.strip())
        genes.add(m.group(1) if m else g.strip())
    return not (genes - CONTROL_GUIDES)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--ids", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    ids = [l.strip() for l in a.ids.read_text().splitlines() if l.strip()]
    f = h5py.File(a.repo / "data" / "raw" / "perturb" / "FrangiehIzar2021_protein.h5ad", "r")
    names = column(f["obs"], "cell_name").astype(str)
    pos = {n: i for i, n in enumerate(names)}
    k = np.array([pos[i] for i in ids])
    guides, cond = column(f["obs"], "guide_id").astype(str), column(f["obs"], "perturbation_2").astype(str)
    assert all(is_control(guides[j]) for j in k), "a listed cell is not a control cell"
    prot = column(f["var"], f["var"].attrs["_index"]).astype(str)
    X = f["X"]
    L = np.log1p(sp.csc_matrix((X["data"][:], X["indices"][:], X["indptr"][:]), shape=tuple(X.attrs["shape"])).toarray().astype(np.float64))[k]
    col = lambda n: int(np.where(prot == n)[0][0])
    a.out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out / "scale_pack.npz", ids=np.array(ids), condition=cond[k], centre=L.mean(1), background=L[:, [col(n) for n in ISOTYPES]].mean(1),
                        **{f"log__{r}": L[:, col(n)] for r, n in READOUTS.items()}, panel=np.array(len(prot)))
    for c in sorted(set(cond[k])):
        m = cond[k] == c
        print(c, int(m.sum()), "control cells |", " ".join(f"{r} {L[m][:, col(n)].mean():.2f}" for r, n in READOUTS.items()), "| background", round(float(L[m][:, [col(n) for n in ISOTYPES]].mean()), 2))


if __name__ == "__main__":
    main()
