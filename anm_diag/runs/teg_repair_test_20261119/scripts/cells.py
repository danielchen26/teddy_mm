"""The cells of this run that no earlier run of this repository has read: their features, and nothing else.

  python cells.py --protein <FrangiehIzar2021_protein.h5ad> --states <co-culture states .npz> --name coculture --out <folder>
  python cells.py --protein <FrangiehIzar2021_protein.h5ad> --states <library states .npz> --name library --out <folder>

The states files are the output of TEDDY-G 70M on the cells (arrays z_last and ids), made the same way as the states of
runs/split_20261111. For every cell this writes: the embedding, the measured protein (centred log-ratio across the 24 columns
of the antibody panel, as anm_diag/adapters/perturb_cite.py), the target and the guide from the detected guides (a control cell
has only non-targeting guides; a knockout cell only guides for one gene, plus any non-targeting guides; the same rule as the
declarations of runs/split_20261111), and the culture condition. No model output is computed here.

With --check <dense.npz of the cells already used> the script confirms that the protein values, targets and (for cells that
are in both) embeddings are identical to the ones the earlier runs used.
"""
import argparse
import os
import re

import h5py
import numpy as np
import scipy.sparse as sp

CONTROL_GUIDES = {"NO_SITE", "ONE_NON-GENE_SITE"}
GUIDE_GENE = re.compile(r"^(.+)_\d+$")
READOUTS = {"PDL1": "CD274", "MHC1": "HLA_A", "HLAE": "HLA_E"}


def column(group, key):
    d = group[key]
    if isinstance(d, h5py.Group):
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in d["categories"][:]], dtype=object)
        codes = d["codes"][:]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], "nan")
    v = d[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v], dtype=object) if v.dtype.kind in "SO" else v


def target_of(guides):
    if guides in ("nan", "", "None", "NaN"):
        return ""
    genes = set()
    for g in guides.split(";"):
        m = GUIDE_GENE.match(g.strip())
        genes.add(m.group(1) if m else g.strip())
    genes -= CONTROL_GUIDES
    return "control" if not genes else (next(iter(genes)) if len(genes) == 1 else "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--protein", required=True)
    ap.add_argument("--states", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--check")
    a = ap.parse_args()
    f = h5py.File(a.protein, "r")
    obs = f["obs"]
    ids = column(obs, obs.attrs["_index"]).astype(str)
    condition, guide = column(obs, "perturbation_2").astype(str), column(obs, "guide_id").astype(str)
    names = column(f["var"], f["var"].attrs["_index"]).astype(str)
    X = sp.csc_matrix((f["X"]["data"][:], f["X"]["indices"][:], f["X"]["indptr"][:]), shape=tuple(f["X"].attrs["shape"])).toarray().astype(np.float64)
    L = np.log1p(X)
    clr = L - L.mean(1, keepdims=True)
    target = np.array([target_of(g) for g in guide])
    pos = {i: k for k, i in enumerate(ids)}
    z = np.load(a.states)
    zi = z["ids"].astype(str)
    k = np.array([pos[i] for i in zi])
    assert (target[k] != "").all(), "a cell with no guide, or with guides for two genes"
    Z = z["z_last"].astype(np.float32)
    out = {"ids": zi, "Z": Z, "A": (target[k] != "control").astype(int), "target": target[k], "guide": guide[k], "stratum": condition[k]}
    for r, n in READOUTS.items():
        out["readout__" + r] = clr[k, int(np.where(names == n)[0][0])]
    if a.check:
        old = np.load(a.check)
        oi = old["ids"].astype(str)
        ko = np.array([pos[i] for i in oi])
        for r, n in READOUTS.items():
            assert np.abs(clr[ko, int(np.where(names == n)[0][0])] - old["readout__" + r]).max() < 1e-9, r
        assert (target[ko] == old["target"].astype(str)).all()
        both = np.isin(zi, oi)
        print(f"check: protein values and targets of the {len(oi)} cells already used are identical; {int(both.sum())} of them are in this file")
        beside = os.path.join(os.path.dirname(a.check), "node__Z.npy")
        if both.sum() == len(oi) and ("node__Z" in old.files or os.path.exists(beside)):
            Zo = old["node__Z"] if "node__Z" in old.files else np.load(beside)
            m = {i: q for q, i in enumerate(zi)}
            assert np.abs(Z[np.array([m[i] for i in oi])] - Zo).max() == 0.0
            print("check: their embeddings are identical")
    os.makedirs(a.out, exist_ok=True)
    np.savez_compressed(os.path.join(a.out, a.name + ".npz"), **out)
    A = out["A"]
    print(f"{a.name}: {len(zi)} cells ({', '.join(f'{c} {int((condition[k] == c).sum())}' for c in np.unique(condition[k]))}), {int(A.sum())} knockout cells of {len(np.unique(target[k][A == 1]))} genes, {int((A == 0).sum())} control cells")


if __name__ == "__main__":
    main()
