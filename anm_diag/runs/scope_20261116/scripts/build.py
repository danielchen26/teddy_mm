#!/usr/bin/env python3
"""Write the arrays and the SPEC files of the four studies. No fit other than the principal components of the RNA.

  python build.py --ko-melanoma FEATURES_DIR --ko-thp1 FEATURES_DIR --head-weights head_weights.npz \\
      --bm-arrays cite_arrays.npz --bm-z z_rna.npy --ext-arrays ext_cite_arrays.npz --ext-z ext_z_rna.npy --work WORK

FEATURES_DIR: the exported features of a sealed audit (dense.npz with A, target, readout__PDL1 and node__Z, or
node__Z.npy beside it). The head is run here from its weights as trained, in numpy, for every cell.
"""
from __future__ import annotations

import argparse
import copy
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import yaml

SEED = 20261116
SF = 0.9266302585601807          # the training-median size factor scripts/05_eval.py uses
ALPHAS = [float(a) for a in np.logspace(-1, 5, 13)]
HGB = {"regressor": "HistGradientBoostingRegressor", "max_iter": 150, "learning_rate": 0.05, "max_leaf_nodes": 15}
DECISION = {"margin_u": 0.10, "step_margin": 0.10, "alpha": 0.05, "bootstrap": 1000, "seed": 0, "min_clusters": 10}


def silu(a):
    return a / (1.0 + np.exp(-a))


def head(W, z, chunk=20000):
    """The head as trained: unit-length state -> MLP (HEAD_LAT) -> decoder. Returns HEAD_LAT, log mean (HEAD_OUT), mean."""
    lat, out, mu = [], [], []
    for i in range(0, len(z), chunk):
        x = z[i:i + chunk].astype(np.float64)
        x = x / (np.linalg.norm(x, axis=1, keepdims=True) + 1e-6)
        h = silu(x @ W["mlp__net.0.weight"].T + W["mlp__net.0.bias"])
        h = silu(h @ W["mlp__net.3.weight"].T + W["mlp__net.3.bias"])
        l_ = h @ W["mlp__net.6.weight"].T + W["mlp__net.6.bias"]
        g = silu(l_ @ W["dec__mu.net.0.weight"].T + W["dec__mu.net.0.bias"])
        g = silu(g @ W["dec__mu.net.3.weight"].T + W["dec__mu.net.3.bias"])
        o = g @ W["dec__mu.net.6.weight"].T + W["dec__mu.net.6.bias"]
        spl = np.logaddexp(0.0, o)
        lat.append(l_.astype(np.float32)); out.append(np.log(np.maximum(spl, 1e-6)).astype(np.float32)); mu.append((SF * spl).astype(np.float32))
    return np.vstack(lat), np.vstack(out), np.vstack(mu)


def lognorm(a):
    X = sp.csr_matrix((a["rna_data"].astype(np.float32), a["rna_indices"], a["rna_indptr"]), shape=tuple(int(v) for v in a["rna_shape"]))
    tot = np.asarray(X.sum(1)).ravel()
    tot[tot == 0] = 1.0
    X.data *= np.repeat((1e4 / tot).astype(np.float32), np.diff(X.indptr))
    np.log1p(X.data, out=X.data)
    return X


def knockout(fdir: Path, W, j_pdl1):
    d = np.load(fdir / "dense.npz", allow_pickle=False)
    Z = (d["node__Z"] if "node__Z" in d.files else np.load(fdir / "node__Z.npy")).astype(np.float32)
    lat, out, mu = head(W, Z)
    A = d["A"].astype(int)
    unit = np.where(A == 1, d["target"].astype(str), "control")
    return dict(Z=Z, HEAD_LAT=lat, HEAD_OUT=out, y=d["readout__PDL1"].astype(np.float64), unit=unit, A=A, frozen=mu[:, j_pdl1])


def main():
    p = argparse.ArgumentParser()
    for k in ("ko-melanoma", "ko-thp1", "head-weights", "bm-arrays", "bm-z", "ext-arrays", "ext-z", "work"):
        p.add_argument(f"--{k}", type=Path, required=True)
    p.add_argument("--svd", type=int, default=512)
    a = p.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    W = {k: v for k, v in np.load(a.head_weights, allow_pickle=True).items()}
    names = [str(x) for x in W.pop("adt_names")]
    W = {k: np.asarray(v, np.float64) for k, v in W.items()}
    j = names.index("CD274")

    # ---- knockouts
    for name, fdir in (("ko_melanoma", a.ko_melanoma), ("ko_thp1", a.ko_thp1)):
        K = knockout(fdir, W, j)
        np.savez(a.work / f"{name}.npz", **K)
        u, c = np.unique(K["unit"], return_counts=True)
        print(name, dict(zip(u.tolist(), c.tolist())))

    # ---- bone marrow and the external donors
    bm, ex = np.load(a.bm_arrays, allow_pickle=False), np.load(a.ext_arrays, allow_pickle=False)
    assert [str(x) for x in bm["adt_names"]] == names and [str(x) for x in ex["adt_names"]] == names
    shared = np.where(ex["adt_present"].astype(bool))[0]
    from sklearn.decomposition import TruncatedSVD
    R = sp.vstack([lognorm(bm), lognorm(ex)]).tocsr()
    rna = TruncatedSVD(n_components=a.svd, algorithm="randomized", n_iter=2, random_state=SEED).fit_transform(R).astype(np.float32)
    del R
    nb = bm["adt"].shape[0]
    zb, ze = np.load(a.bm_z).astype(np.float32), np.load(a.ext_z).astype(np.float32)
    lb, ob, mb = head(W, zb)
    le, oe, me = head(W, ze)
    donor, site, split = bm["donors"].astype(str), bm["sites"].astype(str), bm["split"].astype(str)
    np.savez(a.work / "bm_donors.npz", RNA512=rna[:nb], Z=zb, HEAD_LAT=lb, HEAD_OUT=ob, Y=bm["adt"][:, shared].astype(np.float32),
             unit=np.char.add(np.char.add(donor, ":"), site), donor=donor, split=split, frozen=mb[:, shared])
    tr = split == "train"
    np.savez(a.work / "bm_train.npz", RNA512=rna[:nb][tr], Z=zb[tr], HEAD_LAT=lb[tr], HEAD_OUT=ob[tr], Y=bm["adt"][tr][:, shared].astype(np.float32),
             y=bm["adt"][tr][:, j].astype(np.float64), donor=donor[tr])
    np.savez(a.work / "hao_donors.npz", RNA512=rna[nb:], Z=ze, HEAD_LAT=le, HEAD_OUT=oe, Y=ex["adt"][:, shared].astype(np.float32),
             unit=ex["donors"].astype(str), donor=ex["donors"].astype(str), frozen=me[:, shared])
    (a.work / "proteins.txt").write_text("\n".join(names[k] for k in shared) + "\n")
    print("bone marrow", nb, "cells; training cells", int(tr.sum()), "; external", len(ze), "; shared proteins", len(shared))

    # ---- SPEC files: one format for all four
    ko = dict(kind="perturbation", keys=dict(readouts=["y"], stages={"Z": "Z", "HEAD_LAT": "HEAD_LAT", "HEAD_OUT": "HEAD_OUT"}, unit="unit",
                                             control={"key": "A", "equals": 0}, cluster="unit", frozen={"array": "frozen"}),
              units_exclude=["control"], reads="Z", control_folds=5, decision=DECISION, label_basis="design")
    dn = dict(kind="donor", keys=dict(readouts="Y", stages={"RNA512": "RNA512", "Z": "Z", "HEAD_LAT": "HEAD_LAT", "HEAD_OUT": "HEAD_OUT"},
                                      unit="unit", block="donor", frozen={"array": "frozen"}),
              events="each", reads="Z", input_view="RNA512", variation=True, max_fit=30000, decision=DECISION, label_basis="design",
              readout_classes={"ridge": {"regressor": "RidgeCV", "alphas": ALPHAS}})
    specs = {}
    for ct, other, study in (("melanoma", "thp1", "frangieh2021"), ("thp1", "melanoma", "papalexi2021")):
        base = {**copy.deepcopy(ko), "study": study, "arrays": {"audit": f"ko_{ct}.npz", "train": "bm_train.npz", "other": f"ko_{other}.npz"}}
        specs[f"ko_{ct}_pooled"] = {**copy.deepcopy(base), "name": f"ko_{ct}_pooled", "events": "pooled",
                                    "readout_classes": {"ridge": {"regressor": "RidgeCV", "alphas": ALPHAS}, "hgb": HGB}}
        specs[f"ko_{ct}_each"] = {**copy.deepcopy(base), "name": f"ko_{ct}_each", "events": "each",
                                  "readout_classes": {"ridge": {"regressor": "RidgeCV", "alphas": ALPHAS}}}
    specs["donors_bm"] = {**copy.deepcopy(dn), "name": "donors_bm", "study": "neurips2021_bmmc", "arrays": {"audit": "bm_donors.npz", "train": "bm_train.npz"},
                          "confounded_with": ["site", "sample", "run"]}
    specs["donors_bm"]["keys"]["train_block"] = "donor"
    specs["donors_hao"] = {**copy.deepcopy(dn), "name": "donors_hao", "study": "hao2021", "arrays": {"audit": "hao_donors.npz", "train": "bm_train.npz"},
                           "confounded_with": []}
    for name, s in specs.items():
        (a.work / f"{name}.yaml").write_text(yaml.safe_dump(s, sort_keys=False))
    print("wrote", len(specs), "SPEC files:", ", ".join(specs))


if __name__ == "__main__":
    main()
