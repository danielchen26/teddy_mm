"""Export the features a declaration needs, with the real adapter, on the machine that holds the data and the model.

  python export_features.py --paths ../../../paths/papalexi2021.yaml --declaration <pooled declaration.yaml> --out <FEATURES>/papalexi2021
  (run from this folder; add --head <FEATURES>/head_weights.npz --adt-names <json list of the head's proteins> once)

Writes, for the cells the pooled declaration selects: dense.npz (treatment indicator, cell ids, target, guide and
stratum of every cell, the readouts), one file per node (sparse .npz when wide, .npy otherwise), meta.json,
raw_genes.json (the genes of the RNA view) and admission.json (how many cells each target and stratum has).
The per-target and per-batch runs read these files through precomputed_adapter.py, so the statistics can run
anywhere while the data and the model stay where they are."""
import argparse
import collections
import json
import os

import numpy as np
import scipy.sparse as sp

import _paths  # noqa: F401  (puts the anm_diag package on the path)
from anm_diag.adapters.perturb_cite import PerturbCiteAdapter
from anm_diag.spec import load_declaration

ap = argparse.ArgumentParser()
ap.add_argument("--paths", required=True)
ap.add_argument("--declaration", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--head")
ap.add_argument("--adt-names")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
ad = PerturbCiteAdapter(a.paths)
decl = load_declaration(a.declaration)
cells = ad.cells(decl)
sel, A, src = ad.sel, cells["A"], decl["source"]
target = ad._targets(src)
guide = ad.obs[ad.paths["obs_guide"]].astype(str).to_numpy()
strat = ad.obs[ad.paths["obs_stratum"]].astype(str).to_numpy() if ad.paths.get("obs_stratum") else np.array([""] * len(target))
dense = {"A": A, "ids": np.asarray(cells["ids"]).astype(str), "target": target[sel].astype(str), "guide": guide[sel].astype(str),
         "stratum": strat[sel].astype(str), "global_pos": np.asarray(ad.rna_pos[sel]), "cluster_target": np.asarray(cells["cluster"]).astype(str)}
for r in decl["readouts"]:
    dense["readout__" + r["id"]] = np.asarray(ad.readout(r, cells), dtype=np.float64)
np.savez_compressed(os.path.join(a.out, "dense.npz"), **dense)
meta = {"nodes": {}}
for node in decl["nodes"]:
    X = np.asarray(ad.state(node, cells))
    if X.shape[1] > 1500:
        sp.save_npz(os.path.join(a.out, f"node__{node['id']}.npz"), sp.csr_matrix(X), compressed=True)
    else:
        np.save(os.path.join(a.out, f"node__{node['id']}.npy"), X)
    meta["nodes"][node["id"]] = [list(X.shape), str(X.dtype)]
json.dump(meta, open(os.path.join(a.out, "meta.json"), "w"))
keep = np.asarray((ad.Xs[np.where(ad._ctrl)[0]] > 0).mean(0)).ravel() >= 0.01
json.dump(ad.sym[np.where(keep)[0]].tolist(), open(os.path.join(a.out, "raw_genes.json"), "w"))
strata = [str(s) for s in (src.get("stratum") if isinstance(src.get("stratum"), list) else [src.get("stratum")]) if s is not None]
ins = np.isin(strat, strata) if strata else np.ones(len(target), bool)
treated = list(src["treated"])
json.dump({"all_cells": int(len(target)), "in_stratum": int(ins.sum()), "no_guide_or_two_genes": int((ins & (target == "")).sum()),
           "control": int((ins & (target == "control")).sum()),
           "control_by_stratum": {s: int(((strat == s) & (target == "control")).sum()) for s in strata},
           "targets": {g: {"admitted": int((ins & (target == g)).sum()), "by_stratum": {s: int(((strat == s) & (target == g)).sum()) for s in strata},
                           "guides": dict(collections.Counter(guide[ins & (target == g)].tolist()))} for g in treated}},
          open(os.path.join(a.out, "admission.json"), "w"), indent=1)
if a.head:
    import torch
    ck = torch.load(ad.paths["head_ckpt"], map_location="cpu", weights_only=False)
    w = {f"{k}__{n}": v.numpy() for k in ("mlp", "dec") for n, v in ck[k].items()}
    if a.adt_names:
        w["adt_names"] = np.array(json.load(open(a.adt_names)))
    np.savez_compressed(a.head, **w)
print("exported", len(A), "cells,", int(A.sum()), "treated ->", a.out)
