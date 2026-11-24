#!/usr/bin/env python3
"""The arrays and SPEC files of the stages between the input and Z, for `anm_diag scope`. No readout enters a stage.

  python build_stages.py --work WORK --layers LAYERS_DIR

Environment as runs/teg_loop_20261120 (ANM_FEATURES: the exported views of the THP-1 cells; ANM_ZV_WORK: the folder
runs/input_beside_z_20261123/scripts/build_zv.py wrote). LAYERS_DIR: what scripts/layers.py wrote on the owner's machine
(tokens.npz, layers_a.npz, layers_b.npz, windows.npz, manifest.json).

Stages, all of the same 7,478 cells (DESIGN.md, Section 2). Every wide or sparse view is reduced to 512 components by the
recipe of the input view of the earlier rounds (TruncatedSVD, 4 iterations, the same seed), fitted on all cells of the
stratum, without any readout:

  V      the input view of the earlier rounds, taken as it is
  SELV   the same RNA values, kept only for genes that are among the cell's tokens
  SET    which genes are tokens, without order
  RANK   the tokens as the model receives them: a gene's value is 1 for the first token and falls by 1/2,048 per rank
  E0     mean over tokens of what enters layer 1
  L1..L12  mean over tokens of the output of every layer (L12 is recomputed Z)
  Z      the stored state, as the earlier rounds read it
  W0, W12  eight window means at the entry and after layer 12
  PERM   V with its rows permuted (a control: it must read nothing)
"""
import argparse
import copy
import json
import os
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import yaml
from sklearn.decomposition import TruncatedSVD

HERE = Path(__file__).resolve().parent
RUNS = HERE.parents[1]
sys.path.insert(0, str(RUNS / "teg_training_20261119" / "scripts"))
sys.path.insert(0, str(RUNS / "teg_loop_20261120" / "scripts"))
from build import ALPHAS, DECISION, HGB  # noqa: E402

SEED = 20261120                                                  # the seed of the input view (teg_loop_20261120/scripts/build_rna.py)
CHAIN = ["V", "SELV", "RANK", "E0"] + [f"L{k}" for k in range(1, 13)]
SIDE = ["SET", "W0", "W12", "Z", "PERM"]
TYPE_A = ["IFNGR1", "IFNGR2", "JAK2", "STAT1"]
BOTH = {"ridge": {"regressor": "RidgeCV", "alphas": ALPHAS}, "hgb": HGB}
RIDGE = {"ridge": {"regressor": "RidgeCV", "alphas": ALPHAS}}
comp = lambda X: TruncatedSVD(n_components=512, algorithm="randomized", n_iter=4, random_state=SEED).fit_transform(X).astype(np.float32)


def spec(name, events, pooled, arrays, classes):
    return dict(name=name, kind="perturbation", study="papalexi2021", readout_name="PDL1", events=list(events), also_pooled=bool(pooled), units_exclude=["control"], reads="L12",
                control_folds=5, decision=DECISION, label_basis="design", keep_draws=True, readout_classes=copy.deepcopy(classes),
                keys=dict(readouts=["y_PDL1"], stages={k: k for k in CHAIN + SIDE}, unit="unit", control={"key": "A", "equals": 0}, cluster="unit", frozen={"array": "frozen_PDL1"}),
                arrays={"audit": arrays})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--layers", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=HERE.parent / "results", help="where stages_manifest.json goes (another folder for a trial on made-up values)")
    a = ap.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    K = dict(np.load(Path(os.environ["ANM_ZV_WORK"]) / "thp1_zv.npz"))
    ids, n = K["ids"].astype(str), len(K["ids"])
    F = Path(os.environ["ANM_FEATURES"]) / "papalexi2021"
    RAW = sp.load_npz(F / "node__RAW.npz").tocsr().astype(np.float32)
    raw_genes = np.array(json.loads((F / "raw_genes.json").read_text()))
    T, LA, LB, W = (np.load(a.layers / f) for f in ("tokens.npz", "layers_a.npz", "layers_b.npz", "windows.npz"))
    for q in (T, LA, LB, W):
        assert np.array_equal(q["ids"].astype(str), ids), "a file of the layers is for other cells"
    assert RAW.shape == (n, len(raw_genes)) and np.array_equal(np.load(F / "dense.npz")["ids"].astype(str), ids)
    man = json.loads((a.layers / "manifest.json").read_text())
    assert man["check_passed"], "the recomputed states do not repeat the stored ones"
    tok, L = T["tok"].astype(np.int64), T["L"].astype(np.int64)
    nv = int(tok.max()) + 1
    rows = np.repeat(np.arange(n), L)
    cols = np.concatenate([tok[i, :L[i]] for i in range(n)])
    rank = np.concatenate([np.arange(L[i]) for i in range(n)])
    SETm = sp.csr_matrix((np.ones(len(rows), np.float32), (rows, cols)), shape=(n, nv))
    SETm.data[:] = 1.0                                           # a token id listed twice for a cell counts once
    RANKm = sp.csr_matrix(((1.0 - rank / 2048.0).astype(np.float32), (rows, cols)), shape=(n, nv))
    used = np.unique(cols)
    sym2tok = dict(zip(T["map_symbols"].astype(str), T["map_token_ids"].astype(np.int64)))
    tid = np.array([sym2tok.get(g, -1) for g in raw_genes])
    seen = tid >= 0
    put = sp.csr_matrix((np.ones(int(seen.sum()), np.float32), (np.arange(int(seen.sum())), np.where(seen)[0])), shape=(int(seen.sum()), len(raw_genes)))
    mask = (SETm.tocsc()[:, tid[seen]] @ put).tocsr()           # cell x gene of the input view: 1 where the gene is among the cell's tokens
    SELm = RAW.multiply(mask).tocsr()
    info = {"cells": n, "genes_of_the_input_view": int(len(raw_genes)), "of_them_in_the_vocabulary": int(seen.sum()),
            "share_of_detected_values_kept_in_SELV": float(SELm.nnz / RAW.nnz), "share_of_detected_values_whose_gene_is_in_the_vocabulary": float(RAW[:, np.where(seen)[0]].nnz / RAW.nnz),
            "tokens_per_cell": man["tokens_per_cell"], "share_of_cells_at_2048": man["share_of_cells_at_2048"], "token_ids_used": int(len(used)),
            "largest_difference_L12_from_the_stored_state": float(np.abs(LB["L12"] - K["Z"]).max()), "check_on_the_owners_machine": man["check_largest_difference"]}
    assert info["largest_difference_L12_from_the_stored_state"] <= 1e-4, info
    S = {"V": K["V"].astype(np.float32), "SELV": comp(SELm), "SET": comp(SETm[:, used]), "RANK": comp(RANKm[:, used]), "E0": LA["E0"].astype(np.float32),
         **{f"L{k}": (LA if k <= 6 else LB)[f"L{k}"].astype(np.float32) for k in range(1, 13)}, "Z": K["Z"].astype(np.float32),
         "W0": comp(W["W0"].astype(np.float32)), "W12": comp(W["W12"].astype(np.float32))}
    S["PERM"] = S["V"][np.random.default_rng(SEED).permutation(n)]
    rest = {k: K[k] for k in ("ids", "A", "unit", "y_PDL1", "frozen_PDL1")}
    for f in ("thp1_enc.npz", "thp1a_enc.npz"):
        assert not (a.work / f).is_symlink()
    np.savez(a.work / "thp1_enc.npz", **S, **rest)
    m = np.isin(K["unit"].astype(str), TYPE_A + ["control"])
    np.savez(a.work / "thp1a_enc.npz", **{k: v[m] for k, v in S.items()}, **{k: v[m] for k, v in rest.items()})
    for s in (spec("T_thp1_IRF1", ["IRF1"], False, "thp1_enc.npz", BOTH), spec("C_thp1a_pool", TYPE_A, True, "thp1a_enc.npz", RIDGE)):
        (a.work / f"{s['name']}.yaml").write_text(yaml.safe_dump(s, sort_keys=False))
    (a.results / "stages_manifest.json").write_text(json.dumps({**info, "chain": CHAIN, "side_views": SIDE, "states": {k: man[k] for k in ("device", "torch", "python", "seconds", "cells_file_sha256")}}, indent=1) + "\n")
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
