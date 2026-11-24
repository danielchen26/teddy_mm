#!/usr/bin/env python3
"""The tokens of the cells of the two melanoma strata, on the owner's machine: what the encoder is given.

  python tokens.py --repo TEDDY_MM --ckpt TEDDY_G_70M --strata STRATA.npz --out OUT_DIR

STRATA.npz: for the strata mel_co and mel_ifn the cells (`<stratum>__ids`) and the state the earlier rounds read
(`<stratum>__Z`), in the order of those rounds. Reads RNA counts only: no protein value and no label of an event is
opened here. Steps, each resumable (a file that exists is not made again), as runs/encoder_stages_20261124/scripts/layers.py:

  counts   RNA counts of the listed cells from the screen's RNA file (Frangieh et al. 2021)
  tokens   the official tokeniser: genes mapped to the checkpoint's vocabulary through the study's symbol map, counts per
           10,000 divided by the gene medians, the 2,048 largest values in order

Writes OUT_DIR/tokens.npz (tok, L, ids) and OUT_DIR/recorded.npz (ids, L12: the recorded state of every cell, for the
check of runs/annotation_input_20261124/scripts/annot.py).
"""
import argparse
import csv
import importlib.util
import json
import os
import time
from pathlib import Path

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import h5py
import numpy as np
import scipy.sparse as sp

PAD = 43810
STUDY, RNA, OBS, MAP = "frangieh2021", "FrangiehIzar2021_RNA.h5ad", "cell_name", "gene_map_frangieh2021.csv"
STRATA = ["mel_co", "mel_ifn"]
t0 = time.time()
say = lambda *a: print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)


def column(obs, key):
    """A column of an AnnData obs group, plain or categorical."""
    x = obs[key]
    if isinstance(x, h5py.Group):
        return x["categories"][:].astype(str)[x["codes"][:]]
    return x[:].astype(str)


def counts(src, ids, parts):
    """Rows `ids` of a CSC cells x genes matrix as COO triplets over all genes of the file, read in column blocks."""
    f = h5py.File(src, "r")
    names = column(f["obs"], OBS)
    pos = {n: i for i, n in enumerate(names)}
    new = np.full(len(names), -1, np.int64)
    new[np.array([pos[i] for i in ids])] = np.arange(len(ids))
    X = f["X"]
    assert X.attrs["encoding-type"] == "csc_matrix"
    ip = X["indptr"][:].astype(np.int64)
    ncol = len(ip) - 1
    parts.mkdir(parents=True, exist_ok=True)
    blocks, c = [], 0
    while c < ncol:
        e = min(max(int(np.searchsorted(ip, ip[c] + 30_000_000, side="right")) - 1, c + 1), ncol)
        blocks.append((c, e)); c = e
    for c, e in blocks:
        p = parts / f"part_{c:05d}.npz"
        if p.exists():
            continue
        a, b = ip[c], ip[e]
        r = new[X["indices"][a:b]]
        keep = r >= 0
        np.savez(parts / f"tmp_{c:05d}.npz", rows=r[keep].astype(np.int32), cols=np.repeat(np.arange(c, e), np.diff(ip[c:e + 1]))[keep].astype(np.int32),
                 vals=X["data"][a:b][keep].astype(np.float32))
        (parts / f"tmp_{c:05d}.npz").rename(p)
    P = [np.load(parts / f"part_{c:05d}.npz") for c, _e in blocks]
    genes = f["var"][f["var"].attrs["_index"]][:].astype(str)
    return np.concatenate([q["rows"] for q in P]), np.concatenate([q["cols"] for q in P]), np.concatenate([q["vals"] for q in P]), genes, ncol


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--ckpt", type=Path, required=True)
    ap.add_argument("--strata", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    K = np.load(a.strata)
    ids = [i for g in STRATA for i in K[f"{g}__ids"].astype(str)]
    assert len(set(ids)) == len(ids)
    rf = a.out / "recorded.npz"
    if not rf.exists():
        np.savez(a.out / "tmp_recorded.npz", ids=np.array(ids), L12=np.concatenate([K[f"{g}__Z"] for g in STRATA]).astype(np.float32), stratum=np.array([g for g in STRATA for _ in K[f"{g}__ids"]]))
        (a.out / "tmp_recorded.npz").rename(rf)
    cf = a.out / f"counts_{STUDY}.npz"
    if not cf.exists():
        R, C, V, genes, ncol = counts(a.repo / "data" / "raw" / "perturb" / RNA, ids, a.out / f"parts_{STUDY}")
        np.savez_compressed(a.out / "tmp_counts.npz", rows=R, cols=C, vals=V, ids=np.array(ids), genes=genes, shape=np.array([len(ids), ncol]))
        (a.out / "tmp_counts.npz").rename(cf)
        say("counts", len(ids), "cells,", len(V), "entries")
    spec = importlib.util.spec_from_file_location("teddy_encoder", a.repo / "teddy_mm" / "teddy_encoder.py")
    te = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(te)
    tf = a.out / "tokens.npz"
    if not tf.exists():
        voc = {line.rstrip("\n"): i for i, line in enumerate(open(a.ckpt / "vocab.txt"))}
        meds = te.load_gene_medians(a.repo / "data" / "reference" / "teddy_gene_medians.json")
        c = np.load(cf)
        assert list(c["ids"].astype(str)) == ids, "counts file is for another list of cells"
        X = sp.csr_matrix((c["vals"], (c["rows"].astype(np.int64), c["cols"].astype(np.int64))), shape=tuple(int(v) for v in c["shape"]), dtype=np.float32)
        with open(a.repo / "data" / "reference" / MAP) as f:
            mp = {r["symbol"]: r["ensembl"] for r in csv.DictReader(f)}
        genes = c["genes"].astype(str)
        ens = np.array([mp.get(s, "") for s in genes])
        tk = np.array([voc.get(e, -1) for e in ens])
        mapped = np.where(tk >= 0)[0]
        token_ids = tk[mapped].astype(np.int64)
        fac = np.array([1.0 / meds[e] if e in meds else 1.0 for e in ens[mapped]])
        tok_all, len_all = [], []
        for s in range(0, X.shape[0], 256):
            Xm = np.ascontiguousarray(X[s:s + 256][:, mapped].toarray(), dtype=np.float32)
            t, att = te.rank_encode_official(te.official_values(Xm, fac), token_ids, pad_id=PAD)
            tok_all.append(t.astype(np.int64)); len_all.append(att.sum(1))
        np.savez_compressed(a.out / "tmp_tokens.npz", tok=np.concatenate(tok_all).astype(np.uint16), L=np.concatenate(len_all).astype(np.int32), ids=np.array(ids),
                            genes_in_file=np.array(len(genes)), genes_mapped=np.array(len(mapped)))
        (a.out / "tmp_tokens.npz").rename(tf)
        say("tokens", len(ids), "cells; genes mapped", len(mapped), "of", len(genes))
    T = np.load(tf)
    L = T["L"]
    man = {"study": STUDY, "cells": len(ids), "by_stratum": {g: int(len(K[f"{g}__ids"])) for g in STRATA}, "genes_in_file": int(T["genes_in_file"]), "genes_mapped": int(T["genes_mapped"]),
           "tokens_per_cell": {"min": int(L.min()), "median": float(np.median(L)), "max": int(L.max())}, "share_of_cells_at_2048": float((L == 2048).mean()), "seconds": round(time.time() - t0, 1)}
    (a.out / "tokens_manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    say("done:", json.dumps(man))


if __name__ == "__main__":
    main()
