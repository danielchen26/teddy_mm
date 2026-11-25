#!/usr/bin/env python3
"""States of TEDDY-G 70M for the cells of batch 2 and the cells its readings stand on, on the owner's machine. RNA only.

  python states.py --repo TEDDY_MM --ckpt TEDDY_G_70M --cells batch2_cells.tsv --sets sets.json --out OUT_DIR [--budget SECONDS]

No protein value is opened here. Steps, each resumable (a file that exists is not made again; with --budget the script
stops after that many seconds with status 3 and the same command continues it):

  counts   RNA counts of the listed cells from the two screens' RNA files, written as counts_<study>.npz
  tokens   the official tokeniser, as runs/batch1_round_20261123/scripts/states.py
  states   the encoder as released, all twelve layers, batches of 32 sorted by length. For every cell: the mean of the
           last layer's output over its gene positions (the state), and for every set of genes of its study in sets.json
           the sum of the last layer's output over the positions of those genes, divided by the number of the cell's gene
           positions (a block of the typed state, as runs/typed_repair_20261125/scripts/states.py)

sets.json: per study, {name: [genes]}: the sets of run 40, chosen there on the fitted half of the development cells with
tokens and labels only (a gene is named by its symbol or by its identifier in the checkpoint's vocabulary).
Check: every listed cell that has a stored state (data/processed) is compared with it. The largest absolute difference
is written to states_manifest.json; above 1e-4, or with fewer than 96 cells compared in a study, the script exits with
status 2 and writes no state.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import h5py
import numpy as np
import scipy.sparse as sp

PAD, LIMIT, MIN_CHECK = 43810, 1e-4, 96
FILES = {"frangieh2021": ("FrangiehIzar2021_RNA.h5ad", "cell_name", "gene_map_frangieh2021.csv", "frangieh2021_coculture_teddy_states.npz"),
         "papalexi2021": ("PapalexiSatija2021_eccite_RNA.h5ad", "_index", "gene_map_papalexi2021.csv", "papalexi2021_teddy_states.npz")}
t0 = time.time()
say = lambda *a: print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)


class Stop(Exception):
    pass


def counts(src, obs_key, ids, parts, deadline):
    """Rows `ids` of a CSC cells x genes matrix as COO triplets over all genes of the file, in column blocks kept in `parts`."""
    f = h5py.File(src, "r")
    names = f["obs"][obs_key][:].astype(str)
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
        if time.time() > deadline:
            raise Stop(f"counts: stopped for time at column {c} of {ncol}")
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
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--sets", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--batch", type=int, default=32)
    ap.add_argument("--budget", type=float, default=1e9)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    deadline = t0 + a.budget
    rows = [l.rstrip("\n").split("\t") for l in a.cells.read_text().splitlines()[1:] if l.strip()]
    by = {}
    for study, _cond, _tg, cell, _role in rows:
        by.setdefault(study, []).append(cell)
    sets = json.load(open(a.sets))
    man = {"cells_file_sha256": hashlib.sha256(a.cells.read_bytes()).hexdigest(), "sets_file_sha256": hashlib.sha256(a.sets.read_bytes()).hexdigest(),
           "studies": {s: {"cells": len(v)} for s, v in by.items()}}
    # --- counts
    for study, ids in by.items():
        p = a.out / f"counts_{study}.npz"
        if p.exists():
            continue
        R, C, V, genes, ncol = counts(a.repo / "data" / "raw" / "perturb" / FILES[study][0], FILES[study][1], ids, a.out / f"parts_{study}", deadline)
        np.savez_compressed(a.out / f"tmp_counts_{study}.npz", rows=R, cols=C, vals=V, ids=np.array(ids), genes=genes, shape=np.array([len(ids), ncol]))
        (a.out / f"tmp_counts_{study}.npz").rename(p)
        say("counts", study, len(ids), "cells,", len(V), "entries")
    # --- tokens
    import torch
    spec = importlib.util.spec_from_file_location("teddy_encoder", a.repo / "teddy_mm" / "teddy_encoder.py")
    te = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(te)
    voc = {line.rstrip("\n"): i for i, line in enumerate(open(a.ckpt / "vocab.txt"))}
    meds = te.load_gene_medians(a.repo / "data" / "reference" / "teddy_gene_medians.json")
    maps = {}
    for study in by:
        with open(a.repo / "data" / "reference" / FILES[study][2]) as f:
            maps[study] = {r["symbol"]: r["ensembl"] for r in csv.DictReader(f)}
    tf = a.out / "tokens.npz"
    if not tf.exists():
        tok_all, len_all, ids_all, study_all, mapped_n = [], [], [], [], {}
        for study, ids in by.items():
            c = np.load(a.out / f"counts_{study}.npz")
            assert list(c["ids"].astype(str)) == ids, "counts file is for another list of cells"
            X = sp.csr_matrix((c["vals"], (c["rows"].astype(np.int64), c["cols"].astype(np.int64))), shape=tuple(int(v) for v in c["shape"]), dtype=np.float32)
            ens = np.array([maps[study].get(s, "") for s in c["genes"].astype(str)])
            tk = np.array([voc.get(e, -1) for e in ens])
            mapped = np.where(tk >= 0)[0]
            token_ids = tk[mapped].astype(np.int64)
            fac = np.array([1.0 / meds[e] if e in meds else 1.0 for e in ens[mapped]])
            for s in range(0, X.shape[0], 256):
                if time.time() > deadline:
                    raise Stop("tokens: stopped for time")
                Xm = np.ascontiguousarray(X[s:s + 256][:, mapped].toarray(), dtype=np.float32)
                t, att = te.rank_encode_official(te.official_values(Xm, fac), token_ids, pad_id=PAD)
                tok_all.append(t.astype(np.int32)); len_all.append(att.sum(1))
            ids_all += ids; study_all += [study] * len(ids)
            mapped_n[study] = {"genes_mapped": int(len(mapped)), "genes": int(len(ens))}
            say("tokens", study, len(ids), "cells; genes mapped", len(mapped), "of", len(ens))
        np.savez_compressed(a.out / "tmp_tokens.npz", tok=np.concatenate(tok_all).astype(np.int32), L=np.concatenate(len_all).astype(np.int32), ids=np.array(ids_all), study=np.array(study_all),
                            mapped=json.dumps(mapped_n))
        (a.out / "tmp_tokens.npz").rename(tf)
    T = np.load(tf)
    tok, lens, ids_all, study_all = T["tok"], T["L"].astype(np.int64), T["ids"].astype(str), T["study"].astype(str)
    assert list(ids_all) == [i for st in by for i in by[st]], "tokens file is for another list of cells"
    for k, v in json.loads(str(T["mapped"])).items():
        man["studies"][k].update(v)
    # --- the sets of genes as tokens
    dev = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    m = te.load_teddy(a.ckpt, dev)
    d, nvoc = m.d_model, m.embeddings.weight.shape[0]
    names = sorted({k for s in sets.values() for k in s})
    masks = {}
    for study in by:
        M = torch.zeros((len(names), nvoc), dtype=torch.bool)
        for j, k in enumerate(names):
            ids_k = [voc[g] if g in voc else voc[maps[study][g]] for g in sets[study][k]]
            assert len(set(ids_k)) == len(ids_k) == len(sets[study][k]), (study, k)
            M[j, torch.tensor(ids_k)] = True
        masks[study] = M.to(dev)
    # --- states
    n = len(lens)
    Z, B, done = np.zeros((n, d), np.float32), np.zeros((n, len(names), d), np.float32), np.zeros(n, bool)
    pf = a.out / "states_partial.npz"
    if pf.exists():
        q = np.load(pf); Z, B, done = q["Z"], q["B"], q["done"]

    def keep_partial():
        np.savez(a.out / "tmp_states_partial.npz", Z=Z, B=B, done=done)
        (a.out / "tmp_states_partial.npz").rename(pf)
    layers, last_kept = m.encoder.layers, time.time()
    with torch.no_grad():
        for study in by:
            rows_s = np.where(study_all == study)[0]
            order = rows_s[np.argsort(-lens[rows_s], kind="stable")]
            for s0 in range(0, len(order), a.batch):
                idx = np.sort(order[s0:s0 + a.batch])
                if done[idx].all():
                    continue
                if time.time() > deadline:
                    keep_partial()
                    say(f"states: stopped for time, {int(done.sum())}/{n} done")
                    sys.exit(3)
                Lm = int(lens[idx].max())
                t = torch.from_numpy(tok[idx, :Lm].astype(np.int64)).to(dev)
                am = (torch.arange(Lm)[None, :] < torch.from_numpy(lens[idx])[:, None]).to(dev)
                t = torch.where(am, t, torch.full_like(t, PAD))
                h = m.embeddings(t) + m.position_embeddings(torch.arange(Lm, device=dev))
                for layer in layers:
                    h = layer(src=h, src_key_padding_mask=~am)
                o = h.float() * am.unsqueeze(-1)
                N = am.sum(1, keepdim=True).float()
                Z[idx] = (o.sum(1) / N).cpu().numpy()
                mt = masks[study][:, t]                                                         # [sets, cells, positions]
                B[idx] = (torch.bmm(mt.float().permute(1, 0, 2), o) / N.unsqueeze(1)).cpu().numpy()
                done[idx] = True
                if (s0 // a.batch) % 20 == 0:
                    say(f"states {study} {min(s0 + a.batch, len(order))}/{len(order)} on {dev.type}")
                if time.time() - last_kept > 600:
                    keep_partial(); last_kept = time.time()
    keep_partial()
    worst, ok = 0.0, True
    for study in by:
        z = np.load(a.repo / "data" / "processed" / FILES[study][3])
        pos = {i: k for k, i in enumerate(z["ids"].astype(str))}
        rows_s = [k for k, (i, st) in enumerate(zip(ids_all, study_all)) if st == study and i in pos]
        diff = float(np.abs(Z[rows_s] - z["z_last"][np.array([pos[ids_all[k]] for k in rows_s])]).max()) if rows_s else None
        man["studies"][study].update(cells_compared_with_stored_states=len(rows_s), check_largest_difference=diff)
        say("check", study, len(rows_s), "cells, largest difference", diff)
        ok &= len(rows_s) >= MIN_CHECK and diff <= LIMIT
        worst = max(worst, diff or 0.0)
    man.update(cells=int(n), sets=names, device=dev.type, torch=torch.__version__, numpy=np.__version__, python=sys.version.split()[0], seconds=round(time.time() - t0, 1),
               check_largest_difference=worst, check_passed=bool(ok))
    (a.out / "states_manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    if not ok:
        say("check FAILED", worst)
        sys.exit(2)
    for study in by:
        k = study_all == study
        np.savez(a.out / f"batch2_states_{study}.npz", ids=ids_all[k], z_last=Z[k], L=lens[k].astype(np.int32), sets=np.array(names), blocks=B[k])
    say("done:", n, "cells; check passed", worst)


if __name__ == "__main__":
    try:
        main()
    except Stop as e:
        say(e)
        sys.exit(3)
