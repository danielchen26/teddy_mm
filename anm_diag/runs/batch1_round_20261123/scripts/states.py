#!/usr/bin/env python3
"""The system's states (TEDDY-G 70M, gene-mean of the last layer) for the cells of batch 1, on the owner's machine.

  python states.py --repo TEDDY_MM --ckpt TEDDY_G_70M --cells batch1_cells.tsv --out OUT_DIR [--check 96] [--budget SECONDS]

Reads RNA counts only: no protein value is opened here. Steps, each resumable (a file that exists is not made again; with
--budget the script stops after that many seconds with status 3 and the same command continues it):

  counts   RNA counts of the listed cells from the two screens' RNA files (CSC matrices, read in column blocks), written as
           counts_<study>.npz (rows, cols, vals over all genes; ids; genes; shape)
  tokens   the official tokeniser, as anm_diag/adapters/perturb_cite.py: genes mapped to the checkpoint's vocabulary through the
           study's symbol map, counts per 10,000 divided by the gene medians, the 2,048 largest values in order
  states   the encoder's forward pass, batches of 16 sorted by length, as anm_diag/tools/layer_states.py --z-only

Check: `--check N` adds N cells per study that have stored states (data/processed/<study>_teddy_states.npz, for the co-culture
screen frangieh2021_coculture_teddy_states.npz) and compares. The largest absolute difference is printed and written to
states_manifest.json; above 1e-4 the script exits with status 2.
"""
import argparse
import csv
import hashlib
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import h5py
import numpy as np
import scipy.sparse as sp

PAD = 43810
FILES = {"frangieh2021": ("FrangiehIzar2021_RNA.h5ad", "cell_name", "gene_map_frangieh2021.csv", "frangieh2021_coculture_teddy_states.npz"),
         "papalexi2021": ("PapalexiSatija2021_eccite_RNA.h5ad", "_index", "gene_map_papalexi2021.csv", "papalexi2021_teddy_states.npz")}
t0 = time.time()
say = lambda *a: print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)


class Stop(Exception):
    pass


def counts(src, obs_key, ids, parts, deadline):
    """Rows `ids` of a CSC cells x genes matrix as COO triplets over all genes of the file. Column blocks are kept in
    `parts` so that a run that stops for time continues where it was."""
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
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--check", type=int, default=96)
    ap.add_argument("--budget", type=float, default=1e9)
    ap.add_argument("--maps", type=Path, help="folder with the symbol maps, if not the checkout's data/reference")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    deadline = t0 + a.budget
    rows = [l.rstrip("\n").split("\t") for l in a.cells.read_text().splitlines()[1:] if l.strip()]
    by = {}
    for study, _cond, _tg, cell, _role in rows:
        by.setdefault(study, []).append(cell)
    man = {"cells_file_sha256": hashlib.sha256(a.cells.read_bytes()).hexdigest(), "studies": {}}
    stored = {}
    for study, ids in by.items():
        z = np.load(a.repo / "data" / "processed" / FILES[study][3])
        zi = z["ids"].astype(str)
        chk = [i for i in zi if i not in set(ids)][:a.check]
        stored[study] = {i: z["z_last"][k] for k, i in enumerate(zi) if i in set(chk)}
        by[study] = list(ids) + chk
        man["studies"][study] = {"cells_of_the_batch": len(ids), "cells_of_the_check": len(chk)}
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
    import importlib.util
    import torch
    spec = importlib.util.spec_from_file_location("teddy_encoder", a.repo / "teddy_mm" / "teddy_encoder.py")
    te = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(te)
    voc = {line.rstrip("\n"): i for i, line in enumerate(open(a.ckpt / "vocab.txt"))}
    meds = te.load_gene_medians(a.repo / "data" / "reference" / "teddy_gene_medians.json")
    tok_all, len_all, ids_all, study_all = [], [], [], []
    tf = a.out / "tokens.npz"
    for study, ids in ({} if tf.exists() else by).items():
        c = np.load(a.out / f"counts_{study}.npz")
        assert list(c["ids"].astype(str)) == ids, "counts file is for another list of cells"
        X = sp.csr_matrix((c["vals"], (c["rows"].astype(np.int64), c["cols"].astype(np.int64))), shape=tuple(int(v) for v in c["shape"]), dtype=np.float32)
        ref = a.maps if a.maps and (a.maps / FILES[study][2]).exists() else a.repo / "data" / "reference"
        with open(ref / FILES[study][2]) as f:
            mp = {r["symbol"]: r["ensembl"] for r in csv.DictReader(f)}
        ens = np.array([mp.get(s, "") for s in c["genes"].astype(str)])
        tk = np.array([voc.get(e, -1) for e in ens])
        mapped = np.where(tk >= 0)[0]
        token_ids = tk[mapped].astype(np.int64)
        fac = np.array([1.0 / meds[e] if e in meds else 1.0 for e in ens[mapped]])
        for s in range(0, X.shape[0], 256):
            Xm = np.ascontiguousarray(X[s:s + 256][:, mapped].toarray(), dtype=np.float32)
            t, att = te.rank_encode_official(te.official_values(Xm, fac), token_ids, pad_id=PAD)
            tok_all.append(t.astype(np.int64)); len_all.append(att.sum(1))
        ids_all += ids; study_all += [study] * len(ids)
        man["studies"][study].update(genes_mapped=int(len(mapped)), genes=int(len(ens)))
        say("tokens", study, len(ids), "cells; genes mapped", len(mapped), "of", len(ens))
    if not tf.exists():
        np.savez_compressed(a.out / "tmp_tokens.npz", tok=np.concatenate(tok_all).astype(np.int32), L=np.concatenate(len_all).astype(np.int32), ids=np.array(ids_all), study=np.array(study_all),
                            mapped=json.dumps({k: {q: v[q] for q in ("genes_mapped", "genes")} for k, v in man["studies"].items()}))
        (a.out / "tmp_tokens.npz").rename(tf)
    T = np.load(tf)
    tok, lens, ids_all, study_all = T["tok"].astype(np.int64), T["L"].astype(np.int64), T["ids"].astype(str), T["study"].astype(str)
    assert list(ids_all) == [i for st in by for i in by[st]], "tokens file is for another list of cells"
    for k, v in json.loads(str(T["mapped"])).items():
        man["studies"][k].update(v)
    # --- states
    dev = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    m = te.load_teddy(a.ckpt, dev)
    n = len(lens)
    Z, done = np.zeros((n, m.d_model), np.float32), np.zeros(n, bool)
    pf = a.out / "states_partial.npz"
    if pf.exists():
        q = np.load(pf); Z, done = q["Z"], q["done"]

    def keep_partial():
        np.savez(a.out / "tmp_states_partial.npz", Z=Z, done=done)
        (a.out / "tmp_states_partial.npz").rename(pf)
    order = np.argsort(-lens, kind="stable")
    with torch.inference_mode():
        for s in range(0, n, 16):
            idx = order[s:s + 16]; Lb = int(lens[idx].max())
            if done[idx].all():
                continue
            if time.time() > deadline:
                keep_partial()
                say(f"states: stopped for time, {int(done.sum())}/{n} done")
                sys.exit(3)
            t = torch.from_numpy(tok[idx, :Lb]).to(dev)
            am = (torch.arange(Lb)[None, :] < torch.from_numpy(lens[idx])[:, None]).long().to(dev)
            t = torch.where(am.bool(), t, torch.full_like(t, PAD))
            h, _ = m.hidden_states(t, am, return_layer_means=True)
            Z[idx] = te.pool_hidden(h.float(), am, "gene-mean").cpu().numpy()
            done[idx] = True
            if (s // 16) % 20 == 0:
                say(f"states {s + len(idx)}/{n} on {dev}")
    keep_partial()
    worst = 0.0
    for study in by:
        d = [float(np.abs(Z[k] - stored[study][i]).max()) for k, (i, st) in enumerate(zip(ids_all, study_all)) if st == study and i in stored[study]]
        man["studies"][study]["check_largest_difference"] = max(d) if d else None
        worst = max(worst, max(d) if d else 0.0)
        say("check", study, len(d), "cells, largest difference", max(d) if d else None)
    keep = np.array([i not in stored[st] for i, st in zip(ids_all, study_all)])
    np.savez_compressed(a.out / "batch1_states.npz", ids=ids_all[keep], study=study_all[keep], z_last=Z[keep], L=lens[keep].astype(np.int32))
    man.update(cells=int(keep.sum()), device=str(dev), torch=torch.__version__, numpy=np.__version__, python=sys.version.split()[0], seconds=round(time.time() - t0, 1),
               check_largest_difference=worst, check_passed=bool(worst <= 1e-4))
    (a.out / "states_manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    say("done:", int(keep.sum()), "cells; check", "passed" if worst <= 1e-4 else "FAILED", worst)
    if worst > 1e-4:
        sys.exit(2)


if __name__ == "__main__":
    try:
        main()
    except Stop as e:
        say(e)
        sys.exit(3)
