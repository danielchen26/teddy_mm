#!/usr/bin/env python3
"""What the encoder holds at every stage, for the cells of one stratum, on the owner's machine (TEDDY-G 70M).

  python layers.py --repo TEDDY_MM --ckpt TEDDY_G_70M --cells cells.tsv --out OUT_DIR

Reads RNA counts only: no protein value and no label of an event is opened here. cells.tsv: one cell name per line under
the header `cell` (THP-1 screen, Papalexi et al. 2021). Steps, each resumable (a file that exists is not made again):

  counts   RNA counts of the listed cells from the screen's RNA file (as runs/batch1_round_20261123/scripts/states.py)
  tokens   the official tokeniser: genes mapped to the checkpoint's vocabulary through the study's symbol map, counts per
           10,000 divided by the gene medians, the 2,048 largest values in order
  stages   one forward pass, batches of 16 sorted by length. Kept for every cell:
             E0        mean over its tokens of what enters layer 1 (gene embedding + position embedding)
             L1..L12   mean over its tokens of the output of every layer (L12 is the state Z the head reads)
             W0, W12   the same two places as eight means, one per window of 256 positions (ranks 1-256, 257-512, ...)

Check: L12 against the stored states (data/processed/papalexi2021_teddy_states.npz). The largest absolute difference is
printed and written to manifest.json; above 1e-4 the script exits with status 2.
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

PAD, WIN = 43810, 256
STUDY, RNA, OBS, MAP, STORED = "papalexi2021", "PapalexiSatija2021_eccite_RNA.h5ad", "_index", "gene_map_papalexi2021.csv", "papalexi2021_teddy_states.npz"
t0 = time.time()
say = lambda *a: print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)


def counts(src, ids, parts):
    """Rows `ids` of a CSC cells x genes matrix as COO triplets over all genes of the file, read in column blocks."""
    f = h5py.File(src, "r")
    names = f["obs"][OBS][:].astype(str)
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
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--maps", type=Path, help="folder with the symbol map, if not the checkout's data/reference")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    ids = [l.strip() for l in a.cells.read_text().splitlines()[1:] if l.strip()]
    man = {"cells_file_sha256": hashlib.sha256(a.cells.read_bytes()).hexdigest(), "study": STUDY, "cells": len(ids)}
    # --- counts
    cf = a.out / f"counts_{STUDY}.npz"
    if not cf.exists():
        R, C, V, genes, ncol = counts(a.repo / "data" / "raw" / "perturb" / RNA, ids, a.out / f"parts_{STUDY}")
        np.savez_compressed(a.out / "tmp_counts.npz", rows=R, cols=C, vals=V, ids=np.array(ids), genes=genes, shape=np.array([len(ids), ncol]))
        (a.out / "tmp_counts.npz").rename(cf)
        say("counts", len(ids), "cells,", len(V), "entries")
    # --- tokens
    import torch
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
        ref = a.maps if a.maps and (a.maps / MAP).exists() else a.repo / "data" / "reference"
        with open(ref / MAP) as f:
            mp = {r["symbol"]: r["ensembl"] for r in csv.DictReader(f)}
        genes = c["genes"].astype(str)
        ens = np.array([mp.get(s, "") for s in genes])
        tk = np.array([voc.get(e, -1) for e in ens])
        mapped = np.where(tk >= 0)[0]
        token_ids = tk[mapped].astype(np.int64)
        fac = np.array([1.0 / meds[e] if e in meds else 1.0 for e in ens[mapped]])
        tok_all, len_all, nz_all, nzm_all = [], [], [], []
        for s in range(0, X.shape[0], 256):
            B = X[s:s + 256]
            Xm = np.ascontiguousarray(B[:, mapped].toarray(), dtype=np.float32)
            t, att = te.rank_encode_official(te.official_values(Xm, fac), token_ids, pad_id=PAD)
            tok_all.append(t.astype(np.int64)); len_all.append(att.sum(1))
            nz_all.append(np.diff(B.indptr)); nzm_all.append((Xm > 0).sum(1))
        np.savez_compressed(a.out / "tmp_tokens.npz", tok=np.concatenate(tok_all).astype(np.uint16), L=np.concatenate(len_all).astype(np.int32), ids=np.array(ids),
                            map_symbols=genes[mapped], map_token_ids=token_ids.astype(np.int32), genes_detected=np.concatenate(nz_all).astype(np.int32),
                            genes_detected_in_vocabulary=np.concatenate(nzm_all).astype(np.int32), genes_in_file=np.array(len(genes)), genes_mapped=np.array(len(mapped)))
        (a.out / "tmp_tokens.npz").rename(tf)
        say("tokens", len(ids), "cells; genes mapped", len(mapped), "of", len(genes))
    T = np.load(tf)
    tok, lens = T["tok"].astype(np.int64), T["L"].astype(np.int64)
    assert list(T["ids"].astype(str)) == ids, "tokens file is for another list of cells"
    man.update(genes_in_file=int(T["genes_in_file"]), genes_mapped=int(T["genes_mapped"]), tokens_per_cell={"min": int(lens.min()), "median": float(np.median(lens)), "max": int(lens.max())},
               share_of_cells_at_2048=float((lens == 2048).mean()), genes_detected_per_cell_median=float(np.median(T["genes_detected"])),
               genes_detected_in_vocabulary_per_cell_median=float(np.median(T["genes_detected_in_vocabulary"])))
    # --- stages
    dev = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    m = te.load_teddy(a.ckpt, dev)
    n, d, nl = len(lens), m.d_model, len(m.encoder.layers)
    assert nl == 12
    E0, LM = np.zeros((n, d), np.float32), np.zeros((nl, n, d), np.float32)
    W0, W12, WC = np.zeros((n, 8 * d), np.float32), np.zeros((n, 8 * d), np.float32), np.zeros((n, 8), np.int32)

    def windows(h, am):
        out, cnt = [], []
        for w in range(8):
            lo, hi = w * WIN, min((w + 1) * WIN, h.size(1))
            if lo >= h.size(1):
                out.append(torch.zeros(h.size(0), h.size(2), device=h.device)); cnt.append(torch.zeros(h.size(0), device=h.device))
                continue
            mk = am[:, lo:hi].unsqueeze(-1).to(h.dtype)
            out.append((h[:, lo:hi] * mk).sum(1) / mk.sum(1).clamp(min=1.0)); cnt.append(mk.sum(1).squeeze(-1))
        return torch.cat(out, 1), torch.stack(cnt, 1)

    order = np.argsort(-lens, kind="stable")
    with torch.inference_mode():
        for s in range(0, n, 16):
            idx = order[s:s + 16]; Lb = int(lens[idx].max())
            t = torch.from_numpy(tok[idx, :Lb]).to(dev)
            am = (torch.arange(Lb)[None, :] < torch.from_numpy(lens[idx])[:, None]).long().to(dev)
            t = torch.where(am.bool(), t, torch.full_like(t, PAD))
            pad = ~am.bool()
            h = m.embeddings(t) + m.position_embeddings(torch.arange(Lb, device=dev))
            E0[idx] = te.pool_hidden(h.float(), am, "gene-mean").cpu().numpy()
            w, c = windows(h.float(), am)
            W0[idx], WC[idx] = w.cpu().numpy(), c.cpu().numpy().astype(np.int32)
            for k, layer in enumerate(m.encoder.layers):
                h = layer(src=h, src_key_padding_mask=pad)
                LM[k][idx] = te.pool_hidden(h.float(), am, "gene-mean").cpu().numpy()
            W12[idx] = windows(h.float(), am)[0].cpu().numpy()
            if (s // 16) % 40 == 0:
                say(f"stages {s + len(idx)}/{n} on {dev}")
    z = np.load(a.repo / "data" / "processed" / STORED)
    zp = {i: k for k, i in enumerate(z["ids"].astype(str))}
    have = np.array([i in zp for i in ids])
    worst = float(np.abs(LM[11][have] - z["z_last"][[zp[i] for i in np.array(ids)[have]]]).max()) if have.any() else None
    I = np.array(ids)
    np.savez(a.out / "layers_a.npz", ids=I, L=lens.astype(np.int32), E0=E0, **{f"L{k + 1}": LM[k] for k in range(0, 6)})
    np.savez(a.out / "layers_b.npz", ids=I, **{f"L{k + 1}": LM[k] for k in range(6, 12)})
    np.savez(a.out / "windows.npz", ids=I, W0=W0, W12=W12, window_counts=WC)
    man.update(device=str(dev), torch=torch.__version__, numpy=np.__version__, python=sys.version.split()[0], seconds=round(time.time() - t0, 1), cells_with_a_stored_state=int(have.sum()),
               check_largest_difference=worst, check_passed=bool(worst is not None and worst <= 1e-4), windows="eight windows of 256 positions; a window with no token is zero")
    (a.out / "manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    say("done:", n, "cells; check", "passed" if man["check_passed"] else "FAILED", worst)
    if not man["check_passed"]:
        sys.exit(2)


if __name__ == "__main__":
    main()
