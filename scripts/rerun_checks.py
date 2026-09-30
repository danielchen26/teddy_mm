#!/usr/bin/env python3
"""Checks for the official-preprocessing rerun (scripts/rerun_official.sh).

legacy: a legacy-mode re-embedding (03_embed_rna.py --preprocessing legacy) against the stored
        z_rna.npy: bit-identical rows, max |diff|, cosine. Fails if the minimum cosine is <= 0.9999.
pilot:  a random-cell official pilot (fixed padding, three poolings) against the legacy z,
        descriptive only: cosine per cell by split, 10-NN overlap inside the pilot, the share of
        each cell's expressed genes left outside the context (512 legacy vs 2048 official), and a
        length-bucketing check (gene-mean and token 0 with batches padded to their longest cell vs
        padded to 2048) on 64 pilot cells.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def cos_rows(a, b):
    return (a * b).sum(1) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1) + 1e-12)


def q(x):
    x = np.asarray(x, dtype=np.float64)
    return {"median": float(np.median(x)), "p05": float(np.quantile(x, 0.05)), "p95": float(np.quantile(x, 0.95)),
            "min": float(x.min()), "mean": float(x.mean())}


def knn(z, k):
    zn = z / np.linalg.norm(z, axis=1, keepdims=True)
    S = zn @ zn.T
    np.fill_diagonal(S, -np.inf)
    return np.argpartition(-S, kth=k - 1, axis=1)[:, :k]


def legacy(args):
    d, L = Path(args.run), Path(args.legacy)
    man = json.loads((d / "z_rna_manifest.json").read_text())
    tag = man["cells"]
    assert tag.startswith("range"), tag
    a, b = map(int, tag[len("range"):].split("_"))
    idx = np.arange(a, b)
    z = np.load(d / "z_rna.npy")
    ref = np.asarray(np.load(L / "z_rna.npy", mmap_mode="r")[idx])
    c = cos_rows(z, ref)
    out = {"n_cells": int(idx.size), "cells": [a, b], "bit_identical_rows": int((z == ref).all(1).sum()),
           "max_abs_diff": float(np.abs(z - ref).max()), "cos_min": float(c.min()), "cos_median": float(np.median(c)),
           "device": man["device"], "torch": man["torch"], "mps_fallback_env": man.get("mps_fallback_env")}
    (d / "compare.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out))
    if out["cos_min"] <= 0.9999:
        raise SystemExit(f"legacy re-embedding does not reproduce the stored z (min cosine {out['cos_min']})")


def pilot(args):
    import torch

    from teddy_mm.data import load_prepared
    from teddy_mm.device import resolve_device
    from teddy_mm.teddy_encoder import (load_gene_medians, load_pad_id, load_teddy, load_vocab, median_factors,
                                        official_values, pool_hidden, rank_encode_official)

    d, L = Path(args.run), Path(args.legacy)
    man = json.loads((d / "z_rna_manifest.json").read_text())
    cells = np.load(d / "z_rna_cells.npy")
    zo = {"gene-mean": np.load(d / "z_rna.npy"), "first": np.load(d / "z_rna_first.npy"),
          "all-positions": np.load(d / "z_rna_allpos.npy")}
    zl = np.asarray(np.load(L / "z_rna.npy", mmap_mode="r")[cells])
    pack = load_prepared(L)
    split = pack["split"][cells]
    out = {"n_cells": int(cells.size), "cells": man["cells"], "manifest": man, "cosine_to_legacy": {}, "knn10_overlap_with_legacy": {}}
    for how, z in zo.items():
        c = cos_rows(z, zl)
        out["cosine_to_legacy"][how] = {"all": q(c), **{s: q(c[split == s]) for s in ("train", "val", "test")}}
    out["cosine_between_poolings"] = {"gene-mean_vs_first": q(cos_rows(zo["gene-mean"], zo["first"])),
                                      "gene-mean_vs_all-positions": q(cos_rows(zo["gene-mean"], zo["all-positions"]))}
    nl = knn(zl, 10)
    for how, z in zo.items():
        no = knn(z, 10)
        ov = np.array([len(set(a) & set(b)) / 10 for a, b in zip(nl, no)])
        out["knn10_overlap_with_legacy"][how] = q(ov)
    # share of each cell's expressed (non-zero) genes that fall outside the context
    nnz = np.diff(pack["rna"].indptr)[cells]
    ntok = np.load(d / "z_rna_ntokens.npy")
    drop_off = np.clip(nnz - ntok, 0, None) / np.maximum(nnz, 1)
    drop_leg = np.clip(nnz - 512, 0, None) / np.maximum(nnz, 1)
    out["context"] = {
        "expressed_genes": q(nnz), "official_tokens": q(ntok),
        "share_truncated_official_2048": float(np.mean(nnz > 2048)), "share_truncated_legacy_512": float(np.mean(nnz > 512)),
        "expressed_genes_outside_context": {
            "official_2048": {s: q(drop_off[split == s]) for s in ("train", "val", "test")},
            "legacy_512": {s: q(drop_leg[split == s]) for s in ("train", "val", "test")}},
        "note": "panel-gene dropout not computed: the processed data holds Ensembl ids only, no gene symbols",
    }
    # length-bucketing check on the first 64 pilot cells
    device = resolve_device(args.device)
    ck = ROOT.parent / "teddy_mwe/ckpt/teddy_g_70M" if args.ckpt is None else Path(args.ckpt)
    model = load_teddy(ck, device)
    pad = load_pad_id(ck, load_vocab(ck))
    fac = median_factors(pack["rna_names"], load_gene_medians(Path(args.medians)))
    sub = np.arange(64)
    vals = official_values(pack["rna"][cells[sub]].toarray(), fac)
    tok, att = rank_encode_official(vals, pack["token_ids"], max_len=2048, pad_id=pad)
    nt = att.sum(1); order = np.argsort(nt, kind="stable")
    gm = np.zeros((64, model.d_model), np.float32); f0 = np.zeros_like(gm)
    with torch.no_grad():
        for s in range(0, 64, 16):
            r = order[s:s + 16]; w = int(nt[r].max())
            ids = torch.from_numpy(tok[r, :w]).to(device); mk = torch.from_numpy(att[r, :w]).to(device)
            h = model.hidden_states(ids, mk)
            gm[r] = pool_hidden(h, mk, "gene-mean").cpu().numpy(); f0[r] = pool_hidden(h, mk, "first").cpu().numpy()
    out["length_bucket_check_64_cells"] = {
        "tokens_identical_to_pilot": bool(np.array_equal(nt, ntok[sub])),
        "gene-mean": {"max_abs_diff": float(np.abs(gm - zo["gene-mean"][sub]).max()), "cos_min": float(cos_rows(gm, zo["gene-mean"][sub]).min())},
        "first": {"max_abs_diff": float(np.abs(f0 - zo["first"][sub]).max()), "cos_min": float(cos_rows(f0, zo["first"][sub]).min())},
    }
    (d / "compare.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k != "manifest"}, indent=1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=("legacy", "pilot"))
    ap.add_argument("--run", required=True)
    ap.add_argument("--legacy", required=True, help="dir with the stored legacy z_rna.npy + cite_arrays.npz")
    ap.add_argument("--medians", default=str(ROOT / "data/reference/teddy_gene_medians.json"))
    ap.add_argument("--ckpt", default=None)
    ap.add_argument("--device", default="auto")
    args = ap.parse_args()
    {"legacy": legacy, "pilot": pilot}[args.what](args)


if __name__ == "__main__":
    main()
