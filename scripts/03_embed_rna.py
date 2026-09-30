#!/usr/bin/env python3
"""Embed each cell's RNA with frozen TEDDY-G 70M.

--preprocessing official (default) follows the Merck/TEDDY code and tutorial:
  counts scaled to 10,000 per cell, each gene divided by its TEDDY gene median (1.0 when a gene
  has no median), no log1p; torch.topk over the vocab-mapped genes, k = 2048, zero values dropped,
  no CLS or annotation token, padded with <pad> to 2048; each encoder layer run in turn
  (model.py:256-258). Pooling: gene-mean (mean over real gene tokens, the paper's "averaging over
  output gene embeddings", default), first (token 0) or all-positions (the tutorial's mean over all
  2048 positions, padding included).
--preprocessing legacy reproduces teddy_mm's earlier embedding (data/processed/cite/z_rna.npy):
  counts scaled to 10,000, no gene-median step, stable argsort with ties in column order, 512
  tokens, nn.TransformerEncoder call, mean over real tokens. On MPS it needs
  PYTORCH_ENABLE_MPS_FALLBACK=1.

Work is split into shards of --shard-size cells written to <out-dir>/shards/; a rerun skips the
shards already written, and --time-budget-sec stops cleanly (exit code 75) so the caller can
resume. When every shard exists the script assembles z_rna.npy (+ extra poolings, token counts and
z_rna_manifest.json) in --out-dir and, for a full run, symlinks cite_arrays.npz and meta.json from
--processed so 04_train / 05_eval / the bridge export can point --processed at --out-dir.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from teddy_mm.data import load_prepared
from teddy_mm.device import resolve_device
from teddy_mm.teddy_encoder import (
    POOLINGS,
    load_gene_medians,
    load_pad_id,
    load_teddy,
    load_vocab,
    median_factors,
    official_values,
    pool_hidden,
    rank_encode_matrix,
    rank_encode_official,
)

EXIT_INCOMPLETE = 75
POOL_FILE = {"gene-mean": "genemean", "first": "first", "all-positions": "allpos"}
DEFAULT_MEDIANS = ROOT / "data/reference/teddy_gene_medians.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit() -> str | None:
    try:
        out = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--untracked-files=no"],
                               capture_output=True, text=True, check=True).stdout.strip()
        return out + ("+dirty" if dirty else "")
    except Exception:
        return None


def atomic_savez(path: Path, **arrays) -> None:
    tmp = path.with_name(path.name + ".tmp.npz")
    np.savez(tmp, **arrays)
    os.replace(tmp, path)


def atomic_save_npy(path: Path, arr: np.ndarray) -> None:
    tmp = path.with_name(path.name + ".tmp.npy")
    np.save(tmp, arr)
    os.replace(tmp, path)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite", help="dir with cite_arrays.npz + meta.json")
    p.add_argument("--out-dir", type=Path, default=None,
                   help="where z_rna.npy goes (default: --processed); an existing z_rna.npy is never overwritten without --overwrite")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--ckpt", type=Path, default=ROOT.parent / "teddy_mwe/ckpt/teddy_g_70M")
    p.add_argument("--preprocessing", choices=("official", "legacy"), default="official")
    p.add_argument("--medians", type=Path, default=None,
                   help="TEDDY gene-median JSON (official only; default data/reference/teddy_gene_medians.json, see scripts/01b_fetch_teddy_medians.sh)")
    p.add_argument("--seq-len", type=int, default=None,
                   help="tokens per cell (default 2048 official, 512 legacy; the stored legacy z used 512, the old script default was 1024)")
    p.add_argument("--pooling", choices=POOLINGS, default="gene-mean", help="pooling written to z_rna.npy (official only)")
    p.add_argument("--save-poolings", default="", help="comma list of extra poolings saved from the same forward pass (official only)")
    p.add_argument("--length-buckets", action="store_true",
                   help="sort each shard's cells by length and pad a batch only to its longest cell (official; not with all-positions)")
    p.add_argument("--save-layer-means", action="store_true", help="also save the gene-mean of every layer (float16; official only)")
    p.add_argument("--batch-size", type=int, default=16)
    p.add_argument("--autocast", choices=("none", "fp16", "bf16"), default="none",
                   help="official only: run the layers under torch.autocast in this dtype; states are cast to float32 before pooling")
    p.add_argument("--normalize-total", type=float, default=10000.0)
    p.add_argument("--device", default="auto")
    p.add_argument("--threads", type=int, default=6, help="torch CPU threads")
    p.add_argument("--max-cells", type=int, default=0, help="random subset of this many cells (pilots); 0 = all")
    p.add_argument("--cell-seed", type=int, default=0)
    p.add_argument("--start", type=int, default=None, help="contiguous cell range start (checks)")
    p.add_argument("--stop", type=int, default=None, help="contiguous cell range stop (exclusive)")
    p.add_argument("--shard-size", type=int, default=5000)
    p.add_argument("--shard-dir", type=Path, default=None, help="where shards go (default <out-dir>/shards)")
    p.add_argument("--time-budget-sec", type=float, default=0.0,
                   help="after at least one new shard, stop before the next once this many seconds have passed (exit 75; 0 = none)")
    args = p.parse_args(argv)
    if args.seq_len is None:
        args.seq_len = 2048 if args.preprocessing == "official" else 512
    extra = [s for s in args.save_poolings.split(",") if s]
    for s in extra:
        if s not in POOLINGS:
            p.error(f"--save-poolings: unknown pooling {s!r}")
    args.poolings = [args.pooling] + [s for s in extra if s != args.pooling]
    if args.preprocessing == "legacy":
        if args.poolings != ["gene-mean"] or args.length_buckets or args.save_layer_means or args.autocast != "none":
            p.error("legacy mode is the old code path: gene-mean only, fixed padding, no layer means")
    if args.length_buckets and "all-positions" in args.poolings:
        p.error("all-positions needs every cell padded to --seq-len; drop --length-buckets")
    if args.medians is None:
        args.medians = DEFAULT_MEDIANS
    if args.out_dir is None:
        args.out_dir = args.processed
    return args


def select_cells(n_total: int, args) -> tuple[np.ndarray, str]:
    if args.max_cells and (args.start is not None or args.stop is not None):
        raise SystemExit("use either --max-cells or --start/--stop")
    if args.max_cells:
        rng = np.random.default_rng(args.cell_seed)
        idx = np.sort(rng.choice(n_total, size=min(args.max_cells, n_total), replace=False))
        return idx, f"random{idx.size}_seed{args.cell_seed}"
    start = 0 if args.start is None else args.start
    stop = n_total if args.stop is None else min(args.stop, n_total)
    idx = np.arange(start, stop)
    return idx, ("all" if (start == 0 and stop == n_total) else f"range{start}_{stop}")


def main(argv=None) -> int:
    args = parse_args(argv)
    torch.set_num_threads(max(1, args.threads))
    device = resolve_device(args.device)
    t_start = time.time()
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    z_out = out / "z_rna.npy"
    if z_out.exists() and not args.overwrite:
        print(f"{z_out} exists; nothing to do (pass --overwrite to rebuild)")
        return 0

    pack = load_prepared(args.processed)
    rna = pack["rna"]
    token_ids = pack["token_ids"]
    n_total = rna.shape[0]
    cells, cells_tag = select_cells(n_total, args)
    vocab = load_vocab(args.ckpt)
    pad_id = load_pad_id(args.ckpt, vocab)

    medians_info = None
    factors = None
    if args.preprocessing == "official":
        if not args.medians.exists():
            raise SystemExit(f"missing gene medians {args.medians}; run scripts/01b_fetch_teddy_medians.sh")
        medians = load_gene_medians(args.medians)
        factors = median_factors(pack["rna_names"], medians)
        n_with = int(sum(str(g) in medians for g in pack["rna_names"]))
        medians_info = {"path": str(args.medians), "sha256": sha256_file(args.medians), "n_entries": len(medians),
                        "genes_with_median": n_with, "genes_total": int(len(pack["rna_names"])),
                        "missing_gene_factor": 1.0}
        print(f"medians: {n_with}/{len(pack['rna_names'])} genes have a median; the rest get factor 1.0")
    if args.preprocessing == "legacy" and device.type == "mps" and os.environ.get("PYTORCH_ENABLE_MPS_FALLBACK") != "1":
        print("WARN legacy nn.TransformerEncoder call on MPS usually needs PYTORCH_ENABLE_MPS_FALLBACK=1")

    # legacy shards must keep the old batch boundaries (batches of --batch-size from cell 0)
    shard = args.shard_size if args.shard_size > 0 else len(cells)
    shard = int(np.ceil(shard / args.batch_size) * args.batch_size)
    n_shards = int(np.ceil(len(cells) / shard)) if len(cells) else 0
    config = {
        "preprocessing": args.preprocessing, "seq_len": args.seq_len, "poolings": args.poolings,
        "length_buckets": bool(args.length_buckets), "save_layer_means": bool(args.save_layer_means),
        "batch_size": args.batch_size, "normalize_total": args.normalize_total, "cells": cells_tag,
        "autocast": args.autocast,
        "n_cells": int(len(cells)), "shard_size": shard,
        "medians_sha256": medians_info["sha256"] if medians_info else None,
    }
    cfg_hash = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()[:16]
    shard_dir = args.shard_dir or (out / "shards")
    shard_dir.mkdir(parents=True, exist_ok=True)
    print(f"device={device} preprocessing={args.preprocessing} seq_len={args.seq_len} poolings={args.poolings} "
          f"length_buckets={args.length_buckets} cells={len(cells)} shards={n_shards}x{shard} cfg={cfg_hash}", flush=True)

    model = None
    done_now = 0
    for k in range(n_shards):
        path = shard_dir / f"shard_{k:04d}.npz"
        if path.exists():
            with np.load(path) as old:
                if str(old["cfg_hash"]) != cfg_hash:
                    raise SystemExit(f"{path} was written with another configuration; use a new --out-dir")
            continue
        if args.time_budget_sec and done_now and time.time() - t_start > args.time_budget_sec:
            print(f"time budget reached before shard {k}/{n_shards}; rerun to resume", flush=True)
            return EXIT_INCOMPLETE
        if model is None:
            model = load_teddy(args.ckpt, device)
        t0 = time.time()
        idx = cells[k * shard: (k + 1) * shard]
        res = embed_shard(model, rna, idx, token_ids, pad_id, factors, args, device)
        atomic_savez(path, idx=idx, cfg_hash=np.array(cfg_hash), runtime_sec=np.array(time.time() - t0), **res)
        done_now += 1
        ms = 1000 * (time.time() - t0) / max(len(idx), 1)
        print(f"shard {k + 1}/{n_shards}: {len(idx)} cells in {time.time() - t0:.0f}s ({ms:.1f} ms/cell)", flush=True)

    # ---- assemble
    parts = [np.load(shard_dir / f"shard_{k:04d}.npz") for k in range(n_shards)]
    idx_all = np.concatenate([p["idx"] for p in parts])
    assert np.array_equal(idx_all, cells), "shard cell order mismatch"
    runtime = float(sum(float(p["runtime_sec"]) for p in parts))
    for how in args.poolings:
        arr = np.concatenate([p[f"z_{POOL_FILE[how]}"] for p in parts]).astype(np.float32)
        if how != args.pooling:
            atomic_save_npy(out / f"z_rna_{POOL_FILE[how]}.npy", arr)
    ntok = np.concatenate([p["ntokens"] for p in parts]).astype(np.int32)
    atomic_save_npy(out / "z_rna_ntokens.npy", ntok)
    if args.save_layer_means:
        lm = np.concatenate([p["layer_means"] for p in parts], axis=1).astype(np.float16)
        atomic_save_npy(out / "z_rna_layer_means.npy", lm)
    if cells_tag != "all":
        atomic_save_npy(out / "z_rna_cells.npy", cells.astype(np.int64))
    manifest = {
        **config,
        "config_hash": cfg_hash,
        "z_rna": str(z_out),
        "z_rna_definition": describe(args),
        "extra_poolings": {h: str(out / f"z_rna_{POOL_FILE[h]}.npy") for h in args.poolings if h != args.pooling},
        "pooling": args.pooling if args.preprocessing == "official" else "gene-mean",
        "ntokens": {"file": str(out / "z_rna_ntokens.npy"), "median": float(np.median(ntok)), "min": int(ntok.min()),
                    "max": int(ntok.max()), "share_at_seq_len": float(np.mean(ntok >= args.seq_len))},
        "medians": medians_info,
        "processed_input": str(args.processed),
        "ckpt": str(args.ckpt),
        "torch": torch.__version__, "numpy": np.__version__, "python": platform.python_version(),
        "device": str(device), "threads": args.threads,
        "mps_fallback_env": os.environ.get("PYTORCH_ENABLE_MPS_FALLBACK"),
        "git_commit": git_commit(),
        "embed_runtime_sec": runtime,
        "n_shards": n_shards,
        "shard_dir": str(shard_dir),
    }
    atomic_save_npy(z_out, np.concatenate([p[f"z_{POOL_FILE[args.pooling]}"] for p in parts]).astype(np.float32))
    (out / "z_rna_manifest.json").write_text(json.dumps(manifest, indent=2))
    if cells_tag == "all" and out.resolve() != args.processed.resolve():
        for name in ("cite_arrays.npz", "meta.json"):
            link = out / name
            if not link.exists():
                link.symlink_to((args.processed / name).resolve())
    print(f"saved {z_out} shape=({len(cells)}, {parts[0]['z_' + POOL_FILE[args.pooling]].shape[1]}) runtime {runtime:.0f}s", flush=True)
    return 0


def describe(args) -> str:
    if args.preprocessing == "legacy":
        return (f"legacy teddy_mm: counts/total*1e4 over the vocab-mapped genes, no gene-median step, stable argsort "
                f"(ties in column order), {args.seq_len} tokens, nn.TransformerEncoder, mean over real tokens")
    return (f"official TEDDY-G: counts/total*{args.normalize_total:g}, divided by TEDDY gene medians (1.0 if missing), "
            f"no log1p; torch.topk k={args.seq_len}, zeros dropped, no CLS/annotation tokens, pad <pad> to {args.seq_len}"
            f"{' (length-bucketed batches)' if args.length_buckets else ''}; layer loop with bool key-padding mask"
            f"{'' if args.autocast == 'none' else f', torch.autocast {args.autocast}'}; "
            f"z_rna.npy pooling = {args.pooling}")


def embed_shard(model, rna, idx, token_ids, pad_id, factors, args, device) -> dict:
    n = len(idx)
    bs = args.batch_size
    if args.preprocessing == "legacy":
        z = np.zeros((n, model.d_model), dtype=np.float32)
        ntok = np.zeros(n, dtype=np.int32)
        for s in range(0, n, bs):
            sl = slice(s, min(s + bs, n))
            block = rna[idx[sl]].toarray()
            tot = block.sum(axis=1, keepdims=True)
            tot[tot == 0] = 1.0
            block = block / tot * args.normalize_total
            tokens, attn = rank_encode_matrix(block, token_ids, args.seq_len, pad_id)
            with torch.no_grad():
                z[sl] = model(torch.from_numpy(tokens).to(device), torch.from_numpy(attn).to(device)).detach().cpu().numpy()
            ntok[sl] = attn.sum(1)
        return {"z_genemean": z, "ntokens": ntok}

    # official: tokenise the whole shard first (CPU, float32 torch.topk), then batch
    L = args.seq_len
    tokens = np.full((n, L), pad_id, dtype=np.int64)
    attn = np.zeros((n, L), dtype=np.int64)
    for s in range(0, n, 256):
        e = min(s + 256, n)
        vals = official_values(rna[idx[s:e]].toarray(), factors, args.normalize_total)
        tokens[s:e], attn[s:e] = rank_encode_official(vals, token_ids, max_len=L, pad_id=pad_id, pad_to=L)
    ntok = attn.sum(1).astype(np.int32)
    order = np.argsort(ntok, kind="stable") if args.length_buckets else np.arange(n)
    outs = {POOL_FILE[h]: np.zeros((n, model.d_model), dtype=np.float32) for h in args.poolings}
    amp = {"fp16": torch.float16, "bf16": torch.bfloat16}.get(args.autocast)
    layer_means = np.zeros((len(model.encoder.layers), n, model.d_model), dtype=np.float16) if args.save_layer_means else None
    with torch.no_grad():
        for s in range(0, n, bs):
            rows = order[s: s + bs]
            width = int(max(1, ntok[rows].max())) if args.length_buckets else L
            ids = torch.from_numpy(tokens[rows, :width]).to(device)
            mask = torch.from_numpy(attn[rows, :width]).to(device)
            with (torch.autocast(device_type=device.type, dtype=amp) if amp is not None else contextlib.nullcontext()):
                if args.save_layer_means:
                    h, lm = model.hidden_states(ids, mask, return_layer_means=True)
                    layer_means[:, rows] = lm.float().cpu().numpy().astype(np.float16)
                else:
                    h = model.hidden_states(ids, mask)
            h = h.float()
            for how in args.poolings:
                outs[POOL_FILE[how]][rows] = pool_hidden(h, mask, how).float().cpu().numpy()
    res = {f"z_{k}": v for k, v in outs.items()}
    res["ntokens"] = ntok
    if layer_means is not None:
        res["layer_means"] = layer_means
    return res


if __name__ == "__main__":
    sys.exit(main())
