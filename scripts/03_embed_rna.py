#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import torch
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from teddy_mm.data import load_prepared
from teddy_mm.device import resolve_device
from teddy_mm.teddy_encoder import load_pad_id, load_teddy, load_vocab, rank_encode_matrix


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--ckpt", type=Path, default=ROOT.parent / "teddy_mwe/ckpt/teddy_g_70M")
    p.add_argument("--seq-len", type=int, default=1024)
    p.add_argument("--batch-size", type=int, default=16)
    p.add_argument("--normalize-total", type=float, default=10000.0)
    p.add_argument("--device", default="auto")
    args = p.parse_args()

    device = resolve_device(args.device)
    print("device", device)
    pack = load_prepared(args.processed)
    vocab = load_vocab(args.ckpt)
    pad_id = load_pad_id(args.ckpt, vocab)
    model = load_teddy(args.ckpt, device)

    rna = pack["rna"]
    n = rna.shape[0]
    z_all = np.zeros((n, model.d_model), dtype=np.float32)
    token_ids = pack["token_ids"]

    for start in tqdm(range(0, n, args.batch_size), desc="TEDDY embed"):
        sl = slice(start, min(start + args.batch_size, n))
        block = rna[sl].toarray()
        tot = block.sum(axis=1, keepdims=True)
        tot[tot == 0] = 1.0
        block = block / tot * args.normalize_total
        tokens, attn = rank_encode_matrix(block, token_ids, args.seq_len, pad_id)
        gene_ids = torch.from_numpy(tokens).to(device)
        mask = torch.from_numpy(attn).to(device)
        with torch.no_grad():
            z = model(gene_ids, mask)
        z_all[sl] = z.detach().cpu().numpy()

    out = args.processed / "z_rna.npy"
    np.save(out, z_all)
    print(f"saved {out} shape={z_all.shape}")


if __name__ == "__main__":
    main()
