#!/usr/bin/env python3
"""What enters the last layer of TEDDY-G 70M, for every cell of the THP-1 stratum, on the owner's machine.

  python cache.py --repo TEDDY_MM --ckpt TEDDY_G_70M --tokens TOKENS.npz --recorded THP1_ENC.npz --out CACHE [--batch 16]

The trial trains the last layer only. Everything before it is frozen, so its output is computed once and kept:
CACHE/h.npy, the output of the layer before the last at every position, as 16-bit floats [cells, 2048, d] (positions
after a cell's last token are zero). Reads tokens only; no protein value and no label of an event is opened.

Checks, written to CACHE/manifest.json (the script stops if one fails):
  the mean over gene tokens of what is kept equals the recorded state of that layer (THP1_ENC.npz, L11)      limit 2e-3
  the last layer, as released, applied to what is kept, gives the recorded state Z (L12)                      limit 2e-3
"""
import argparse
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import numpy as np

PAD, LMAX, LIMIT = 43810, 2048, 2e-3
t0 = time.time()
say = lambda *a: print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--ckpt", type=Path, required=True)
    ap.add_argument("--tokens", type=Path, required=True)
    ap.add_argument("--recorded", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--batch", type=int, default=16)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    import torch
    spec = importlib.util.spec_from_file_location("teddy_encoder", a.repo / "teddy_mm" / "teddy_encoder.py")
    te = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(te)
    T = np.load(a.tokens)
    tok, lens, ids = T["tok"].astype(np.int64), T["L"].astype(np.int64), T["ids"].astype(str)
    R = np.load(a.recorded)
    assert np.array_equal(R["ids"].astype(str), ids), "the recorded states are of other cells, or in another order"
    dev = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    m = te.load_teddy(a.ckpt, dev)
    n, d, lmax = len(ids), m.d_model, min(LMAX, tok.shape[1])
    f = a.out / "h.npy"
    H = np.lib.format.open_memmap(a.out / "tmp_h.npy", mode="w+", dtype=np.float16, shape=(n, lmax, d))
    worst = {"layer_before_the_last": 0.0, "last_layer_from_what_is_kept": 0.0, "last_layer_from_full_precision": 0.0}
    order = np.argsort(-lens, kind="stable")
    with torch.no_grad():
        for s in range(0, n, a.batch):
            idx = order[s:s + a.batch]
            Lb = int(lens[idx].max())
            t = torch.from_numpy(tok[idx, :Lb]).to(dev)
            am = (torch.arange(Lb)[None, :] < torch.from_numpy(lens[idx])[:, None]).long().to(dev)
            t = torch.where(am.bool(), t, torch.full_like(t, PAD))
            pad = ~am.bool()
            h = m.embeddings(t) + m.position_embeddings(torch.arange(Lb, device=dev))
            for layer in m.encoder.layers[:-1]:
                h = layer(src=h, src_key_padding_mask=pad)
            h = h.float() * am.unsqueeze(-1)
            kept = h.to(torch.float16)
            H[idx, :Lb] = kept.cpu().numpy()
            last = m.encoder.layers[-1]
            z_full = te.pool_hidden(last(src=h, src_key_padding_mask=pad).float(), am, "gene-mean").cpu().numpy()
            z_kept = te.pool_hidden(last(src=kept.float(), src_key_padding_mask=pad).float(), am, "gene-mean").cpu().numpy()
            worst["layer_before_the_last"] = max(worst["layer_before_the_last"], float(np.abs(te.pool_hidden(kept.float(), am, "gene-mean").cpu().numpy() - R["L11"][idx]).max()))
            worst["last_layer_from_what_is_kept"] = max(worst["last_layer_from_what_is_kept"], float(np.abs(z_kept - R["L12"][idx]).max()))
            worst["last_layer_from_full_precision"] = max(worst["last_layer_from_full_precision"], float(np.abs(z_full - R["L12"][idx]).max()))
            if (s // a.batch) % 40 == 0:
                say(f"{min(s + a.batch, n)}/{n} on {dev.type}; largest differences so far {json.dumps({k: round(v, 6) for k, v in worst.items()})}")
    H.flush()
    del H
    hold = bool(worst["layer_before_the_last"] <= LIMIT and worst["last_layer_from_what_is_kept"] <= LIMIT)
    man = {"cells": int(n), "positions": int(lmax), "d_model": int(d), "bytes": int(n * lmax * d * 2), "dtype": "float16", "layers_frozen": len(m.encoder.layers) - 1,
           "largest_difference_from_the_recorded_states": worst, "limit": LIMIT, "hold": hold, "device": dev.type, "torch": torch.__version__, "numpy": np.__version__,
           "python": sys.version.split()[0], "seconds": round(time.time() - t0, 1)}
    (a.out / "manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    if not hold:
        raise SystemExit("what is kept does not give the recorded states: stopped (the file tmp_h.npy is left for inspection)")
    (a.out / "tmp_h.npy").rename(f)
    np.savez(a.out / "index.npz", ids=ids, L=lens)
    say("done:", json.dumps(man))


if __name__ == "__main__":
    main()
