#!/usr/bin/env python3
"""Per-layer TEDDY states for the formal ANM test (runs on the Mac, teddy_mm .venv, MPS or CPU).

Input : a tokens file with arrays tok [n, 2048] (official TEDDY tokens, pad 43810) and L [n].
Output: <out>.npz with layer_means [12, n, 512] (gene-mean of each encoder layer's output) and
        z_last [n, 512] (gene-mean of the last layer = the embedding used in KS-2).
        --z-only keeps only z_last (with L, rows, ids): a small file (about 2 KB per cell) that can be
        copied back to the cloud workspace; use it when the declaration names no layer or window state.
Check : extra npz arguments (rows/z, plus the input's 'rows') report max |z_last - z_ref|, which confirms
        that the official tokens reproduce the stored KS-2 embeddings.

Usage (on the Mac; tokens are recomputed there so they match the stored embeddings exactly):
  cd <this folder>
  <teddy_mm repo>/.venv/bin/python layer_states.py counts:rank_need.npy layers_dev.npz \
      <teddy_mm repo>/ks2/pipe/out/z_V.npz <teddy_mm repo>/ks2/pipe/out/z_controls.npz
  rank_need.npy holds row indices into ks2/pipe/out/counts.npz (the 2,716 development cells).
"""
import os, sys, time
from pathlib import Path

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import numpy as np
import torch

ROOT = Path(os.environ.get("TEDDY_ROOT", Path(__file__).resolve().parents[2]))   # the teddy_mm checkout
sys.path.insert(0, str(ROOT))
from teddy_mm.teddy_encoder import load_teddy, pool_hidden  # noqa: E402

CKPT = Path(os.environ.get("TEDDY_CKPT", (ROOT / "../teddy_mwe/ckpt/teddy_g_70M").resolve()))
PAD = 43810
WINDOWS = [(0, 64), (64, 256), (256, 2048)]


def tokens_on_this_machine(cell_idx_file: str) -> dict:
    """Official tokens computed here with KS-2's own GeneMap (same torch as the stored embeddings)."""
    import csv, gzip
    from scipy import sparse
    sys.path.insert(0, str(ROOT / "ks2/pipe"))
    from common import GeneMap  # noqa: E402  (KS-2 pipeline: checkpoint vocab + TEDDY medians)
    out_dir = ROOT / "ks2/pipe/out"
    c = np.load(out_dir / "counts.npz")
    X = sparse.csr_matrix((c["data"], c["indices"], c["indptr"]), shape=tuple(c["shape"]))
    rows_all = np.load(out_dir / "count_rows.npy")
    with gzip.open(ROOT / "ks2/meta/var.csv.gz", "rt") as f:
        r = csv.reader(f); h = next(r); ens = [row[h.index("ensembl_id")] for row in r]
    gm = GeneMap(np.array(ens))
    idx = np.load(cell_idx_file)
    toks = gm.tokens(X[idx].toarray().astype(np.float32))
    tok = np.full((len(toks), 2048), PAD, np.int64); L = np.array([len(t) for t in toks])
    for i, t in enumerate(toks):
        tok[i, :len(t)] = t
    return {"tok": tok, "L": L, "rows": rows_all[idx]}


def main(inp: str, out: str, z_only: bool = False, zrefs=()) -> None:
    d = tokens_on_this_machine(inp[len("counts:"):]) if inp.startswith("counts:") else dict(np.load(inp))
    if "rows" in d:
        d["rows"] = np.asarray(d["rows"])
    tok = d["tok"].astype(np.int64); lens = d["L"].astype(np.int64)
    dev = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    m = load_teddy(CKPT, dev)
    n = len(lens); nl = len(m.encoder.layers)
    Z = np.zeros((n, m.d_model), np.float32)
    LM = None if z_only else np.zeros((nl, n, m.d_model), np.float32)
    WM = None if z_only else np.zeros((nl + 1, len(WINDOWS), n, m.d_model), np.float32)   # layer 0 = input embeddings
    order = np.argsort(-lens, kind="stable"); bs = 16; t0 = time.time()
    with torch.inference_mode():
        for s in range(0, n, bs):
            idx = order[s:s + bs]; Lb = int(lens[idx].max())
            t = torch.from_numpy(tok[idx, :Lb]).to(dev)
            a = (torch.arange(Lb)[None, :] < torch.from_numpy(lens[idx])[:, None]).long().to(dev)
            t = torch.where(a.bool(), t, torch.full_like(t, PAD))
            # walk the layers ourselves to also keep window means (positions 0-63, 64-255, 256+),
            # which separate what the encoder keeps from what whole-sequence mean pooling removes
            if z_only:
                h, _ = m.hidden_states(t, a, return_layer_means=True)
                Z[idx] = pool_hidden(h.float(), a, "gene-mean").cpu().numpy()
                if (s // bs) % 20 == 0:
                    print(f"{s + len(idx)}/{n}  {time.time() - t0:.0f}s on {dev}", flush=True)
                continue
            pos = torch.arange(Lb, device=dev)
            hcur = m.embeddings(t) + m.position_embeddings(pos)
            pad_mask = ~a.bool()
            for li, layer in enumerate([None] + list(m.encoder.layers)):
                if layer is not None:
                    hcur = layer(src=hcur, src_key_padding_mask=pad_mask)
                af = a.to(hcur.dtype).unsqueeze(-1)
                for wi, (s0, s1) in enumerate(WINDOWS):
                    seg = hcur[:, s0:s1] * af[:, s0:s1]
                    WM[li, wi, idx] = (seg.sum(1) / af[:, s0:s1].sum(1).clamp(min=1.0)).float().cpu().numpy()
                if layer is not None:
                    LM[li - 1, idx] = pool_hidden(hcur, a, "gene-mean").float().cpu().numpy()
            h, lm = m.hidden_states(t, a, return_layer_means=True)
            Z[idx] = pool_hidden(h.float(), a, "gene-mean").cpu().numpy()
            assert np.allclose(lm.float().cpu().numpy(), LM[:, idx], atol=1e-4), "layer walk differs from hidden_states"
            if (s // bs) % 20 == 0:
                print(f"{s + len(idx)}/{n}  {time.time() - t0:.0f}s on {dev}", flush=True)
    keep = {"rows": d["rows"] if "rows" in d else np.arange(n), **({"ids": d["ids"]} if "ids" in d else {})}
    if z_only:
        np.savez_compressed(out, z_last=Z, L=lens, **keep)
    else:
        np.savez_compressed(out, layer_means=LM, window_means=WM, windows=np.array(WINDOWS), z_last=Z,
                            tok=tok, L=lens, **keep)
    if zrefs:   # optional check against stored KS-2 embeddings (npz files with rows, z)
        zr = {}
        for f in zrefs:
            e = np.load(f); zr.update({int(r): v for r, v in zip(e["rows"], e["z"])})
        rows = d["rows"]; ok = [i for i in range(n) if int(rows[i]) in zr]
        diff = max(float(np.abs(Z[i] - zr[int(rows[i])]).max()) for i in ok) if ok else float("nan")
        print(f"check vs stored KS-2 embeddings: {len(ok)} cells, max |diff| = {diff:.2e}", flush=True)
    print(f"done: {n} cells, {nl} layers, {time.time() - t0:.0f}s -> {out}", flush=True)


if __name__ == "__main__":
    pos_args = [x for x in sys.argv[1:] if not x.startswith("--")]
    main(pos_args[0], pos_args[1], z_only="--z-only" in sys.argv[1:], zrefs=pos_args[2:])
