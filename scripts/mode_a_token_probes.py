#!/usr/bin/env python3
"""Mode A, step 3 (token level): is the NK-T difference in TEDDY's per-token states even where the
gene-mean loses it?

mode_a_layer_probes.py (linear) and mode_a_nonlinear_probes.py (MLP, kNN) read the NK-T protein gap
off the gene-mean of each depth. The gene-mean averages ~1,300 token states; a difference carried by a
few gene tokens can be diluted there. This script re-embeds a small fixed set of site4 cells with the
official TEDDY-G preprocessing and keeps, for three depths only (input embedding, layer 6, layer 12):

  * the hidden state at the position of each of 13 relevant gene tokens, where the gene is among the
    cell's tokens (float16; a zero vector and presence flag 0 where it is not):
      4 coding genes of the probed proteins: NCAM1 (CD56), KLRD1 (CD94), NCR1 (CD335), CD3E (CD3);
      9 NK-vs-T genes: CD3D, CD3G, CD5, CD6, CD28 (T), KLRF1, FCGR3A, KLRC1, SH2D1B (NK);
    (Ensembl ids below, GRCh38. All 13 are in the processed genes and the TEDDY vocab; checked on this
    dataset: each is expressed in a much larger share of one class than the other, and the four coding
    genes correlate with their measured protein. The TCR constant genes TRAC / TRBC1 / TRBC2 / TRDC
    are NOT in the processed genes or the vocab, so they cannot be tokens. The script logs presence and
    writes the per-cell rank of each gene; --genes takes another set);
  * the gene-mean of that depth (same forward; checked against z_rna_layer_means.npy);
  * a CLS-free attention-pooled summary: TEDDY-G has no CLS token, so the summary is a parameter-free
    attention pooling, sum_t w_t h_t / sum_t w_t with w_t = attention that token t receives, averaged
    over heads and over the cell's real query tokens, taken from the layer's own self-attention
    (layer 6 and layer 12: that layer's attention; input: layer 1's attention over the input states).

Cells (<= 1,500): every cell of the site4 NK-T look-alike pairs without gdT CD158b+ cells (the
'no_gdT158' set of mode_a_layer_probes.py, rebuilt with the same code and checked against the linear
run), plus a random set of other site4 NK and T cells (no gdT CD158b+, not in any look-alike pair),
equal numbers of NK and T, stratified to the token-count quintiles of the pair cells of the same class.

Probes (--stage fit, CPU): because only site4 cells are embedded, the probes are trained on the random
cells (80 % fit / 20 % val, seeded) and scored on the pair cells, which none of them saw. That is a
within-site protocol, easier than the training-site -> site4 protocol of the gene-mean runs, so every
token feature set is compared with the gene-mean of the SAME forward pass under the SAME protocol:
feature sets genemean, attnpool, tokens (13 x 512 states + 13 presence flags), tokens+genemean;
read-outs ridge (alpha on val), MLP (2 hidden layers, early stopping on val, 3 seeds, seed-mean
prediction) and kNN (cosine on standardised features, k on val). Reported: Pearson on the pair cells and,
on the no_gdT158 pairs, the gap ratio median |pred NK - pred T| / median |meas NK - meas T| per protein
with the NK-cell bootstrap interval, plus paired differences (token set - genemean) per read-out.

Stages: --stage embed (needs the model; GPU/MPS recommended), --stage fit (CPU only), --stage all.
Outputs in --out-dir: token_states.npz, mode_a_token_results.json, README_token.md.

Command (from the repo root; not run yet):
  PY=<venv>/bin/python
  $PY scripts/mode_a_token_probes.py --stage all \
      --processed data/processed/cite --embed-dir data/processed/cite_official \
      --ckpt ../teddy_mwe/ckpt/teddy_g_70M --medians data/reference/teddy_gene_medians.json \
      --linear-dir outputs/mode_a_official --out-dir outputs/mode_a_official/token_probes \
      --device auto --autocast fp16 --batch-size 4 --threads 4
The embed stage prints its own time estimate before it starts (from the official run's measured
ms/cell and the selected cells' token counts).
"""
from __future__ import annotations

import os
import sys


def _early_threads(argv) -> int:
    n = 4
    for i, a in enumerate(argv):
        if a == "--threads" and i + 1 < len(argv):
            n = int(argv[i + 1])
        elif a.startswith("--threads="):
            n = int(a.split("=", 1)[1])
    return max(1, n)


_THREADS = _early_threads(sys.argv[1:])
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = str(_THREADS)

import argparse  # noqa: E402
import contextlib  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


lp = _load("mode_a_layer_probes")
nl = _load("mode_a_nonlinear_probes")

SCRIPT_VERSION = "mode_a_token_probes v1"
PROTEINS, GDT158 = lp.PROTEINS, lp.GDT158
rnd, say, pearson = lp.rnd, lp.say, lp.pearson
DEPTHS = ("input", "layer_6", "layer_12")          # which states are kept
ATTN_FROM = {"input": 1, "layer_6": 6, "layer_12": 12}  # whose self-attention gives the pooling weights
GENES = {  # symbol -> Ensembl (GRCh38)
    "NCAM1": "ENSG00000149294", "KLRD1": "ENSG00000134539", "NCR1": "ENSG00000189430", "CD3E": "ENSG00000198851",
    "CD3D": "ENSG00000167286", "CD3G": "ENSG00000160654", "CD5": "ENSG00000110448", "CD6": "ENSG00000013725",
    "CD28": "ENSG00000178562", "KLRF1": "ENSG00000150045", "FCGR3A": "ENSG00000203747", "KLRC1": "ENSG00000134545",
    "SH2D1B": "ENSG00000198574",
}
FEATURE_SETS = ("genemean", "attnpool", "tokens", "tokens+genemean")
# measured by the official run (z_rna_manifest.json: 10615 s for 90,261 cells, MPS, fp16, batch 32, length buckets)
OFFICIAL_SEC_PER_CELL = 10615.457 / 90261


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", choices=("embed", "fit", "all"), default="all")
    p.add_argument("--processed", type=Path, required=True)
    p.add_argument("--embed-dir", type=Path, required=True, help="official embedding (z_rna.npy for the pairs, layer means for the check)")
    p.add_argument("--ckpt", type=Path, default=ROOT.parent / "teddy_mwe/ckpt/teddy_g_70M")
    p.add_argument("--medians", type=Path, default=ROOT / "data/reference/teddy_gene_medians.json")
    p.add_argument("--linear-dir", type=Path, required=True, help="mode_a_layer_probes.py out-dir (pairs check)")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--genes", default=",".join(f"{s}={e}" for s, e in GENES.items()), help="SYMBOL=ENSG,... ")
    p.add_argument("--max-cells", type=int, default=1500)
    p.add_argument("--cell-seed", type=int, default=0)
    p.add_argument("--device", default="auto")
    p.add_argument("--autocast", choices=("none", "fp16", "bf16"), default="fp16", help="as the official run (fp16)")
    p.add_argument("--batch-size", type=int, default=4)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--limit", type=int, default=0, help="smoke only: embed just the first N selected cells")
    p.add_argument("--mlp-seeds", default="0,1,2")
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--smoke", default=None, metavar="NOTE")
    a = p.parse_args(argv)
    a.genes = [tuple(x.split("=", 1)) for x in a.genes.split(",") if x]
    a.mlp_seeds = [int(s) for s in a.mlp_seeds.split(",") if s != ""]
    if a.limit and not a.smoke:
        p.error("--limit is for smoke runs only")
    return a


# ============================================================================ cells

def select_cells(a, meta, emb):
    """Pair cells (no_gdT158) + stratified random site4 NK / T cells. Returns global indices + roles."""
    key = np.array([lp.ANNOTATION_MAP[c] or "" for c in meta["cell_types"]])
    ctype, split = meta["cell_types"], meta["split"]
    cells = emb.cells
    if emb.partial:
        raise SystemExit("needs the full official embedding (pairs are built on all site4 cells)")
    LR = json.loads((a.linear_dir / "mode_a_results.json").read_text())
    cfg = LR["config"]
    ev = np.where(split[cells] == "test")[0]
    z = np.asarray(emb.z, dtype=np.float32)
    pairs = lp.build_pairs(z[ev], ev, key[cells], ctype[cells], int(cfg["k"]), int(cfg["n_random_per_seed"]), list(cfg["random_seeds"]))
    LG = np.load(a.linear_dir / "mode_a_pair_gaps.npz")
    nk, t = pairs["neighbour"]["no_gdT158"]
    if not (np.array_equal(LG["nk_no_gdT158"], cells[nk]) and np.array_equal(LG["t_no_gdT158"], cells[t])):
        raise SystemExit("rebuilt no_gdT158 pairs differ from the linear run's")
    nk_all, t_all = pairs["neighbour"]["all"]
    in_any_pair = np.zeros(len(cells), bool)
    in_any_pair[nk_all] = True
    in_any_pair[t_all] = True
    pair_cells = np.unique(np.concatenate([nk, t]))
    n_rand = max(0, a.max_cells - pair_cells.size)
    ntok = emb.ntokens
    rng = np.random.default_rng(a.cell_seed)
    rand = []
    for cls in ("NK", "T"):
        want = n_rand // 2
        pool = ev[(key[cells][ev] == cls) & (ctype[cells][ev] != GDT158) & ~in_any_pair[ev]]
        ref = ntok[pair_cells[key[cells][pair_cells] == cls]]
        edges = np.quantile(ref, [0.2, 0.4, 0.6, 0.8])
        sp = np.digitize(ntok[pool], edges)
        sr = np.digitize(ref, edges)
        for q in range(5):
            share = np.mean(sr == q)
            cand = pool[sp == q]
            m = min(cand.size, int(round(want * share)))
            if m:
                rand.append(rng.choice(cand, m, replace=False))
    rand = np.sort(np.concatenate(rand)) if rand else np.zeros(0, np.int64)
    sel = np.concatenate([pair_cells, rand])
    role = np.array(["pair"] * pair_cells.size + ["random"] * rand.size)
    o = np.argsort(sel)
    return cells[sel[o]], role[o], {"nk": cells[nk], "t": cells[t]}


# ============================================================================ embed

def layer_forward_with_attn(layer, x, pad_mask):
    """nn.TransformerEncoderLayer forward (eval; dropout inactive) that also returns the head-averaged
    attention weights [B, L, L]. Same arithmetic as the module's slow path."""
    if layer.norm_first:
        xn = layer.norm1(x)
        sa, w = layer.self_attn(xn, xn, xn, key_padding_mask=pad_mask, need_weights=True, average_attn_weights=True)
        x = x + layer.dropout1(sa)
        x = x + layer._ff_block(layer.norm2(x))
    else:
        sa, w = layer.self_attn(x, x, x, key_padding_mask=pad_mask, need_weights=True, average_attn_weights=True)
        x = layer.norm1(x + layer.dropout1(sa))
        x = layer.norm2(x + layer._ff_block(x))
    return x, w


def attn_pool(h, w, mask):
    """sum_t w_t h_t / sum_t w_t, w_t = mean over real queries of the attention token t receives."""
    q = mask.to(w.dtype)
    recv = (w * q[:, :, None]).sum(1) / q.sum(1, keepdim=True).clamp(min=1.0)  # [B, L]
    recv = recv * mask.to(recv.dtype)
    return (recv[..., None] * h).sum(1) / recv.sum(1, keepdim=True).clamp(min=1e-12)


def embed(a, meta_rna, sel, gene_tok, out_path, emb):
    import torch

    from teddy_mm.device import resolve_device
    from teddy_mm.teddy_encoder import (load_gene_medians, load_pad_id, load_teddy, load_vocab, median_factors,
                                        official_values, rank_encode_official)
    torch.set_num_threads(_THREADS)
    device = resolve_device(a.device)
    seq_len = int(emb.manifest.get("seq_len", 2048))
    norm = float(emb.manifest.get("normalize_total", 1e4))
    vocab = load_vocab(a.ckpt)
    pad_id = load_pad_id(a.ckpt, vocab)
    factors = median_factors(meta_rna["rna_names"], load_gene_medians(a.medians))
    n = len(sel)
    tokens = np.full((n, seq_len), pad_id, np.int64)
    attn = np.zeros((n, seq_len), np.int64)
    for s in range(0, n, 256):
        e = min(s + 256, n)
        vals = official_values(meta_rna["rna"][sel[s:e]].toarray(), factors, norm)
        tokens[s:e], attn[s:e] = rank_encode_official(vals, meta_rna["token_ids"], max_len=seq_len, pad_id=pad_id, pad_to=seq_len)
    ntok = attn.sum(1).astype(np.int32)
    ref_ntok = emb.ntokens[np.searchsorted(emb.cells, sel)]
    mism = int(np.sum(ntok != ref_ntok))
    rel = float(np.sum(ntok.astype(np.float64) ** 2) / np.mean(emb.ntokens.astype(np.float64) ** 2))
    say(f"embed: {n} cells on {device}, autocast {a.autocast}; token-count mismatch vs official run: {mism}; "
        f"rough estimate {OFFICIAL_SEC_PER_CELL * n / 60:.1f} min at the official run's mean speed, "
        f"{OFFICIAL_SEC_PER_CELL * rel / 60:.1f} min scaled by these cells' tokens^2 (+ attention-weight overhead)")
    G = len(gene_tok)
    pos = np.full((n, G), -1, np.int32)
    for g, tid in enumerate(gene_tok):
        if tid is None:
            continue
        r, c = np.nonzero((tokens == tid) & (attn == 1))
        pos[r, g] = c
    model = load_teddy(a.ckpt, device)
    d = model.d_model
    tok_states = np.zeros((len(DEPTHS), n, G, d), np.float16)
    genemean = np.zeros((len(DEPTHS), n, d), np.float16)
    attnpool = np.zeros((len(DEPTHS), n, d), np.float16)
    amp = {"fp16": torch.float16, "bf16": torch.bfloat16}.get(a.autocast)
    order = np.argsort(ntok, kind="stable")
    check_diff = None
    t0 = time.time()
    with torch.no_grad():
        for s in range(0, n, a.batch_size):
            rows = order[s:s + a.batch_size]
            width = int(max(1, ntok[rows].max()))
            ids = torch.from_numpy(tokens[rows, :width]).to(device)
            mask = torch.from_numpy(attn[rows, :width]).to(device).bool()
            pad = ~mask
            with (torch.autocast(device_type=device.type, dtype=amp) if amp is not None else contextlib.nullcontext()):
                posv = torch.arange(width, device=device)
                h = model.embeddings(ids) + model.position_embeddings(posv)
                states = {"input": h}
                for li, layer in enumerate(model.encoder.layers, start=1):
                    need = li in ATTN_FROM.values()
                    if need:
                        h_new, w = layer_forward_with_attn(layer, h, pad)
                        if s == 0 and li == 1:  # the manual forward must equal the module's own forward (real tokens)
                            ref_h = layer(src=h, src_key_padding_mask=pad).float()
                            check_diff = float(((ref_h - h_new.float()).abs() * mask.unsqueeze(-1)).max())
                        for dn, src in ATTN_FROM.items():
                            if src == li:
                                base = states["input"] if dn == "input" else h_new
                                attnpool[DEPTHS.index(dn), rows] = attn_pool(base.float(), w.float(), mask).cpu().numpy().astype(np.float16)
                    else:
                        h_new = layer(src=h, src_key_padding_mask=pad)
                    h = h_new
                    if f"layer_{li}" in DEPTHS:
                        states[f"layer_{li}"] = h
            m = mask.unsqueeze(-1).float()
            for di, dn in enumerate(DEPTHS):
                hs = states[dn].float()
                genemean[di, rows] = ((hs * m).sum(1) / m.sum(1).clamp(min=1.0)).cpu().numpy().astype(np.float16)
                hc = hs.cpu().numpy()
                for g in range(G):
                    ok = pos[rows, g] >= 0
                    if ok.any():
                        tok_states[di, rows[ok], g] = hc[np.where(ok)[0], pos[rows[ok], g]].astype(np.float16)
            if (s // a.batch_size) % 50 == 0:
                say(f"  {min(s + a.batch_size, n)}/{n} cells, {time.time() - t0:.0f}s")
    runtime = time.time() - t0
    # check the gene-mean of layers 6 and 12 against the official stored layer means
    ix = np.searchsorted(emb.cells, sel)
    checks = {"manual_layer_forward_max_abs_diff_first_batch": check_diff, "ntokens_mismatch_vs_official": mism}
    for dn in ("layer_6", "layer_12"):
        li = int(dn.split("_")[1]) - 1
        ref = np.asarray(emb.layer_means[li][ix], dtype=np.float32)
        checks[f"{dn}_genemean_max_abs_diff_vs_official"] = float(np.abs(ref - genemean[DEPTHS.index(dn)].astype(np.float32)).max())
    np.savez(out_path, cells=sel, ntokens=ntok, gene_symbols=np.array([g for g, _ in a.genes]),
             gene_ensembl=np.array([e for _, e in a.genes]), gene_pos=pos, depths=np.array(DEPTHS),
             token_states=tok_states, genemean=genemean, attnpool=attnpool,
             runtime_sec=np.array(runtime), device=np.array(str(device)), autocast=np.array(a.autocast),
             checks=np.array(json.dumps(checks)))
    say(f"embed done in {runtime:.0f}s; checks {checks}")


# ============================================================================ fit

RIDGE_ALPHAS = lp.RIDGE_ALPHAS + (1e7, 1e8)  # the linear run's grid, widened for the 6,669-dim token sets


def dual_ridge(Xs, tr, va, ev, Y):
    """Ridge on train-standardised features in the dual (n_fit << n_features): W = X^T (X X^T + a I)^-1 Yc;
    alpha per protein on val MSE of the clipped prediction, as lp.ridge_fit_predict."""
    Xtr = Xs[tr].astype(np.float64)
    ym = Y[tr].mean(0, dtype=np.float64)
    Yc = Y[tr] - ym
    S, U = np.linalg.eigh(Xtr @ Xtr.T)
    S = np.clip(S, 0, None)
    Kva, Kev = Xs[va].astype(np.float64) @ Xtr.T, Xs[ev].astype(np.float64) @ Xtr.T
    UtY = U.T @ Yc
    best = [None] * Y.shape[1]
    for al in RIDGE_ALPHAS:
        A = U @ (UtY / (S + al)[:, None])
        mse = ((np.clip(Kva @ A + ym, 0, 1) - Y[va]) ** 2).mean(0)
        for j in range(Y.shape[1]):
            if best[j] is None or mse[j] < best[j][0]:
                best[j] = (float(mse[j]), al, A[:, j].copy())
    A = np.stack([b[2] for b in best], 1)
    return np.clip(Kev @ A + ym, 0, 1).astype(np.float32), [b[1] for b in best]

def fit(a, meta, emb, sel_expected, sel_role, pairs_global, npz_path, out):
    D = np.load(npz_path)
    sel = D["cells"]
    if not np.array_equal(sel, sel_expected):
        raise SystemExit(f"{npz_path} holds other cells than the current selection")
    role = sel_role
    names = meta["adt_names"]
    col = {n: j for j, n in enumerate(names)}
    p95 = lp.p95_replicate(meta)
    p95v = np.array([p95[n] for n in names], dtype=np.float32)
    Yall = np.clip(meta["adt"] / p95v, 0, 1)[:, [col[p] for p in PROTEINS]].astype(np.float32)
    Y = Yall[sel]
    loc = {int(c): i for i, c in enumerate(sel)}
    nk = np.array([loc[int(c)] for c in pairs_global["nk"] if int(c) in loc])
    t = np.array([loc[int(c)] for c in pairs_global["t"] if int(c) in loc])
    if nk.size != len(pairs_global["nk"]):
        raise SystemExit("not every pair cell was embedded")
    empty = np.zeros(0, np.int64)
    pairs = {"neighbour": {"all": (empty, empty), "no_gdT158": (nk, t)}, "random": {"all": [], "no_gdT158": []}}
    boot = lp.cluster_boot_indices(nk, a.n_boot, 0)
    boots = {"all": [], "no_gdT158": boot}
    rng = np.random.default_rng(a.cell_seed)
    rand = np.where(role == "random")[0]
    rng.shuffle(rand)
    nva = max(1, int(round(0.2 * rand.size)))
    rows = {"val": np.sort(rand[:nva]), "eval": np.where(role == "pair")[0]}
    tr = np.sort(rand[nva:])
    mg = np.abs(Y[nk] - Y[t])
    G = D["gene_pos"].shape[1]
    present = (D["gene_pos"] >= 0)
    res = {"feature_sets": {}, "gene_presence": {}}
    for g, sym in enumerate(D["gene_symbols"]):
        res["gene_presence"][str(sym)] = {"share_pair_NK": rnd(present[nk, g].mean()), "share_pair_T": rnd(present[t, g].mean()),
                                          "share_random": rnd(present[rand, g].mean())}
    margs = argparse.Namespace(mlp_hidden=256, mlp_dropout=0.1, mlp_lr=1e-3, mlp_wd=1e-4, mlp_batch=64, mlp_max_epochs=200, mlp_patience=10)
    for di, dn in enumerate(D["depths"]):
        dn = str(dn)
        feats = {"genemean": D["genemean"][di].astype(np.float32), "attnpool": D["attnpool"][di].astype(np.float32)}
        tok = D["token_states"][di].astype(np.float32).reshape(len(sel), -1)
        feats["tokens"] = np.concatenate([tok, present.astype(np.float32)], 1)
        feats["tokens+genemean"] = np.concatenate([feats["tokens"], feats["genemean"]], 1)
        res["feature_sets"][dn] = {}
        cache = {}
        for fs in FEATURE_SETS:
            X = feats[fs]
            Xs = lp.standardise(X, tr)
            blk = {}
            p_r, alphas = dual_ridge(Xs, tr, rows["val"], rows["eval"], Y)
            preds = {"ridge": p_r}
            blk["ridge_alpha"] = dict(zip(PROTEINS, alphas))
            blk["ridge_alpha_at_grid_edge"] = [p for p, al in zip(PROTEINS, alphas) if al in (min(RIDGE_ALPHAS), max(RIDGE_ALPHAS))]
            ps = [nl.fit_mlp(Xs, Y, tr, rows["val"], rows["eval"], sd, margs)[0] for sd in a.mlp_seeds]
            preds["mlp"] = np.mean(ps, 0).astype(np.float32)
            p_k, kinfo, _ = nl.knn_probe(X, Xs, "cos_std", Y, rows, tr) if tr.size >= max(nl.KNN_GRID) else (None, None, None)
            if p_k is not None:
                preds["knn"] = p_k
                blk["knn_k"] = kinfo["k"]
            for ro, pe in preds.items():
                full = np.full((len(sel), len(PROTEINS)), np.nan, np.float32)
                full[rows["eval"]] = pe
                g, gc = lp.gap_block(full, Y, pairs, boots)
                cache[(fs, ro)] = gc["no_gdT158"]
                blk[ro] = {"pair_cell_pearson": {p: rnd(pearson(pe[:, j], Y[rows["eval"], j])) for j, p in enumerate(PROTEINS)},
                           "gap_no_gdT158": {p: g["no_gdT158"][p] for p in PROTEINS}}
            res["feature_sets"][dn][fs] = blk
            say(f"{dn} {fs}: " + " | ".join(f"{ro} " + " ".join(f"{p}={blk[ro]['gap_no_gdT158'][p]['ratio']}" for p in PROTEINS)
                                            for ro in preds))
        res["feature_sets"][dn]["paired_vs_genemean"] = {
            f"{fs} - genemean ({ro})": nl.paired(cache[(fs, ro)], cache[("genemean", ro)], mg, boot)
            for fs in FEATURE_SETS[1:] for ro in ("ridge", "mlp", "knn") if (fs, ro) in cache}
    R = {"smoke_test": bool(a.smoke), "smoke_note": a.smoke, "created": time.strftime("%Y-%m-%d %H:%M:%S"),
         "script_version": SCRIPT_VERSION, "protocol": "probes fit on random site4 NK/T cells (80/20 fit/val), scored on the no_gdT158 pair cells",
         "n_cells": int(len(sel)), "n_pair_cells": int(rows["eval"].size), "n_fit": int(tr.size), "n_val": int(rows["val"].size),
         "n_pairs": int(nk.size), "embed_checks": json.loads(str(D["checks"])), "embed_runtime_sec": float(D["runtime_sec"]),
         "device": str(D["device"]), "genes": dict(zip(map(str, D["gene_symbols"]), map(str, D["gene_ensembl"]))), **res}
    (out / "mode_a_token_results.json").write_text(json.dumps(R, indent=1))
    lines = [f"# Mode A step 3: token-level probes ({SCRIPT_VERSION})", "", R["protocol"] + ".", "",
             "| depth | features | read-out | " + " | ".join(PROTEINS) + " |", "|---|---|---|" + "---|" * len(PROTEINS)]
    for dn, fsb in res["feature_sets"].items():
        for fs in FEATURE_SETS:
            for ro in ("ridge", "mlp", "knn"):
                if ro in fsb[fs]:
                    lines.append(f"| {dn} | {fs} | {ro} | " + " | ".join(nl.cell_ratio(fsb[fs][ro]["gap_no_gdT158"], p) for p in PROTEINS) + " |")
    (out / "README_token.md").write_text("\n".join(lines) + "\n")
    say(f"wrote {out / 'mode_a_token_results.json'}")


# ============================================================================ main

def main(argv=None) -> int:
    a = parse_args(argv)
    from threadpoolctl import threadpool_limits
    threadpool_limits(_THREADS)
    out = a.out_dir
    out.mkdir(parents=True, exist_ok=True)
    meta = lp.load_meta(a.processed, need_rna=a.stage in ("embed", "all"))
    emb = lp.Embedding(a.embed_dir, None, len(meta["split"]))
    sel, role, pairs_global = select_cells(a, meta, emb)
    if a.limit:
        keep = np.arange(min(a.limit, len(sel)))
        sel, role = sel[keep], role[keep]
    say(f"cells: {len(sel)} ({int(np.sum(role == 'pair'))} pair cells, {int(np.sum(role == 'random'))} random)")
    np.savez(out / "token_cells.npz", cells=sel, role=role, nk=pairs_global["nk"], t=pairs_global["t"])
    npz = out / "token_states.npz"
    if a.stage in ("embed", "all"):
        from teddy_mm.teddy_encoder import load_vocab
        vocab = load_vocab(a.ckpt)
        rna_names = [str(x) for x in meta["rna_names"]]
        gene_tok = []
        for sym, ens in a.genes:
            tid = vocab.get(ens)
            ok = ens in rna_names
            say(f"  gene {sym} {ens}: in processed genes {ok}, vocab id {tid}")
            gene_tok.append(tid if ok else None)
        embed(a, meta, sel, gene_tok, npz, emb)
    if a.stage in ("fit", "all"):
        if a.limit:
            say("smoke run with --limit: skipping the fit stage (the pairs are incomplete)")
            return 0
        fit(a, meta, emb, sel, role, pairs_global, npz, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
