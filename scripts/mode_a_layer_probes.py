#!/usr/bin/env python3
"""Mode A, first question: does TEDDY's representation lose the NK-T difference, or does our head?

For every depth of the frozen TEDDY-G encoder -- the input embedding (token + position embedding,
"input"), the output of each encoder layer 1..L, and the final z that the head reads -- this script
fits simple linear read-outs on training-site cells and scores them on the held-out site (site4):

  1. NK vs T probe: logistic regression on the standardised layer (L2 penalty, C picked on the val
     split by log-loss, class_weight balanced, fixed seed), site4 AUC / accuracy / balanced accuracy,
     and on the NK-T look-alike pairs the share of pairs in which the NK cell gets the higher P(NK).
  2. Protein probes: ridge regression from the standardised layer to measured CD56, CD94, CD335 and
     CD3 (target = measured ADT / train p95, clipped to [0, 1], the scale of the Experiment 6 gap
     tables; alpha per protein picked on val by MSE), site4 Pearson, and on the look-alike pairs the
     ratio  median |predicted NK - predicted T| / median |measured NK - measured T|  per protein.
  3. Geometry: the median cosine of the look-alike pairs and of random NK x T pairs at that layer
     (raw and centred on the layer's train mean), and how many site4 k-NN pairs at that layer are NK-T.

Look-alike pairs are built exactly as redesign/experiment-6/s8_pair_gap_nogdt.py and
s9_random_pair_control.py: site4 cells (split == 'test'), k = 10 cosine neighbours on the raw final
z, every neighbour edge kept once, NK-T pairs under the fixed annotation key of s2_lookalike.py
(gdT CD158b+ counts as T); the "no_gdT158" variant drops pairs that involve a 'gdT CD158b+' cell.
The pairs are fixed at z and followed back through the layers. Random pairs: 20,000 NK x T site4
pairs per seed (seeds 0-4), neighbour pairs excluded.

With --head-ckpt (04_train.py best.pt) or --head-pred (s1_predict.py pred_all.npy) the same
numbers are computed for the head's predictions (the reference the layer curves are read against),
plus an NK vs T probe on the head's 134 predicted proteins and, as a ceiling, on the 134 measured ones.

How to read the output -- see the README.md this script writes next to its results (the text is
READING_GUIDE below). In one line: if a linear probe on z (or on an earlier layer) keeps the NK-T gap
while the head's predictions compress it, the loss is in the head; if the probe on z compresses it
too, the loss is in the representation (and the layer curve shows at which depth it goes).

Inputs
  --processed   folder with cite_arrays.npz + meta.json (split, sites, cell_types, adt, adt_names,
                adt_size_factor; rna + token_ids + rna_names only for the input layer)
  --embed-dir   folder written by 03_embed_rna.py --save-layer-means: z_rna.npy,
                z_rna_layer_means.npy [L, n, d] float16, z_rna_ntokens.npy, z_rna_manifest.json
  --shard-dir   alternative to --embed-dir: read the embedding shards directly (read-only). A shard
                set that does not cover every cell is a partial input and is only accepted with --smoke.
  --ckpt, --medians   TEDDY-G checkpoint + gene medians: needed for the input layer only. Its gene-mean
                is recomputed on the CPU with the same official tokenisation (03_embed_rna.py stores
                layer outputs 1..L only); without --ckpt the input layer is left out and that is recorded.
Outputs (--out-dir)
  mode_a_probes.json       layer results (cached; a rerun with the same inputs reuses them)
  mode_a_pair_gaps.npz     per-pair gaps, for the paired bootstrap against the head
  mode_a_head_reference.json   head reference (only with --head-ckpt / --head-pred)
  mode_a_results.json      everything above in one file
  fig_*.svg                one figure per metric vs layer
  README.md                the numbers as tables + how to read them
CPU only; --threads caps BLAS / torch threads (default 3).
"""
from __future__ import annotations

import os
import sys


def _early_threads(argv) -> int:
    n = 3
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
import hashlib  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SCRIPT_VERSION = "mode_a_layer_probes v1"
PROTEINS = ["CD56", "CD94", "CD335", "CD3"]
VARIANTS = ["all", "no_gdT158"]
GDT158 = "gdT CD158b+"
LR_C_GRID = (0.01, 0.1, 1.0)
RIDGE_ALPHAS = (1e-2, 1e-1, 1.0, 10.0, 1e2, 1e3, 1e4, 1e5, 1e6)

# ---- fixed annotation key: copied verbatim from redesign/experiment-6/s2_lookalike.py (KEY_VARIANT main)
T_TYPES = ['CD4+ T activated', 'CD4+ T naive', 'CD8+ T naive', 'CD8+ T CD57+ CD45RO+',
           'CD8+ T CD57+ CD45RA+', 'CD8+ T TIGIT+ CD45RO+', 'CD4+ T activated integrinB7+',
           'CD8+ T TIGIT+ CD45RA+', 'CD8+ T CD49f+', 'CD8+ T CD69+ CD45RO+', 'MAIT',
           'CD8+ T CD69+ CD45RA+', 'gdT CD158b+', 'T reg', 'gdT TCRVD2+',
           'CD4+ T CD314+ CD45RA+', 'dnT', 'CD8+ T naive CD127+ CD26- CD101-', 'T prog cycling']
ANNOTATION_MAP = {**{t: 'T' for t in T_TYPES},
                  'NK': 'NK', 'NK CD158e1+': 'NK',
                  'Naive CD20+ B IGKC+': 'B', 'Naive CD20+ B IGKC-': 'B', 'Transitional B': 'B',
                  'B1 B IGKC+': 'B', 'B1 B IGKC-': 'B',
                  'CD14+ Mono': 'myeloid', 'CD16+ Mono': 'myeloid', 'cDC2': 'myeloid', 'cDC1': 'myeloid',
                  'Reticulocyte': 'erythroid', 'Erythroblast': 'erythroid', 'Normoblast': 'erythroid',
                  'Proerythroblast': 'erythroid',
                  'ILC': None, 'ILC1': None, 'pDC': None, 'HSC': None, 'G/M prog': None,
                  'Lymph prog': None, 'MK/E prog': None, 'Plasma cell IGKC+': None,
                  'Plasma cell IGKC-': None, 'Plasmablast IGKC+': None, 'Plasmablast IGKC-': None}


# ============================================================================ args

def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--processed", type=Path, required=True, help="dir with cite_arrays.npz + meta.json")
    p.add_argument("--embed-dir", type=Path, default=None, help="dir with z_rna.npy + z_rna_layer_means.npy")
    p.add_argument("--shard-dir", type=Path, default=None, help="read embedding shards directly instead of --embed-dir")
    p.add_argument("--ckpt", type=Path, default=None, help="TEDDY-G checkpoint dir (input layer only)")
    p.add_argument("--medians", type=Path, default=None, help="TEDDY gene medians JSON (input layer only)")
    p.add_argument("--head-ckpt", type=Path, default=None, help="04_train.py best.pt (head reference)")
    p.add_argument("--head-pred", type=Path, default=None, help="alternative: [n_cells, n_adt] head predictions (s1_predict pred_all.npy)")
    p.add_argument("--p95-json", type=Path, default=None, help="optional p95_all.json to check the p95 replicate against")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--threads", type=int, default=3)
    p.add_argument("--k", type=int, default=10)
    p.add_argument("--n-random", type=int, default=20000)
    p.add_argument("--random-seeds", default="0,1,2,3,4")
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--seq-len", type=int, default=None, help="input-layer tokenisation length (default: manifest, else 2048)")
    p.add_argument("--normalize-total", type=float, default=None, help="default: manifest, else 1e4")
    p.add_argument("--no-own-knn", action="store_true", help="skip the per-layer site4 k-NN count")
    p.add_argument("--force", action="store_true", help="recompute the layer probes even if a matching cache exists")
    p.add_argument("--smoke", default=None, metavar="NOTE",
                   help="mark every output as a SMOKE TEST (synthetic or partial inputs); NOTE says what was used")
    args = p.parse_args(argv)
    if (args.embed_dir is None) == (args.shard_dir is None):
        p.error("give exactly one of --embed-dir / --shard-dir")
    if args.head_ckpt is not None and args.head_pred is not None:
        p.error("give at most one of --head-ckpt / --head-pred")
    if args.ckpt is not None and args.medians is None:
        p.error("--ckpt (input layer) needs --medians")
    args.random_seeds = [int(s) for s in str(args.random_seeds).split(",") if s != ""]
    return args


# ============================================================================ small utils

def file_sig(path: Path | None) -> dict | None:
    if path is None:
        return None
    path = Path(path)
    if not path.exists():
        return {"path": str(path), "exists": False}
    st = path.stat()
    return {"path": str(path.resolve()), "size": st.st_size, "mtime_ns": st.st_mtime_ns}


def rnd(x, n=4):
    if x is None:
        return None
    x = float(x)
    if not np.isfinite(x):
        return None
    return round(x, n)


def pearson(a: np.ndarray, b: np.ndarray) -> float | None:
    a = np.asarray(a, np.float64)
    b = np.asarray(b, np.float64)
    if a.size < 3 or a.std() == 0 or b.std() == 0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def say(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ============================================================================ data

def load_meta(processed: Path, need_rna: bool) -> dict:
    try:
        npz = np.load(processed / "cite_arrays.npz", allow_pickle=False)
        _ = npz["adt_names"]
    except ValueError:  # object-dtype names, as teddy_mm.data.load_prepared handles
        npz = np.load(processed / "cite_arrays.npz", allow_pickle=True)
    out = {k: npz[k] for k in ("split", "sites", "cell_types", "adt", "adt_size_factor")}
    out["adt_names"] = [str(x) for x in npz["adt_names"]]
    out["split"] = out["split"].astype(str)
    out["sites"] = out["sites"].astype(str)
    out["cell_types"] = out["cell_types"].astype(str)
    out["adt"] = out["adt"].astype(np.float32)
    if need_rna:
        from scipy import sparse
        out["rna"] = sparse.csr_matrix((npz["rna_data"], npz["rna_indices"], npz["rna_indptr"]), shape=tuple(npz["rna_shape"]))
        out["token_ids"] = npz["token_ids"]
        out["rna_names"] = npz["rna_names"]
    return out


def p95_replicate(meta: dict) -> dict:
    """p95 normalisers exactly as redesign/experiment-6/s1_predict.py (and export_cite_events.py):
    rng seed 17; one draw of <= 2000 val cells outside site4 first, then 8,000 train cells."""
    split, sites, adt, names = meta["split"], meta["sites"], meta["adt"], meta["adt_names"]
    rng = np.random.default_rng(17)
    ood_idx = np.where((split == "val") & (sites != "site4"))[0]
    _ = rng.choice(ood_idx, size=min(2000, ood_idx.size), replace=False)
    train_idx = np.where(split == "train")[0]
    train_sub = rng.choice(train_idx, size=min(8000, train_idx.size), replace=False)
    p95 = {}
    for j, n in enumerate(names):
        v = float(np.percentile(adt[train_sub, j], 95))
        p95[n] = v if v >= 1e-4 else 1.0
    return p95


class Embedding:
    """Final z and the per-layer gene-means, for the cells the embedding run covered."""

    def __init__(self, embed_dir: Path | None, shard_dir: Path | None, n_total: int):
        self.manifest = {}
        self.source = {}
        if embed_dir is not None:
            lm_path = embed_dir / "z_rna_layer_means.npy"
            z_path = embed_dir / "z_rna.npy"
            if not lm_path.exists():
                raise SystemExit(f"missing {lm_path}: the embedding was run without --save-layer-means")
            self.layer_means = np.load(lm_path, mmap_mode="r")
            self.z = np.load(z_path, mmap_mode="r")
            nt = embed_dir / "z_rna_ntokens.npy"
            self.ntokens = np.load(nt) if nt.exists() else None
            cp = embed_dir / "z_rna_cells.npy"
            self.cells = np.load(cp).astype(np.int64) if cp.exists() else np.arange(self.z.shape[0], dtype=np.int64)
            mp = embed_dir / "z_rna_manifest.json"
            self.manifest = json.loads(mp.read_text()) if mp.exists() else {}
            self.source = {"kind": "embed_dir", "z_rna": file_sig(z_path), "layer_means": file_sig(lm_path)}
        else:
            files = sorted(shard_dir.glob("shard_[0-9][0-9][0-9][0-9].npz"))
            if not files:
                raise SystemExit(f"no shard_NNNN.npz in {shard_dir}")
            nums = [int(f.stem.split("_")[1]) for f in files]
            idx, z, lm, nt, cfg = [], [], [], [], set()
            for f in files:
                with np.load(f) as d:
                    if "layer_means" not in d.files:
                        raise SystemExit(f"{f} has no layer_means (embedding run without --save-layer-means)")
                    idx.append(d["idx"].astype(np.int64))
                    z.append(d["z_genemean"].astype(np.float32))
                    lm.append(d["layer_means"])
                    nt.append(d["ntokens"].astype(np.int32))
                    cfg.add(str(d["cfg_hash"]))
            if len(cfg) != 1:
                raise SystemExit(f"shards in {shard_dir} come from different configurations: {sorted(cfg)}")
            self.cells = np.concatenate(idx)
            self.z = np.concatenate(z)
            self.layer_means = np.concatenate(lm, axis=1)
            self.ntokens = np.concatenate(nt)
            self.source = {"kind": "shard_dir", "shard_dir": str(shard_dir.resolve()), "shard_numbers": nums,
                           "cfg_hash": sorted(cfg)[0],
                           "shards": [{"name": f.name, "size": f.stat().st_size} for f in files]}
        if self.layer_means.shape[1] != len(self.cells) or self.z.shape[0] != len(self.cells):
            raise SystemExit(f"shape mismatch: layer_means {self.layer_means.shape}, z {self.z.shape}, cells {len(self.cells)}")
        if len(np.unique(self.cells)) != len(self.cells) or self.cells.max(initial=-1) >= n_total:
            raise SystemExit("embedding cell index is not a unique subset of the processed cells")
        self.n_layers = int(self.layer_means.shape[0])
        self.d = int(self.layer_means.shape[2])
        self.partial = len(self.cells) < n_total

    def layer(self, i: int) -> np.ndarray:
        return np.asarray(self.layer_means[i], dtype=np.float32)

    def z_head_input(self) -> np.ndarray:
        z = np.asarray(self.z, dtype=np.float32)
        return z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)  # as 04_train.py / s1_predict.py


def input_layer_means(meta: dict, cells: np.ndarray, ckpt: Path, medians: Path, seq_len: int,
                      normalize_total: float, check_ntokens: np.ndarray | None) -> tuple[np.ndarray, dict]:
    """Gene-mean of the input to encoder layer 1 (token embedding + position embedding), recomputed on
    the CPU with the official tokenisation of 03_embed_rna.py. mean_t (E[id_t] + P[t]) over the n real
    tokens = (sum of the embeddings of the cell's top genes + sum_{t<n} P[t]) / n."""
    import torch
    from safetensors.numpy import load_file
    from scipy import sparse

    from teddy_mm.teddy_encoder import (load_gene_medians, load_pad_id, load_vocab, median_factors,
                                        official_values, rank_encode_official)

    t0 = time.time()
    w = load_file(str(ckpt / "model.safetensors"))
    E = w["embeddings.weight"].astype(np.float32)
    Pm = w["position_embeddings.weight"].astype(np.float32)
    E64 = E.astype(np.float64)
    cumP = np.cumsum(Pm.astype(np.float64), axis=0)
    vocab = load_vocab(ckpt)
    pad_id = load_pad_id(ckpt, vocab)
    factors = median_factors(meta["rna_names"], load_gene_medians(medians))
    rna, token_ids = meta["rna"], meta["token_ids"]
    V, d = E.shape
    out = np.zeros((len(cells), d), dtype=np.float32)
    ntok = np.zeros(len(cells), dtype=np.int32)
    direct_diff = None
    for s in range(0, len(cells), 256):
        e = min(s + 256, len(cells))
        vals = official_values(rna[cells[s:e]].toarray(), factors, normalize_total)
        tokens, attn = rank_encode_official(vals, token_ids, max_len=seq_len, pad_id=pad_id, pad_to=seq_len)
        n = attn.sum(1).astype(np.int64)
        r, c = np.nonzero(attn)
        M = sparse.csr_matrix((np.ones(r.size, np.float64), (r, tokens[r, c])), shape=(e - s, V))
        gsum = np.asarray(M @ E64)
        psum = np.where(n[:, None] > 0, cumP[np.maximum(n - 1, 0)], 0.0)
        out[s:e] = ((gsum + psum) / np.maximum(n, 1)[:, None]).astype(np.float32)
        ntok[s:e] = n
        if s == 0:  # direct check on the first cells: embeddings(ids) + position_embeddings(pos), masked mean
            b = min(8, e - s)
            ids = torch.from_numpy(tokens[:b])
            h = torch.from_numpy(E)[ids] + torch.from_numpy(Pm)[torch.arange(seq_len)][None]
            m = torch.from_numpy(attn[:b]).unsqueeze(-1).float()
            ref = ((h * m).sum(1) / m.sum(1).clamp(min=1.0)).numpy()
            direct_diff = float(np.abs(ref - out[:b]).max())
    info = {"runtime_sec": round(time.time() - t0, 1), "seq_len": seq_len, "normalize_total": normalize_total,
            "ckpt": str(ckpt), "medians": str(medians), "max_abs_diff_vs_direct_torch_first8": direct_diff}
    if check_ntokens is not None:
        mism = int(np.sum(check_ntokens.astype(np.int64) != ntok))
        info["ntokens_mismatch_vs_embedding_run"] = mism
        if mism:
            say(f"WARN input layer: {mism} cells get a different token count than the embedding run")
    return out, info


# ============================================================================ pairs

def knn_pairs(Z: np.ndarray, k: int) -> np.ndarray:
    """Undirected k-NN edge set (each i -> its k cosine neighbours, stored once as (lo, hi)); local indices."""
    zn = Z / (np.linalg.norm(Z, axis=1, keepdims=True) + 1e-8)
    n = zn.shape[0]
    k = min(k, n - 1)
    codes = []
    for s in range(0, n, 1024):
        e = min(s + 1024, n)
        S = zn[s:e] @ zn.T
        S[np.arange(e - s), np.arange(s, e)] = -np.inf
        part = np.argpartition(-S, kth=k - 1, axis=1)[:, :k]
        a = np.repeat(np.arange(s, e, dtype=np.int64), k)
        b = part.ravel().astype(np.int64)
        codes.append(np.minimum(a, b) * n + np.maximum(a, b))
    codes = np.unique(np.concatenate(codes))
    return np.stack([codes // n, codes % n], 1)


def build_pairs(z_eval: np.ndarray, eval_rows: np.ndarray, key: np.ndarray, ctype: np.ndarray,
                k: int, n_random: int, seeds: list[int]) -> dict:
    P = eval_rows[knn_pairs(z_eval, k)]
    ka, kb = key[P[:, 0]], key[P[:, 1]]
    m = ((ka == "NK") & (kb == "T")) | ((ka == "T") & (kb == "NK"))
    Q = P[m]
    nk = np.where(key[Q[:, 0]] == "NK", Q[:, 0], Q[:, 1])
    t = np.where(key[Q[:, 0]] == "NK", Q[:, 1], Q[:, 0])
    gd = (ctype[nk] == GDT158) | (ctype[t] == GDT158)
    nbr_codes = set((P[:, 0].astype(np.int64) * (1 << 32) + P[:, 1]).tolist())
    out = {"n_knn_pairs": int(len(P)), "neighbour": {}, "random": {}}
    nk_all = eval_rows[key[eval_rows] == "NK"]
    t_all = eval_rows[key[eval_rows] == "T"]
    for v in VARIANTS:
        sel = np.ones(len(nk), bool) if v == "all" else ~gd
        out["neighbour"][v] = (nk[sel], t[sel])
        nk_pool = nk_all if v == "all" else nk_all[ctype[nk_all] != GDT158]
        t_pool = t_all if v == "all" else t_all[ctype[t_all] != GDT158]
        rl = []
        for sd in seeds:
            if nk_pool.size == 0 or t_pool.size == 0:
                rl.append((np.zeros(0, np.int64), np.zeros(0, np.int64)))
                continue
            rng = np.random.default_rng(sd)
            o = rng.choice(nk_pool, n_random)
            tt = rng.choice(t_pool, n_random)
            lo, hi = np.minimum(o, tt), np.maximum(o, tt)
            keep = np.array([int(a) * (1 << 32) + int(b) not in nbr_codes for a, b in zip(lo, hi)], dtype=bool)
            rl.append((o[keep], tt[keep]))
        out["random"][v] = rl
    return out


def cluster_boot_indices(nk: np.ndarray, n_boot: int, seed: int) -> list[np.ndarray]:
    """Bootstrap over NK cells: each resample draws NK cells with replacement and keeps all their pairs."""
    if nk.size == 0:
        return []
    u, inv = np.unique(nk, return_inverse=True)
    groups = [np.where(inv == g)[0] for g in range(len(u))]
    rng = np.random.default_rng(seed)
    return [np.concatenate([groups[g] for g in rng.integers(0, len(u), len(u))]) for _ in range(n_boot)]


# ============================================================================ probes

def standardise(X: np.ndarray, rows: np.ndarray) -> np.ndarray:
    mu = X[rows].mean(0, dtype=np.float64)
    sd = X[rows].std(0, dtype=np.float64)
    sd[sd < 1e-8] = 1.0
    return ((X - mu) / sd).astype(np.float32)


def nkt_probe(X: np.ndarray, key: np.ndarray, ctype: np.ndarray, rows: dict, pairs: dict, seed: int) -> dict:
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import log_loss, roc_auc_score

    tr = rows["train_nkt"]
    Xs = standardise(X, tr)
    ytr = (key[tr] == "NK").astype(int)
    best, n_warn = None, 0
    va = rows["val_nkt"]
    grid = LR_C_GRID if (va.size and len(set(key[va])) == 2) else (1.0,)
    for C in grid:
        with warnings.catch_warnings(record=True) as wl:
            warnings.simplefilter("always", ConvergenceWarning)
            clf = LogisticRegression(C=C, max_iter=3000, class_weight="balanced", random_state=seed)
            clf.fit(Xs[tr], ytr)
        n_warn += sum(issubclass(w.category, ConvergenceWarning) for w in wl)
        ll = log_loss((key[va] == "NK").astype(int), clf.predict_proba(Xs[va])[:, 1], labels=[0, 1]) if len(grid) > 1 else 0.0
        if best is None or ll < best[0]:
            best = (ll, C, clf)
    clf = best[2]
    ev = rows["eval_nkt"]
    pnk = np.full(X.shape[0], np.nan)
    pnk[rows["eval"]] = clf.predict_proba(Xs[rows["eval"]])[:, 1]
    y = key[ev] == "NK"
    p = pnk[ev]
    res = {"C": best[1], "C_at_grid_edge": len(grid) > 1 and best[1] in (min(grid), max(grid)),
           "val_logloss": rnd(best[0]) if len(grid) > 1 else None, "convergence_warnings": n_warn,
           "site4_n_NK": int(y.sum()), "site4_n_T": int((~y).sum()),
           "site4_auc": rnd(roc_auc_score(y, p)) if 0 < y.sum() < y.size else None,
           "site4_accuracy": rnd(np.mean((p >= 0.5) == y)),
           "site4_balanced_accuracy": rnd(0.5 * (np.mean(p[y] >= 0.5) + np.mean(p[~y] < 0.5))) if 0 < y.sum() < y.size else None}
    ng = ctype[ev] != GDT158
    yg, pg = y[ng], p[ng]
    res["site4_auc_without_gdT158_cells"] = rnd(roc_auc_score(yg, pg)) if 0 < yg.sum() < yg.size else None
    res["pairs"] = {}
    for v in VARIANTS:
        nk, t = pairs["neighbour"][v]
        if nk.size == 0:
            res["pairs"][v] = {"n_pairs": 0}
            continue
        a, b = pnk[nk], pnk[t]
        res["pairs"][v] = {"n_pairs": int(nk.size),
                           "order_rate": rnd(np.mean(a > b) + 0.5 * np.mean(a == b)),
                           "both_correct": rnd(np.mean((a >= 0.5) & (b < 0.5)))}
    return res


def ridge_fit_predict(X: np.ndarray, rows_tr: np.ndarray, rows_va: np.ndarray, rows_ev: np.ndarray,
                      Y: np.ndarray) -> tuple[np.ndarray, list]:
    """Ridge per target (alpha on val MSE of the clipped prediction); returns clipped predictions for rows_ev."""
    mu = X[rows_tr].mean(0, dtype=np.float64)
    sd = X[rows_tr].std(0, dtype=np.float64)
    sd[sd < 1e-8] = 1.0
    Xtr = (X[rows_tr] - mu) / sd
    ym = Y[rows_tr].mean(0, dtype=np.float64)
    XtX = Xtr.T @ Xtr
    Xty = Xtr.T @ (Y[rows_tr] - ym)
    del Xtr
    Xva = (X[rows_va] - mu) / sd if rows_va.size else None
    I = np.eye(XtX.shape[0])
    best = [None] * Y.shape[1]
    for a in RIDGE_ALPHAS:
        W = np.linalg.solve(XtX + a * I, Xty)
        if Xva is None:
            mse = np.zeros(Y.shape[1])
        else:
            pv = np.clip(Xva @ W + ym, 0, 1)
            mse = ((pv - Y[rows_va]) ** 2).mean(0)
        for j in range(Y.shape[1]):
            if best[j] is None or mse[j] < best[j][0]:
                best[j] = (float(mse[j]), a, W[:, j].copy())
        if Xva is None:
            break
    W = np.stack([b[2] for b in best], 1)
    pred = np.clip(((X[rows_ev] - mu) / sd) @ W + ym, 0, 1).astype(np.float32)
    return pred, [b[1] for b in best]


def gap_block(pred_full: np.ndarray, meas_full: np.ndarray, pairs: dict, boot: dict) -> tuple[dict, dict]:
    """Gap statistics per variant and protein; returns (json block, per-pair abs pred gaps for the cache)."""
    out, cache = {}, {}
    for v in VARIANTS:
        nk, t = pairs["neighbour"][v]
        blk = {"n_pairs": int(nk.size)}
        if nk.size:
            pg = np.abs(pred_full[nk] - pred_full[t])
            mg = np.abs(meas_full[nk] - meas_full[t])
            cache[v] = pg.astype(np.float32)
            bp = np.array([np.median(pg[ix], axis=0) for ix in boot[v]]).reshape(-1, len(PROTEINS))
            bm = np.array([np.median(mg[ix], axis=0) for ix in boot[v]]).reshape(-1, len(PROTEINS))
            for j, p in enumerate(PROTEINS):
                pmed, mmed = float(np.median(pg[:, j])), float(np.median(mg[:, j]))
                ok = bm[:, j] > 0
                bs = bp[ok, j] / bm[ok, j]
                blk[p] = {"pred_median_gap": rnd(pmed), "meas_median_gap": rnd(mmed),
                          "ratio": rnd(pmed / mmed) if mmed > 0 else None,
                          "ratio_ci95": [rnd(np.percentile(bs, 2.5)), rnd(np.percentile(bs, 97.5))] if bs.size else None,
                          "pred_share_NK_higher": rnd(np.mean(pred_full[nk, j] > pred_full[t, j]))}
        rr = []
        for o, tt in pairs["random"][v]:
            if o.size == 0:
                continue
            pgr = np.median(np.abs(pred_full[o] - pred_full[tt]), axis=0)
            mgr = np.median(np.abs(meas_full[o] - meas_full[tt]), axis=0)
            rr.append(np.where(mgr > 0, pgr / np.where(mgr > 0, mgr, 1), np.nan))
        if rr:
            rr = np.array(rr)
            blk["random_pairs_ratio_median_over_seeds"] = {p: rnd(np.nanmedian(rr[:, j])) for j, p in enumerate(PROTEINS)}
            blk["random_pairs_ratio_range_over_seeds"] = {p: [rnd(np.nanmin(rr[:, j])), rnd(np.nanmax(rr[:, j]))]
                                                          for j, p in enumerate(PROTEINS)}
        out[v] = blk
    return out, cache


def protein_probe(X: np.ndarray, meas4: np.ndarray, rows: dict, pairs: dict, boot: dict) -> tuple[dict, dict]:
    pred_ev, alphas = ridge_fit_predict(X, rows["train"], rows["val"], rows["eval"], meas4)
    full = np.full((X.shape[0], len(PROTEINS)), np.nan, dtype=np.float32)
    full[rows["eval"]] = pred_ev
    res = {"alpha": dict(zip(PROTEINS, alphas)),
           "alpha_at_grid_edge": [p for p, a in zip(PROTEINS, alphas) if a in (min(RIDGE_ALPHAS), max(RIDGE_ALPHAS))],
           "site4_pearson": {p: rnd(pearson(pred_ev[:, j], meas4[rows["eval"], j])) for j, p in enumerate(PROTEINS)}}
    res["gap"], cache = gap_block(full, meas4, pairs, boot)
    return res, cache


def cosine_block(X: np.ndarray, rows: dict, pairs: dict, own_knn: bool, key: np.ndarray, k: int) -> dict:
    mu = X[rows["train"]].mean(0, dtype=np.float64).astype(np.float32)
    res = {}
    for geo, Y in (("raw", X), ("centred", X - mu)):
        Yn = Y / (np.linalg.norm(Y, axis=1, keepdims=True) + 1e-8)
        cos = lambda a, b: np.einsum("ij,ij->i", Yn[a], Yn[b])  # noqa: E731
        g = {}
        for v in VARIANTS:
            nk, t = pairs["neighbour"][v]
            nb = float(np.median(cos(nk, t))) if nk.size else None
            rs = [float(np.median(cos(o, tt))) for o, tt in pairs["random"][v] if o.size]
            rm = float(np.median(rs)) if rs else None
            g[v] = {"neighbour_median_cos": rnd(nb), "random_median_cos": rnd(rm),
                    "random_median_cos_range_over_seeds": [rnd(min(rs)), rnd(max(rs))] if rs else None,
                    "distance_ratio_neighbour_over_random": rnd((1 - nb) / (1 - rm)) if (nb is not None and rm is not None and rm < 1) else None}
        res[geo] = g
    if own_knn:
        ev = rows["eval"]
        P = ev[knn_pairs(X[ev], k)]
        ka, kb = key[P[:, 0]], key[P[:, 1]]
        graded = (ka != "") & (kb != "")
        nkt = ((ka == "NK") & (kb == "T")) | ((ka == "T") & (kb == "NK"))
        nk_cells = ev[key[ev] == "NK"]
        with_t = np.zeros(len(key), bool)
        with_t[P[nkt, 0]] = True
        with_t[P[nkt, 1]] = True
        res["own_knn"] = {"k": k, "n_pairs_graded": int(graded.sum()), "n_NK_T": int(nkt.sum()),
                          "share_NK_T_of_graded": rnd(nkt.sum() / max(graded.sum(), 1)),
                          "share_NK_cells_with_a_T_neighbour": rnd(with_t[nk_cells].mean()) if nk_cells.size else None}
    return res


# ============================================================================ head reference

def head_predictions(args, emb: Embedding, meta: dict) -> tuple[np.ndarray, dict]:
    n_adt = meta["adt"].shape[1]
    split = meta["split"]
    if args.head_pred is not None:
        pred = np.load(args.head_pred).astype(np.float32)
        if pred.shape != (len(split), n_adt):
            raise SystemExit(f"--head-pred shape {pred.shape} != ({len(split)}, {n_adt})")
        return pred[emb.cells], {"kind": "head_pred", "file": file_sig(args.head_pred)}
    import torch

    from teddy_mm.models import MLP, AdtDecoder
    torch.set_num_threads(_THREADS)
    blob = torch.load(args.head_ckpt, map_location="cpu", weights_only=False)
    w0 = blob["mlp"]["net.0.weight"]
    hidden, z_dim = int(w0.shape[0]), int(w0.shape[1])
    mlp = MLP(z_dim, z_dim, hidden=hidden)
    dec = AdtDecoder(z_dim, n_adt, hidden=hidden)
    mlp.load_state_dict(blob["mlp"])
    dec.load_state_dict(blob["dec"])
    mlp.eval()
    dec.eval()
    sf_med = float(np.median(meta["adt_size_factor"][split == "train"]))
    z = emb.z_head_input()
    pred = np.zeros((z.shape[0], n_adt), dtype=np.float32)
    with torch.no_grad():
        for s in range(0, z.shape[0], 4096):
            e = min(s + 4096, z.shape[0])
            pred[s:e] = dec(mlp(torch.from_numpy(z[s:e])), torch.full((e - s,), sf_med))[0].numpy()
    return pred, {"kind": "head_ckpt", "file": file_sig(args.head_ckpt), "hidden": hidden, "z_dim": z_dim,
                  "size_factor": "train median", "train_median_adt_size_factor": sf_med,
                  "input": "z_rna.npy, L2-normalised per cell (as 04_train.py)"}


# ============================================================================ SVG figures

SERIES_LIGHT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
SERIES_DARK = ["#3987e5", "#d95926", "#199e70", "#c98500"]


def svg_line_chart(path: Path, title: str, subtitle: str, x_labels: list, series: list, y_label: str,
                   y_range=None, ref_lines=(), head_points=None, smoke=False) -> None:
    """Static line chart (one y axis). series: [{"name", "values" (None = missing), "dash"}]; head_points:
    {series name: value} drawn as hollow markers in a separate 'head' column; ref_lines: [(y, label)]."""
    W, H = 760, 400
    ml, mr, mt, mb = 60, 150, 92, 56
    pw, ph = W - ml - mr, H - mt - mb
    cats = list(x_labels) + (["head"] if head_points else [])
    ncat = len(cats)
    gap = 0.8 if head_points else 0.0

    def xpos(i):
        span = (ncat - 1) + gap
        xi = i + (gap if (head_points and i == ncat - 1) else 0.0)
        return ml + (pw * xi / span if span > 0 else pw / 2)

    vals = [v for s in series for v in s["values"] if v is not None]
    vals += [v for v in (head_points or {}).values() if v is not None] + [r[0] for r in ref_lines]
    if y_range is None:
        lo, hi = (min(vals), max(vals)) if vals else (0.0, 1.0)
        pad = 0.08 * (hi - lo if hi > lo else 1.0)
        lo, hi = lo - pad, hi + pad
    else:
        lo, hi = y_range
    ypos = lambda v: mt + ph * (1 - (v - lo) / (hi - lo))  # noqa: E731
    step = (hi - lo) / 5
    mag = 10 ** np.floor(np.log10(step)) if step > 0 else 1
    step = min((m * mag for m in (1, 2, 2.5, 5, 10) if m * mag >= step), default=step)
    ticks = np.arange(np.ceil(lo / step) * step, hi + 1e-12, step)

    css_l = "".join(f".s{i}{{stroke:{c};fill:{c}}}" for i, c in enumerate(SERIES_LIGHT))
    css_d = "".join(f".s{i}{{stroke:{c};fill:{c}}}" for i, c in enumerate(SERIES_DARK))
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
         f'font-family="system-ui,-apple-system,Segoe UI,sans-serif" role="img" aria-label="{_esc(title)}">',
         "<style>",
         ":root{--bg:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--muted:#8a8984;--grid:#e6e5e1}",
         css_l,
         "@media (prefers-color-scheme: dark){:root{--bg:#1a1a19;--ink:#ffffff;--ink2:#c3c2b7;--muted:#8f8e86;--grid:#33332f}"
         + css_d + "}",
         ".bg{fill:var(--bg)} .t{fill:var(--ink);font-size:15px;font-weight:600} .st{fill:var(--ink2);font-size:12px}"
         " .ax{fill:var(--ink2);font-size:11px} .lg{fill:var(--ink);font-size:12px} .gr{stroke:var(--grid);stroke-width:1}"
         " .ref{stroke:var(--muted);stroke-width:1.5;stroke-dasharray:5 4;fill:none} .refl{fill:var(--ink2);font-size:11px}"
         " .ln{fill:none;stroke-width:2} .ring{stroke:var(--bg);stroke-width:2} .hollow{fill:var(--bg);stroke-width:2.5}"
         " .smk{fill:var(--ink);font-size:13px;font-weight:700}",
         "</style>",
         f'<rect class="bg" x="0" y="0" width="{W}" height="{H}"/>']
    ty = 22
    if smoke:
        o.append(f'<text class="smk" x="{ml}" y="{ty}">SMOKE TEST: synthetic / partial inputs; these numbers are not results</text>')
        ty += 19
    o.append(f'<text class="t" x="{ml}" y="{ty}">{_esc(title)}</text>')
    o.append(f'<text class="st" x="{ml}" y="{ty + 17}">{_esc(subtitle)}</text>')
    # legend
    lx = ml
    for i, s in enumerate(series):
        dash = ' stroke-dasharray="6 4"' if s.get("dash") else ""
        o.append(f'<line class="ln s{i}" x1="{lx}" y1="{mt - 16}" x2="{lx + 18}" y2="{mt - 16}"{dash}/>')
        o.append(f'<text class="lg" x="{lx + 23}" y="{mt - 12}">{_esc(s["name"])}</text>')
        lx += 30 + 7.0 * len(s["name"])
    # grid + y axis
    for tv in ticks:
        tv = 0.0 if abs(tv) < 1e-9 else float(tv)
        y = ypos(tv)
        o.append(f'<line class="gr" x1="{ml}" y1="{y:.1f}" x2="{ml + pw}" y2="{y:.1f}"/>')
        o.append(f'<text class="ax" x="{ml - 8}" y="{y + 4:.1f}" text-anchor="end">{tv:.3g}</text>')
    o.append(f'<text class="ax" transform="translate(14,{mt + ph / 2:.0f}) rotate(-90)" text-anchor="middle">{_esc(y_label)}</text>')
    for i, c in enumerate(cats):
        o.append(f'<text class="ax" x="{xpos(i):.1f}" y="{mt + ph + 18}" text-anchor="middle">{_esc(c)}</text>')
    o.append(f'<text class="ax" x="{ml + pw / 2:.0f}" y="{H - 12}" text-anchor="middle">TEDDY depth (in = input embedding, 1..L = encoder layer, z = head input)</text>')
    if head_points:
        xd = (xpos(ncat - 2) + xpos(ncat - 1)) / 2
        o.append(f'<line class="gr" x1="{xd:.1f}" y1="{mt}" x2="{xd:.1f}" y2="{mt + ph}" stroke-dasharray="3 3"/>')
    for y0, lab in ref_lines:
        y = ypos(y0)
        o.append(f'<line class="ref" x1="{ml}" y1="{y:.1f}" x2="{ml + pw}" y2="{y:.1f}"/>')
        o.append(f'<text class="refl" x="{ml + pw + 6}" y="{y + 4:.1f}">{_esc(lab)}</text>')
    # series
    ends = []
    for i, s in enumerate(series):
        pts = [(xpos(j), ypos(v), v, x_labels[j]) for j, v in enumerate(s["values"]) if v is not None]
        dash = ' stroke-dasharray="6 4"' if s.get("dash") else ""
        if len(pts) > 1:
            o.append(f'<polyline class="ln s{i}"{dash} points="' + " ".join(f"{x:.1f},{y:.1f}" for x, y, _, _ in pts) + '"/>')
        for x, y, v, lab in pts:
            o.append(f'<circle class="s{i} ring" cx="{x:.1f}" cy="{y:.1f}" r="4"><title>{_esc(s["name"])} at {lab}: {v:.3f}</title></circle>')
        if head_points and head_points.get(s["name"]) is not None:
            hv = head_points[s["name"]]
            o.append(f'<circle class="s{i} hollow" cx="{xpos(ncat - 1):.1f}" cy="{ypos(hv):.1f}" r="5"><title>{_esc(s["name"])}, head: {hv:.3f}</title></circle>')
        if pts:
            ends.append([pts[-1][1], pts[-1][1], s["name"]])
    # direct labels right of the plot (short names only, no reference-line labels there), spaced >= 14 px,
    # with a leader line when a label had to move away from its line's last point
    if not ref_lines and ends and max(len(e[2]) for e in ends) <= 10:
        ends.sort()
        for j in range(1, len(ends)):
            ends[j][0] = max(ends[j][0], ends[j - 1][0] + 14)
        for y, y_true, name in ends:
            if abs(y - y_true) > 3:
                o.append(f'<line class="gr" x1="{ml + pw + 8}" y1="{y_true:.1f}" x2="{ml + pw + 20}" y2="{y:.1f}"/>')
            o.append(f'<text class="lg" x="{ml + pw + 24}" y="{y + 4:.1f}">{_esc(name)}</text>')
    o.append("</svg>")
    path.write_text("\n".join(o))


def _esc(s) -> str:
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


# ============================================================================ README

READING_GUIDE = """\
## How to read this

**The question.** On site4, NK cells and T cells that TEDDY places next to each other (the NK-T
look-alike pairs of Experiment 6) differ clearly in measured CD56, CD94, CD335 and CD3, but our head's
predictions for the two cells of a pair are much closer than the measurements. Either TEDDY's
representation no longer carries the difference, or it does and our head (MLP + NB decoder trained on z)
throws it away. The probes here read the difference straight off each depth with a linear model and
compare that with the head.

**The numbers.**
- *Gap ratio* (per protein, per pair set): median over pairs of |predicted NK - predicted T| divided by
  the median of |measured NK - measured T|, both on the p95-scaled, clipped scale of the Experiment 6
  tables. 1 = the predictions keep the measured gap; 0 = the two cells of a pair get the same prediction;
  above 1 = the read-out exaggerates the gap (possible where the measured gap is small).
  The interval is a bootstrap over NK cells (each resample keeps all pairs of the drawn NK cells); T cells
  shared between pairs are not clustered, so the interval is on the narrow side.
- *Probe minus head*: probe gap ratio at a depth minus the head's gap ratio, on the same pairs, with a
  paired bootstrap interval (same resamples for both).
- *NK vs T probe*: site4 AUC / balanced accuracy of a logistic regression trained on training-site
  cells; *pair order rate* = share of look-alike pairs in which the NK cell gets the higher P(NK)
  (0.5 = the probe cannot tell the two cells of a pair apart).
- *Cosine*: median cosine of the look-alike pairs and of random NK x T pairs at each depth; the
  distance ratio (1 - cos neighbour) / (1 - cos random) is below 1 when the pairs are closer than chance.

**Reading rule.** Take the gap ratios on the look-alike pairs (both pair sets) and compare the linear
probe with the head:
1. *Loss in the head*: at z (or at an earlier depth) the linear probe's gap ratio is clearly above the
   head's (the probe-minus-head interval lies above 0) and the NK vs T probe orders the look-alike pairs
   well (order rate well above 0.5), while the head's predictions compress the gap. The information is in
   what the head is given; the head does not use it.
2. *Loss in the representation*: the linear probe on z compresses the gap as much as the head does, or
   more (the interval includes 0 or lies below it), and the order rate at z is near 0.5. Then read the
   curve backwards: a depth where the probe still keeps the gap and a later depth where it no longer does
   is where TEDDY (as pooled here) loses it. If no depth keeps it, the difference is not linearly readable
   from any gene-mean layer, including the input embedding, which is only the set of the cell's top genes.
3. *Mixed*: proteins or pair sets can disagree (for example CD3 kept, CD56 lost). Report per protein;
   do not average across them.

**What this does not show.**
- A linear probe is a lower bound on what a depth contains. "Probe compresses too" means the difference
  is not *linearly* available in the gene-mean of that depth; a non-linear read-out or the per-token
  states (not saved) could still hold it. "Probe keeps it" is the stronger statement.
- Ridge shrinks predictions toward the mean, so any probe with imperfect fit compresses gaps. That is why
  the probe is compared with the head on the same pairs, and why the random NK x T pairs are given: a
  look-alike ratio far below the random-pair ratio for the same read-out means the compression is
  specific to look-alikes, not a general loss of the NK-T axis.
- The probes and the head are trained differently (probes: squared error on the p95-scaled, clipped
  scale, one protein at a time; head: NB likelihood on counts for all proteins with a constant size
  factor); both are scored on the same scale and the same pairs.
- The pairs are chosen as neighbours at z, so their cosine at z is high by construction; at earlier
  depths the pairs are not selected and their cosine is only meaningful next to the random-pair cosine.
- Layer L and z are the same states (z is float32 and L2-normalised as the head sees it; layer means
  are stored in float16); the two columns should agree closely, which is a check, not a finding.
- The embedding run used the settings recorded under `inputs.embedding_manifest` (e.g. fp16 autocast);
  the input embedding is recomputed here in float32 from the same tokens.
- One held-out site, one head seed, one pair construction (k = 10). The C and alpha of the probes are
  picked on the val split (a training site), never on site4.
"""


def fmt(x, n=3) -> str:
    return "" if x is None else f"{x:.{n}f}"


def write_readme(out_dir: Path, R: dict) -> None:
    L = R["layers"]
    head = R.get("head_reference")
    lines = []
    if R["smoke_test"]:
        lines += ["> **SMOKE TEST. These numbers are not results.** Inputs: " + _md((R["smoke_note"] or "").rstrip(".")) +
                  ". The run only checks that the script works end to end.", ""]
    lines += ["# Mode A: layer-wise NK-T probes", "",
              f"Written by `scripts/mode_a_layer_probes.py` ({R['config']['script_version']}), "
              f"{R['created']}. Status: **{R['status']}**.", ""]
    c = R["cells"]
    lines += [f"Cells: {c['n_present']} of {c['n_total']} with an embedding"
              + (" (partial input)" if c["partial"] else "") +
              f"; site4 NK {c['eval_NK']}, T {c['eval_T']}; train NK {c['train_NK']}, T {c['train_T']}.",
              "Look-alike pairs (k = %d): all %s, no_gdT158 %s." % (R["config"]["k"], R["pairs"]["all"]["n_pairs"],
                                                                     R["pairs"]["no_gdT158"]["n_pairs"]), ""]
    if not R["config"]["input_layer"]:
        lines += ["The input embedding is not included (no --ckpt given).", ""]
    lines += ["## NK vs T linear probe (site4)", "",
              "| depth | AUC | balanced acc. | AUC without gdT CD158b+ | pair order rate (all) | pair order rate (no_gdT158) |",
              "|---|---|---|---|---|---|"]
    for l in L:
        q = l["nk_t_probe"]
        lines.append(f"| {l['label']} | {fmt(q['site4_auc'])} | {fmt(q['site4_balanced_accuracy'])} | "
                     f"{fmt(q['site4_auc_without_gdT158_cells'])} | {fmt(q['pairs']['all'].get('order_rate'))} | "
                     f"{fmt(q['pairs']['no_gdT158'].get('order_rate'))} |")
    if head:
        for nm, q in ((f"head outputs ({head['n_proteins']} predicted proteins)", head["nk_t_probe_on_head_outputs"]),
                      ("measured proteins (ceiling)", head["nk_t_probe_on_measured_proteins"])):
            lines.append(f"| {nm} | {fmt(q['site4_auc'])} | {fmt(q['site4_balanced_accuracy'])} | "
                         f"{fmt(q['site4_auc_without_gdT158_cells'])} | {fmt(q['pairs']['all'].get('order_rate'))} | "
                         f"{fmt(q['pairs']['no_gdT158'].get('order_rate'))} |")
    lines += ["", "## Protein probes: site4 Pearson (p95-scaled, clipped)", "",
              "| depth | " + " | ".join(PROTEINS) + " |", "|---|" + "---|" * len(PROTEINS)]
    for l in L:
        lines.append(f"| {l['label']} | " + " | ".join(fmt(l["protein_probe"]["site4_pearson"][p]) for p in PROTEINS) + " |")
    if head:
        lines.append("| head | " + " | ".join(fmt(head["site4_pearson"][p]) for p in PROTEINS) + " |")
    for v in VARIANTS:
        lines += ["", f"## Gap ratio on the look-alike pairs: {v} ({R['pairs'][v]['n_pairs']} pairs)", "",
                  "Median predicted gap / median measured gap; [95% bootstrap interval]; random NK x T pairs in the last columns.", "",
                  "| depth | " + " | ".join(PROTEINS) + " | " + " | ".join(f"{p} random" for p in PROTEINS) + " |",
                  "|---|" + "---|" * (2 * len(PROTEINS))]
        rows = [(l["label"], l["protein_probe"]["gap"][v]) for l in L] + ([("head", head["gap"][v])] if head else [])
        for lab, g in rows:
            cells_ = []
            for p in PROTEINS:
                e = g.get(p) or {}
                ci = e.get("ratio_ci95")
                cells_.append(fmt(e.get("ratio")) + (f" [{fmt(ci[0])}, {fmt(ci[1])}]" if ci else ""))
            rr = g.get("random_pairs_ratio_median_over_seeds") or {}
            lines.append(f"| {lab} | " + " | ".join(cells_) + " | " + " | ".join(fmt(rr.get(p)) for p in PROTEINS) + " |")
        g0 = L[0]["protein_probe"]["gap"][v] if L else {}
        lines += ["", "Measured median gap: " + ", ".join(f"{p} {fmt((g0.get(p) or {}).get('meas_median_gap'))}" for p in PROTEINS) + "."]
    if R.get("comparisons"):
        lines += ["", "## Probe minus head (gap ratio, paired bootstrap)", "",
                  "| pair set | depth | " + " | ".join(PROTEINS) + " |", "|---|---|" + "---|" * len(PROTEINS)]
        for v in VARIANTS:
            for l in L:
                cmp_ = R["comparisons"][v].get(l["name"], {})
                cells_ = []
                for p in PROTEINS:
                    e = cmp_.get(p)
                    cells_.append("" if not e or e["diff"] is None else
                                  f"{e['diff']:+.3f}" + (f" [{e['ci95'][0]:+.3f}, {e['ci95'][1]:+.3f}]" if e.get("ci95") else ""))
                lines.append(f"| {v} | {l['label']} | " + " | ".join(cells_) + " |")
    lines += ["", "## Pair cosine", "",
              "| depth | raw: neighbour (all) | raw: random | raw: distance ratio | centred: neighbour (all) | centred: random | centred: distance ratio | own k-NN: NK-T share |",
              "|---|---|---|---|---|---|---|---|"]
    for l in L:
        cr, cc = l["cosine"]["raw"]["all"], l["cosine"]["centred"]["all"]
        ok = l["cosine"].get("own_knn") or {}
        lines.append(f"| {l['label']} | {fmt(cr['neighbour_median_cos'])} | {fmt(cr['random_median_cos'])} | "
                     f"{fmt(cr['distance_ratio_neighbour_over_random'])} | {fmt(cc['neighbour_median_cos'])} | "
                     f"{fmt(cc['random_median_cos'])} | {fmt(cc['distance_ratio_neighbour_over_random'])} | "
                     f"{fmt(ok.get('share_NK_T_of_graded'), 4)} |")
    edge = [(l["label"], l["protein_probe"]["alpha_at_grid_edge"]) for l in L if l["protein_probe"]["alpha_at_grid_edge"]]
    if edge:
        lines += ["", "Ridge alpha at the edge of its grid (the probe may be under- or over-regularised there): "
                  + "; ".join(f"{a}: {', '.join(b)}" for a, b in edge) + "."]
    cedge = [l["label"] for l in L if l["nk_t_probe"].get("C_at_grid_edge")]
    if cedge:
        lines += ["", f"Logistic-regression C at the edge of its grid {list(LR_C_GRID)} at depth: " + ", ".join(cedge) + "."]
    lines += ["", "Figures: " + ", ".join(f"`{f}`" for f in R["figures"]) + ".", "", READING_GUIDE]
    (out_dir / "README.md").write_text("\n".join(lines) + "\n")


def _md(s: str) -> str:
    return s.replace("|", "/")


# ============================================================================ main

def config_hash(args, emb_source: dict, extra: dict) -> str:
    script_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    cfg = {"version": SCRIPT_VERSION, "script_sha256": script_sha, "processed": file_sig(args.processed / "cite_arrays.npz"), "embedding": emb_source,
           "ckpt": file_sig(args.ckpt / "model.safetensors") if args.ckpt else None, "medians": file_sig(args.medians),
           "k": args.k, "n_random": args.n_random, "seeds": args.random_seeds, "n_boot": args.n_boot, "seed": args.seed,
           "own_knn": not args.no_own_knn, "smoke": args.smoke, **extra}
    return hashlib.sha256(json.dumps(cfg, sort_keys=True, default=str).encode()).hexdigest()[:16]


def main(argv=None) -> int:
    args = parse_args(argv)
    from threadpoolctl import threadpool_limits
    threadpool_limits(_THREADS)
    t_start = time.time()
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    say(f"mode_a_layer_probes: threads={_THREADS} out={out}")

    meta = load_meta(args.processed, need_rna=args.ckpt is not None)
    n_total = len(meta["split"])
    missing = sorted(set(meta["cell_types"]) - set(ANNOTATION_MAP))
    if missing:
        raise SystemExit(f"cell types not in the fixed key: {missing[:10]}")
    key = np.array([ANNOTATION_MAP[c] or "" for c in meta["cell_types"]])
    ctype = meta["cell_types"]
    split = meta["split"]
    names = meta["adt_names"]
    col = {n: j for j, n in enumerate(names)}
    lack = [p for p in PROTEINS if p not in col]
    if lack:
        raise SystemExit(f"proteins missing from adt_names: {lack}")
    p95 = p95_replicate(meta)
    p95v = np.array([p95[n] for n in names], dtype=np.float32)
    MEAS = np.clip(meta["adt"] / p95v, 0, 1)
    meas4 = MEAS[:, [col[p] for p in PROTEINS]]
    checks = {}
    if args.p95_json is not None and args.p95_json.exists():
        ref = json.loads(args.p95_json.read_text())["p95"]
        checks["p95_max_abs_diff_vs_p95_json"] = float(max(abs(ref[n] - p95[n]) for n in names if n in ref))

    emb = Embedding(args.embed_dir, args.shard_dir, n_total)
    if emb.partial and not args.smoke:
        raise SystemExit(f"the embedding covers {len(emb.cells)} of {n_total} cells; a partial run is only allowed with --smoke")
    present = np.zeros(n_total, bool)
    present[emb.cells] = True
    # work in "present-cell" coordinates: every array below is indexed by position in emb.cells
    cells = emb.cells
    keyp, ctp, spp = key[cells], ctype[cells], split[cells]
    measp, meas4p = MEAS[cells], meas4[cells]
    rows = {"train": np.where(spp == "train")[0], "val": np.where(spp == "val")[0], "eval": np.where(spp == "test")[0]}
    nkt = (keyp == "NK") | (keyp == "T")
    rows["train_nkt"] = rows["train"][nkt[rows["train"]]]
    rows["val_nkt"] = rows["val"][nkt[rows["val"]]]
    rows["eval_nkt"] = rows["eval"][nkt[rows["eval"]]]
    if rows["eval"].size == 0 or rows["train_nkt"].size == 0:
        raise SystemExit(f"no site4 (split == 'test') cells or no training NK/T cells among the {len(cells)} embedded cells")

    z_raw = np.asarray(emb.z, dtype=np.float32)
    t0 = time.time()
    pairs = build_pairs(z_raw[rows["eval"]], rows["eval"], keyp, ctp, args.k, args.n_random, args.random_seeds)
    say(f"pairs: knn {pairs['n_knn_pairs']}, NK-T all {pairs['neighbour']['all'][0].size}, "
        f"no_gdT158 {pairs['neighbour']['no_gdT158'][0].size} ({time.time() - t0:.0f}s)")
    boot = {v: cluster_boot_indices(pairs["neighbour"][v][0], args.n_boot, args.seed) for v in VARIANTS}

    manifest = emb.manifest
    seq_len = args.seq_len or int(manifest.get("seq_len", 2048))
    norm_total = args.normalize_total or float(manifest.get("normalize_total", 1e4))
    if args.ckpt is not None and manifest and manifest.get("preprocessing", "official") != "official":
        raise SystemExit("the input-layer recomputation follows the official tokenisation; this embedding is not official")
    chash = config_hash(args, emb.source, {"seq_len": seq_len, "normalize_total": norm_total})

    # ---------------------------------------------------------------- layer probes (cached)
    probes_path, gaps_path = out / "mode_a_probes.json", out / "mode_a_pair_gaps.npz"
    cached = None
    if probes_path.exists() and gaps_path.exists() and not args.force:
        c = json.loads(probes_path.read_text())
        if c.get("config_hash") == chash:
            with np.load(gaps_path) as g:
                same = all(np.array_equal(g[f"nk_{v}"], cells[pairs["neighbour"][v][0]]) and
                           np.array_equal(g[f"t_{v}"], cells[pairs["neighbour"][v][1]]) for v in VARIANTS)
            if same:
                cached = c
                say(f"reusing layer probes from {probes_path} (config {chash})")
            else:
                say("cached probes were computed on other pairs; recomputing")
    if cached is None:
        layer_specs = []
        if args.ckpt is not None:
            layer_specs.append(("input", "in"))
        layer_specs += [(f"layer_{i + 1}", str(i + 1)) for i in range(emb.n_layers)]
        layer_specs.append(("z", "z"))
        layers, gap_cache = [], {v: [] for v in VARIANTS}
        input_info = None
        for name, label in layer_specs:
            tl = time.time()
            timing = {}
            if name == "input":
                X, input_info = input_layer_means(meta, cells, args.ckpt, args.medians, seq_len, norm_total, emb.ntokens)
                timing["load"] = input_info["runtime_sec"]
            elif name == "z":
                X = emb.z_head_input()
                timing["load"] = round(time.time() - tl, 1)
            else:
                X = emb.layer(int(name.split("_")[1]) - 1)
                timing["load"] = round(time.time() - tl, 1)
            t1 = time.time()
            q = nkt_probe(X, keyp, ctp, rows, pairs, args.seed)
            timing["nk_t_probe"] = round(time.time() - t1, 1)
            t1 = time.time()
            pp, gc = protein_probe(X, meas4p, rows, pairs, boot)
            timing["protein_probe"] = round(time.time() - t1, 1)
            t1 = time.time()
            cb = cosine_block(X, rows, pairs, not args.no_own_knn, keyp, args.k)
            timing["cosine_and_knn"] = round(time.time() - t1, 1)
            for v in VARIANTS:
                gap_cache[v].append(gc.get(v, np.zeros((0, len(PROTEINS)), np.float32)))
            layers.append({"name": name, "label": label, "nk_t_probe": q, "protein_probe": pp, "cosine": cb,
                           "timing_sec": timing})
            say(f"{name:9s} AUC {q['site4_auc']} order(all) {q['pairs']['all'].get('order_rate')} "
                f"gap(all) " + " ".join(f"{p}={pp['gap']['all'].get(p, {}).get('ratio')}" for p in PROTEINS) +
                f" ({time.time() - tl:.0f}s)")
            del X
        z16 = emb.layer(emb.n_layers - 1)
        checks["last_layer_mean_vs_z_max_abs_diff"] = float(np.abs(z16 - z_raw).max())
        checks["last_layer_mean_vs_z_max_rel_diff"] = float(np.abs(z16 - z_raw).max() / (np.abs(z_raw).max() + 1e-12))
        del z16
        cached = {"config_hash": chash, "layers": layers, "input_layer": input_info, "checks": checks,
                  "probe_runtime_sec": round(time.time() - t_start, 1)}
        np.savez(gaps_path, **{f"pred_gap_{v}": np.stack(gap_cache[v]) for v in VARIANTS},
                 **{f"nk_{v}": cells[pairs["neighbour"][v][0]] for v in VARIANTS},
                 **{f"t_{v}": cells[pairs["neighbour"][v][1]] for v in VARIANTS},
                 layer_names=np.array([l["name"] for l in layers]))
        probes_path.write_text(json.dumps(cached, indent=1))
    layers = cached["layers"]
    checks = {**cached.get("checks", {}), **checks}

    # ---------------------------------------------------------------- head reference
    head = None
    comparisons = None
    if args.head_ckpt is not None or args.head_pred is not None:
        src = args.head_ckpt or args.head_pred
        if not Path(src).exists():
            raise SystemExit(f"head reference {src} does not exist")
        pred, hinfo = head_predictions(args, emb, meta)
        PRED = np.clip(pred / p95v, 0, 1)
        pred4 = PRED[:, [col[p] for p in PROTEINS]]
        full = np.full_like(pred4, np.nan)
        full[rows["eval"]] = pred4[rows["eval"]]
        gblk, hcache = gap_block(full, meas4p, pairs, boot)
        head = {"source": hinfo, "n_proteins": int(PRED.shape[1]),
                "site4_pearson": {p: rnd(pearson(pred4[rows["eval"], j], meas4p[rows["eval"], j])) for j, p in enumerate(PROTEINS)},
                "gap": gblk,
                "nk_t_probe_on_head_outputs": nkt_probe(PRED, keyp, ctp, rows, pairs, args.seed),
                "nk_t_probe_on_measured_proteins": nkt_probe(measp, keyp, ctp, rows, pairs, args.seed)}
        g = np.load(gaps_path)
        meas_gap = {v: np.abs(meas4p[pairs["neighbour"][v][0]] - meas4p[pairs["neighbour"][v][1]]) for v in VARIANTS}
        comparisons = {}
        for v in VARIANTS:
            comparisons[v] = {}
            if v not in hcache:
                continue
            pg_layers = g[f"pred_gap_{v}"]
            hg, mg = hcache[v], meas_gap[v]
            bm = np.array([np.median(mg[ix], axis=0) for ix in boot[v]]).reshape(-1, len(PROTEINS))
            bh = np.array([np.median(hg[ix], axis=0) for ix in boot[v]]).reshape(-1, len(PROTEINS))
            m0, h0 = np.median(mg, axis=0), np.median(hg, axis=0)
            for li, l in enumerate(layers):
                pl = pg_layers[li]
                bl = np.array([np.median(pl[ix], axis=0) for ix in boot[v]]).reshape(-1, len(PROTEINS))
                p0 = np.median(pl, axis=0)
                comparisons[v][l["name"]] = {}
                for j, p in enumerate(PROTEINS):
                    if m0[j] <= 0:
                        comparisons[v][l["name"]][p] = {"diff": None}
                        continue
                    ok = bm[:, j] > 0
                    bs = (bl[ok, j] - bh[ok, j]) / bm[ok, j]
                    comparisons[v][l["name"]][p] = {"diff": rnd((p0[j] - h0[j]) / m0[j]),
                                                    "ci95": [rnd(np.percentile(bs, 2.5)), rnd(np.percentile(bs, 97.5))] if bs.size else None}
        (out / "mode_a_head_reference.json").write_text(json.dumps({"config_hash": chash, "head_reference": head,
                                                                   "comparisons": comparisons}, indent=1))

    # ---------------------------------------------------------------- results, figures, README
    cnt = lambda r, c: int(np.sum(keyp[r] == c))  # noqa: E731
    R = {"smoke_test": bool(args.smoke), "smoke_note": args.smoke,
         "status": "complete (layer probes + head reference)" if head else "layer probes only (no head reference yet)",
         "created": time.strftime("%Y-%m-%d %H:%M:%S"),
         "config": {"script_version": SCRIPT_VERSION, "config_hash": chash, "k": args.k, "proteins": PROTEINS,
                    "variants": VARIANTS, "n_random_per_seed": args.n_random, "random_seeds": args.random_seeds,
                    "n_boot": args.n_boot, "seed": args.seed, "threads": _THREADS, "lr_C_grid": LR_C_GRID,
                    "ridge_alphas": RIDGE_ALPHAS, "input_layer": args.ckpt is not None, "own_knn": not args.no_own_knn,
                    "target_scale": "measured ADT / train p95 (s1_predict replicate), clipped to [0, 1]",
                    "probe_features": "each depth standardised on train cells; z = z_rna.npy L2-normalised (head input)",
                    "pairs": "site4 k-NN on raw z, NK-T under the s2_lookalike key; no_gdT158 drops pairs with a gdT CD158b+ cell",
                    "annotation_map": ANNOTATION_MAP},
         "inputs": {"processed": str(args.processed), "embedding": emb.source,
                    "embedding_manifest": {k_: manifest.get(k_) for k_ in ("preprocessing", "seq_len", "pooling", "autocast",
                                                                            "batch_size", "length_buckets", "config_hash",
                                                                            "git_commit")} if manifest else None,
                    "ckpt": str(args.ckpt) if args.ckpt else None, "medians": str(args.medians) if args.medians else None,
                    "head": str(args.head_ckpt or args.head_pred) if (args.head_ckpt or args.head_pred) else None},
         "cells": {"n_total": n_total, "n_present": int(len(cells)), "partial": bool(emb.partial),
                   "n_train": int(rows["train"].size), "n_val": int(rows["val"].size), "n_eval_site4": int(rows["eval"].size),
                   "train_NK": cnt(rows["train"], "NK"), "train_T": cnt(rows["train"], "T"),
                   "eval_NK": cnt(rows["eval"], "NK"), "eval_T": cnt(rows["eval"], "T")},
         "pairs": {v: {"n_pairs": int(pairs["neighbour"][v][0].size),
                       "n_NK_cells": int(np.unique(pairs["neighbour"][v][0]).size),
                       "n_T_cells": int(np.unique(pairs["neighbour"][v][1]).size),
                       "n_random_pairs_per_seed": [int(o.size) for o, _ in pairs["random"][v]]} for v in VARIANTS},
         "input_layer": cached.get("input_layer"), "checks": checks, "layers": layers,
         "head_reference": head, "comparisons": comparisons,
         "runtime_sec": {"layer_probes": cached.get("probe_runtime_sec"), "this_call": round(time.time() - t_start, 1)}}
    R["figures"] = write_figures(out, R)
    (out / "mode_a_results.json").write_text(json.dumps(R, indent=1))
    write_readme(out, R)
    say(f"wrote {out / 'mode_a_results.json'} ({R['status']}) in {time.time() - t_start:.0f}s")
    return 0


def write_figures(out: Path, R: dict) -> list:
    L = R["layers"]
    xl = [l["label"] for l in L]
    head = R.get("head_reference")
    smoke = R["smoke_test"]
    figs = []
    sub_pairs = f"site4 look-alike NK-T pairs (k = {R['config']['k']})"

    q = [l["nk_t_probe"] for l in L]
    series = [{"name": "site4 AUC", "values": [x["site4_auc"] for x in q]},
              {"name": "site4 balanced acc.", "values": [x["site4_balanced_accuracy"] for x in q]},
              {"name": "pair order rate (all)", "values": [x["pairs"]["all"].get("order_rate") for x in q]},
              {"name": "pair order rate (no gdT158)", "values": [x["pairs"]["no_gdT158"].get("order_rate") for x in q], "dash": True}]
    hp = None
    refs = [(0.5, "chance")]
    if head:
        h = head["nk_t_probe_on_head_outputs"]
        hp = {"site4 AUC": h["site4_auc"], "site4 balanced acc.": h["site4_balanced_accuracy"],
              "pair order rate (all)": h["pairs"]["all"].get("order_rate"),
              "pair order rate (no gdT158)": h["pairs"]["no_gdT158"].get("order_rate")}
    svg_line_chart(out / "fig_nkt_probe.svg", "NK vs T linear probe at each depth",
                   "logistic regression trained on training sites; head column = same probe on the head's "
                   + (f"{head['n_proteins']} " if head else "") + "predicted proteins",
                   xl, series, "score", y_range=(0.4, 1.0), ref_lines=refs, head_points=hp, smoke=smoke)
    figs.append("fig_nkt_probe.svg")

    series = [{"name": p, "values": [l["protein_probe"]["site4_pearson"][p] for l in L]} for p in PROTEINS]
    hp = {p: head["site4_pearson"][p] for p in PROTEINS} if head else None
    svg_line_chart(out / "fig_protein_pearson.svg", "Protein probes: site4 Pearson at each depth",
                   "ridge from the layer to measured protein (p95-scaled, clipped); hollow = our head", xl, series,
                   "Pearson r", head_points=hp, smoke=smoke)
    figs.append("fig_protein_pearson.svg")

    for v in VARIANTS:
        series = [{"name": p, "values": [(l["protein_probe"]["gap"][v].get(p) or {}).get("ratio") for l in L]} for p in PROTEINS]
        hp = {p: (head["gap"][v].get(p) or {}).get("ratio") for p in PROTEINS} if head else None
        nm = "all pairs" if v == "all" else "without gdT CD158b+ pairs"
        svg_line_chart(out / f"fig_gap_ratio_{v}.svg", f"NK-T gap kept by a linear read-out ({nm})",
                       f"median |pred NK - pred T| / median |measured gap|, {sub_pairs}; hollow = our head",
                       xl, series, "gap ratio", ref_lines=[(1.0, "measured gap"), (0.0, "no gap")], head_points=hp, smoke=smoke)
        figs.append(f"fig_gap_ratio_{v}.svg")

    for geo in ("raw", "centred"):
        c = [l["cosine"][geo] for l in L]
        series = [{"name": "look-alike pairs (all)", "values": [x["all"]["neighbour_median_cos"] for x in c]},
                  {"name": "look-alike pairs (no gdT158)", "values": [x["no_gdT158"]["neighbour_median_cos"] for x in c]},
                  {"name": "random NK x T pairs", "values": [x["all"]["random_median_cos"] for x in c], "dash": True}]
        svg_line_chart(out / f"fig_pair_cosine_{geo}.svg", f"Median cosine of NK-T pairs at each depth ({geo})",
                       "pairs are neighbours at z (so high there by construction); random pairs are the baseline"
                       + ("; centred on the layer's train mean" if geo == "centred" else ""),
                       xl, series, "median cosine", smoke=smoke)
        figs.append(f"fig_pair_cosine_{geo}.svg")

    if all(l["cosine"].get("own_knn") for l in L):
        series = [{"name": "NK-T share of graded pairs", "values": [l["cosine"]["own_knn"]["share_NK_T_of_graded"] for l in L]}]
        svg_line_chart(out / "fig_own_knn_nkt_share.svg", "How often NK and T cells are k-NN neighbours at each depth",
                       f"site4, k = {R['config']['k']} cosine neighbours computed at that depth", xl, series,
                       "share of graded neighbour pairs", smoke=smoke)
        figs.append("fig_own_knn_nkt_share.svg")
    return figs


if __name__ == "__main__":
    sys.exit(main())
