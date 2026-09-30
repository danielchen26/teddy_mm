#!/usr/bin/env python3
"""Mode A, step 2: can a NON-linear read-out get the NK-T difference out of TEDDY's pooled layers?

scripts/mode_a_layer_probes.py showed that a linear (ridge) read-out of the gene-mean at every depth
keeps only part of the measured NK-T gap on the site4 look-alike pairs. A linear probe is a lower bound
on what a depth contains, so this script repeats the protein read-out with two non-linear probes at a
few depths (default: input embedding, layers 3, 6, 9, 12 and z):

  (a) MLP: 2 hidden layers (ReLU, dropout), one network for the 4 proteins, AdamW on squared error,
      early stopping on the val split (a training site; never site4), 3 fixed seeds; the seed-mean of
      the 3 networks' predictions is reported next to each seed.
  (b) k-nearest-neighbour regression: mean target of the k nearest TRAINING cells (cosine), k picked per
      protein on val MSE. Two geometries: cos_std (train-standardised features, what the ridge and the
      MLP see) and cos_raw (raw features; at z this is exactly the geometry the pairs were selected in).

Everything else is reused from mode_a_layer_probes.py unchanged (imported, not copied): the fixed NK/T
annotation key, the p95 target scale (measured ADT / train p95, clipped to [0, 1]), the splits
(train / val / site4 = 'test'), the look-alike pairs (site4 k = 10 cosine neighbours on raw z; 'all' and
'no_gdT158'), the random NK x T control pairs, the NK-cell cluster bootstrap (same seed, so the same
resamples), the gap-ratio statistic, the input-layer recomputation and the head (04_train.py best.pt on
L2-normalised z). The script checks that the rebuilt pairs are identical to the ones stored by the linear
run and that the head's gap ratios are reproduced, and it takes the linear-probe numbers and the
per-pair linear gaps from that run (mode_a_results.json, mode_a_pair_gaps.npz) for the paired
comparisons  nonlinear - linear  and  nonlinear - head.

A built-in bias of kNN: the look-alike pairs are neighbours at z, so at z (cos_raw exactly, cos_std
nearly) the two cells of a pair share many of their nearest training cells and a kNN read-out gives them
similar predictions BY CONSTRUCTION. The script measures that directly (share of shared training
neighbours per pair, next to random NK x T pairs) and reports kNN at earlier depths, where the pairs were
not selected.

Outputs (--out-dir): mode_a_nonlinear_results.json (written after every depth, so a partial run leaves
partial results marked as such) and README_nonlinear.md (tables; an optional --reading markdown file is
appended as the reading section). CPU only (the GPU is never used); --threads caps BLAS / torch threads.
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
import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("mode_a_layer_probes", HERE / "mode_a_layer_probes.py")
lp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(lp)

SCRIPT_VERSION = "mode_a_nonlinear_probes v1"
PROTEINS, VARIANTS, GDT158 = lp.PROTEINS, lp.VARIANTS, lp.GDT158
rnd, say, pearson = lp.rnd, lp.say, lp.pearson
KNN_GRID = (5, 10, 20, 50, 100, 200)
KNN_METRICS = ("cos_std", "cos_raw")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--processed", type=Path, required=True)
    p.add_argument("--embed-dir", type=Path, required=True)
    p.add_argument("--ckpt", type=Path, default=None, help="TEDDY-G checkpoint (input layer only)")
    p.add_argument("--medians", type=Path, default=None)
    p.add_argument("--head-ckpt", type=Path, required=True, help="04_train.py best.pt used by the linear run")
    p.add_argument("--linear-dir", type=Path, required=True, help="out-dir of mode_a_layer_probes.py (mode_a_results.json + mode_a_pair_gaps.npz)")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--cache-dir", type=Path, default=None, help="where the recomputed input-layer gene-mean is cached (float32)")
    p.add_argument("--layers", default="input,3,6,9,12,z", help="depths: input, 1..L, z")
    p.add_argument("--knn-layers", default=None, help="depths for kNN (default: --layers)")
    p.add_argument("--mlp-layers", default=None, help="depths for the MLP (default: --layers)")
    p.add_argument("--knn-metrics", default=",".join(KNN_METRICS))
    p.add_argument("--mlp-seeds", default="0,1,2")
    p.add_argument("--mlp-hidden", type=int, default=256)
    p.add_argument("--mlp-dropout", type=float, default=0.1)
    p.add_argument("--mlp-lr", type=float, default=1e-3)
    p.add_argument("--mlp-wd", type=float, default=1e-4)
    p.add_argument("--mlp-batch", type=int, default=256)
    p.add_argument("--mlp-max-epochs", type=int, default=60)
    p.add_argument("--mlp-patience", type=int, default=6)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--n-boot", type=int, default=None, help="default: the linear run's n_boot")
    p.add_argument("--reading", type=Path, default=None, help="markdown appended to README_nonlinear.md as the reading")
    p.add_argument("--max-train", type=int, default=0, help="smoke only: subsample training cells")
    p.add_argument("--smoke", default=None, metavar="NOTE")
    p.add_argument("--readme-only", action="store_true", help="rebuild README_nonlinear.md from the existing JSON (no fitting)")
    a = p.parse_args(argv)
    a.layers = [s.strip() for s in a.layers.split(",") if s.strip()]
    a.knn_layers = a.layers if a.knn_layers is None else [s.strip() for s in a.knn_layers.split(",") if s.strip()]
    a.mlp_layers = a.layers if a.mlp_layers is None else [s.strip() for s in a.mlp_layers.split(",") if s.strip()]
    a.knn_metrics = [s for s in a.knn_metrics.split(",") if s]
    for m in a.knn_metrics:
        if m not in KNN_METRICS:
            p.error(f"unknown kNN metric {m}")
    a.mlp_seeds = [int(s) for s in a.mlp_seeds.split(",") if s != ""]
    if a.max_train and not a.smoke:
        p.error("--max-train is for smoke runs only")
    return a


# ============================================================================ probes

def fit_mlp(Xs, Y, tr, va, ev, seed, a) -> tuple[np.ndarray, dict]:
    import torch
    from torch import nn
    torch.set_num_threads(_THREADS)
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    d, h = Xs.shape[1], a.mlp_hidden
    net = nn.Sequential(nn.Linear(d, h), nn.ReLU(), nn.Dropout(a.mlp_dropout),
                        nn.Linear(h, h), nn.ReLU(), nn.Dropout(a.mlp_dropout), nn.Linear(h, Y.shape[1]))
    opt = torch.optim.AdamW(net.parameters(), lr=a.mlp_lr, weight_decay=a.mlp_wd)
    Xtr, Ytr = torch.from_numpy(Xs[tr]), torch.from_numpy(Y[tr])
    Xva, Yva = torch.from_numpy(Xs[va]), torch.from_numpy(Y[va])
    best_mse, best_state, best_ep, bad, hist = np.inf, None, -1, 0, []
    t0 = time.time()
    for ep in range(a.mlp_max_epochs):
        net.train()
        perm = torch.from_numpy(rng.permutation(len(tr)))
        for s in range(0, len(tr), a.mlp_batch):
            ix = perm[s:s + a.mlp_batch]
            loss = ((net(Xtr[ix]) - Ytr[ix]) ** 2).mean()
            opt.zero_grad(set_to_none=True)
            loss.backward()
            opt.step()
        net.eval()
        with torch.no_grad():
            mse = float(((net(Xva).clamp(0, 1) - Yva) ** 2).mean())
        hist.append(round(mse, 6))
        if mse < best_mse - 1e-6:
            best_mse, best_ep, bad = mse, ep, 0
            best_state = {k: v.detach().clone() for k, v in net.state_dict().items()}
        else:
            bad += 1
            if bad >= a.mlp_patience:
                break
    net.load_state_dict(best_state)
    net.eval()
    with torch.no_grad():
        pred = net(torch.from_numpy(Xs[ev])).clamp(0, 1).numpy().astype(np.float32)
        pva = net(Xva).clamp(0, 1).numpy()
    info = {"seed": seed, "epochs_run": len(hist), "best_epoch": best_ep + 1, "stopped_early": len(hist) < a.mlp_max_epochs,
            "val_mse_mean": rnd(best_mse, 6), "val_mse": {p: rnd(((pva[:, j] - Y[va, j]) ** 2).mean(), 6) for j, p in enumerate(PROTEINS)},
            "val_mse_curve": hist, "fit_sec": round(time.time() - t0, 1)}
    return pred, info


def knn_neighbours(F: np.ndarray, tr: np.ndarray, queries: np.ndarray, kmax: int, block: int = 1024) -> np.ndarray:
    """Sorted indices (into tr) of the kmax most cosine-similar training cells for each query row. F is
    row-normalised; query cells are never training cells (val / site4), so no self-match exists."""
    Ftr = np.ascontiguousarray(F[tr])
    out = np.empty((len(queries), kmax), np.int32)
    for s in range(0, len(queries), block):
        e = min(s + block, len(queries))
        S = F[queries[s:e]] @ Ftr.T
        part = np.argpartition(-S, kmax - 1, axis=1)[:, :kmax]
        sub = np.take_along_axis(S, part, 1)
        out[s:e] = np.take_along_axis(part, np.argsort(-sub, axis=1, kind="stable"), 1)
    return out


def knn_probe(X, Xs, metric, Y, rows, tr) -> tuple[np.ndarray, dict, dict]:
    F = Xs if metric == "cos_std" else X
    F = (F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-8)).astype(np.float32)
    kmax = max(KNN_GRID)
    t0 = time.time()
    nva = knn_neighbours(F, tr, rows["val"], kmax)
    nev = knn_neighbours(F, tr, rows["eval"], kmax)
    Ytr = Y[tr]
    cva = np.cumsum(Ytr[nva], axis=1, dtype=np.float64)  # [n_val, kmax, 4]
    ks, mses = [], {}
    for j, p in enumerate(PROTEINS):
        m = {k: float(((cva[:, k - 1, j] / k - Y[rows["val"], j]) ** 2).mean()) for k in KNN_GRID}
        mses[p] = {str(k): rnd(v, 6) for k, v in m.items()}
        ks.append(min(m, key=m.get))
    cev = np.cumsum(Ytr[nev], axis=1, dtype=np.float64)
    pred = np.stack([cev[:, k - 1, j] / k for j, k in enumerate(ks)], 1).astype(np.float32)
    info = {"k": dict(zip(PROTEINS, ks)), "k_at_grid_edge": [p for p, k in zip(PROTEINS, ks) if k in (min(KNN_GRID), max(KNN_GRID))],
            "val_mse_by_k": mses, "fit_sec": round(time.time() - t0, 1)}
    return pred, info, {"nev": nev, "k": ks}


def neighbour_overlap(nev, ks, pos_eval, pairs) -> dict:
    """Share of the k training neighbours two cells have in common, median over pairs (per protein's k)."""
    out = {}
    for v in VARIANTS:
        nk, t = pairs["neighbour"][v]
        o = pairs["random"][v][0]
        blk = {}
        for p, k in zip(PROTEINS, ks):
            def share(a, b):
                A, B = nev[pos_eval[a], :k], nev[pos_eval[b], :k]
                return np.array([np.intersect1d(x, y, assume_unique=True).size / k for x, y in zip(A, B)])
            sh = share(nk, t) if nk.size else np.zeros(0)
            r = min(2000, o[0].size)
            shr = share(o[0][:r], o[1][:r]) if r else np.zeros(0)
            blk[p] = {"k": int(k), "lookalike_median_shared": rnd(np.median(sh)) if sh.size else None,
                      "lookalike_mean_shared": rnd(sh.mean()) if sh.size else None,
                      "random_pairs_mean_shared": rnd(shr.mean()) if shr.size else None}
        out[v] = blk
    return out


def paired(pg_a: np.ndarray, pg_b: np.ndarray, mg: np.ndarray, boot: list) -> dict:
    """(median gap a - median gap b) / median measured gap, with the same NK-cell resamples."""
    m0, a0, b0 = np.median(mg, 0), np.median(pg_a, 0), np.median(pg_b, 0)
    bm = np.array([np.median(mg[ix], 0) for ix in boot]).reshape(-1, len(PROTEINS))
    ba = np.array([np.median(pg_a[ix], 0) for ix in boot]).reshape(-1, len(PROTEINS))
    bb = np.array([np.median(pg_b[ix], 0) for ix in boot]).reshape(-1, len(PROTEINS))
    res = {}
    for j, p in enumerate(PROTEINS):
        if m0[j] <= 0:
            res[p] = {"diff": None}
            continue
        ok = bm[:, j] > 0
        bs = (ba[ok, j] - bb[ok, j]) / bm[ok, j]
        res[p] = {"diff": rnd((a0[j] - b0[j]) / m0[j]),
                  "ci95": [rnd(np.percentile(bs, 2.5)), rnd(np.percentile(bs, 97.5))] if bs.size else None}
    return res


# ============================================================================ README

def fmt(x, n=3):
    return "" if x is None else f"{x:.{n}f}"


def cell_ratio(g, p):
    e = (g or {}).get(p) or {}
    ci = e.get("ratio_ci95")
    return fmt(e.get("ratio")) + (f" [{fmt(ci[0])}, {fmt(ci[1])}]" if ci else "")


def cell_diff(e):
    if not e or e.get("diff") is None:
        return ""
    ci = e.get("ci95")
    return f"{e['diff']:+.3f}" + (f" [{ci[0]:+.3f}, {ci[1]:+.3f}]" if ci else "")


def write_readme(out: Path, R: dict, reading: str | None) -> None:
    L, H = R["layers"], R["head_reference"]
    lines = []
    if R["smoke_test"]:
        lines += [f"> **SMOKE TEST. These numbers are not results.** {R['smoke_note']}", ""]
    lines += ["# Mode A step 2: non-linear probes on TEDDY's pooled layers", "",
              f"Written by `scripts/mode_a_nonlinear_probes.py` ({SCRIPT_VERSION}), {R['created']}. Status: **{R['status']}**.", "",
              f"Same cells, targets, splits, look-alike pairs (all {R['pairs']['all']['n_pairs']}, no_gdT158 "
              f"{R['pairs']['no_gdT158']['n_pairs']}), bootstrap and head as `README.md` (the linear run); checks: pairs identical to "
              f"the linear run = {R['checks'].get('pairs_identical_to_linear_run')}, head gap ratios reproduced (max abs diff "
              f"{fmt(R['checks'].get('head_ratio_max_abs_diff_vs_linear_run'), 4)}).", "",
              f"MLP: 2 x {R['config']['mlp']['hidden']} ReLU, dropout {R['config']['mlp']['dropout']}, AdamW lr "
              f"{R['config']['mlp']['lr']}, early stopping on val (patience {R['config']['mlp']['patience']}, max "
              f"{R['config']['mlp']['max_epochs']} epochs), seeds {R['config']['mlp']['seeds']}; 'MLP' rows = mean of the seeds' "
              f"predictions, per-seed ratios in the JSON. kNN: k per protein on val from {list(KNN_GRID)}; cos_std = cosine on "
              "train-standardised features, cos_raw = cosine on raw features (at z: the pair-selection geometry).", ""]
    # Pearson
    lines += ["## site4 Pearson (p95-scaled, clipped)", "", "| depth | read-out | " + " | ".join(PROTEINS) + " |",
              "|---|---|" + "---|" * len(PROTEINS)]
    for l in L:
        for nm, blk in readouts(l):
            lines.append(f"| {l['label']} | {nm} | " + " | ".join(fmt((blk.get('site4_pearson') or {}).get(p)) for p in PROTEINS) + " |")
    lines.append("| - | head | " + " | ".join(fmt(H["site4_pearson"][p]) for p in PROTEINS) + " |")
    for v in VARIANTS:
        lines += ["", f"## Gap ratio on the look-alike pairs: {v} ({R['pairs'][v]['n_pairs']} pairs)", "",
                  "Median |pred NK - pred T| / median |measured NK - measured T|, [95% NK-cell bootstrap]; 1 = gap kept, 0 = gap gone. "
                  "Last column: the same read-out on random NK x T pairs (median over 5 seeds), CD56 / CD94 / CD335 / CD3.", "",
                  "| depth | read-out | " + " | ".join(PROTEINS) + " | random pairs |", "|---|---|" + "---|" * (len(PROTEINS) + 1)]
        for l in L:
            for nm, blk in readouts(l):
                g = (blk.get("gap") or {}).get(v) or {}
                rr = g.get("random_pairs_ratio_median_over_seeds") or {}
                lines.append(f"| {l['label']} | {nm} | " + " | ".join(cell_ratio(g, p) for p in PROTEINS) +
                             " | " + " / ".join(fmt(rr.get(p), 2) for p in PROTEINS) + " |")
        g = H["gap"][v]
        rr = g.get("random_pairs_ratio_median_over_seeds") or {}
        lines.append("| - | head | " + " | ".join(cell_ratio(g, p) for p in PROTEINS) + " | " + " / ".join(fmt(rr.get(p), 2) for p in PROTEINS) + " |")
        lines += ["", "MLP per-seed range of the gap ratio: " + "; ".join(
            f"{l['label']}: " + ", ".join(f"{p} {fmt(l['mlp']['seed_range'][v][p][0])}-{fmt(l['mlp']['seed_range'][v][p][1])}" for p in PROTEINS)
            for l in L if l.get("mlp")) + "."]
        lines += ["", f"### Paired differences in gap ratio ({v}); interval above 0 = the first read-out keeps more of the gap", "",
                  "| depth | comparison | " + " | ".join(PROTEINS) + " |", "|---|---|" + "---|" * len(PROTEINS)]
        for l in L:
            for cname, c in l.get("comparisons", {}).items():
                lines.append(f"| {l['label']} | {cname} | " + " | ".join(cell_diff(c.get(v, {}).get(p)) for p in PROTEINS) + " |")
    # look-alike ratio relative to the same read-out's random-pair ratio
    lines += ["", "## Look-alike gap ratio divided by the same read-out's random-pair ratio", "",
              "A read-out that shrinks every prediction toward the mean lowers both ratios; this quotient asks how much MORE the "
              "look-alike gap is compressed than a random NK x T gap by the same read-out (1 = no look-alike-specific loss). "
              "Point values only (derived from the two tables above).", "",
              "| pair set | depth | read-out | " + " | ".join(PROTEINS) + " |", "|---|---|---|" + "---|" * len(PROTEINS)]
    for v in VARIANTS:
        for l in L:
            for nm, blk in readouts(l):
                g = (blk.get("gap") or {}).get(v) or {}
                lines.append(f"| {v} | {l['label']} | {nm} | " + " | ".join(fmt(_norm(g, p)) for p in PROTEINS) + " |")
        lines.append(f"| {v} | - | head | " + " | ".join(fmt(_norm(H['gap'][v], p)) for p in PROTEINS) + " |")
    # kNN bias
    lines += ["", "## kNN construction bias: shared training neighbours", "",
              "Mean share of the k nearest training cells that the two cells of a pair have in common (k = the protein's chosen k; "
              "CD56 shown, other proteins in the JSON). High for look-alike pairs at z by construction.", "",
              "| depth | metric | k (CD56) | look-alike all | look-alike no_gdT158 | random NK x T |", "|---|---|---|---|---|---|"]
    for l in L:
        for m, blk in (l.get("knn") or {}).items():
            o = blk["neighbour_overlap"]
            lines.append(f"| {l['label']} | {m} | {o['all']['CD56']['k']} | {fmt(o['all']['CD56']['lookalike_mean_shared'])} | "
                         f"{fmt(o['no_gdT158']['CD56']['lookalike_mean_shared'])} | {fmt(o['all']['CD56']['random_pairs_mean_shared'])} |")
    edge = [(l["label"], m, b["k_at_grid_edge"]) for l in L for m, b in (l.get("knn") or {}).items() if b["k_at_grid_edge"]]
    if edge:
        lines += ["", "kNN k at the edge of its grid: " + "; ".join(f"{a} {m}: {', '.join(e)}" for a, m, e in edge) + "."]
    if reading:
        lines += ["", reading.strip()]
    (out / "README_nonlinear.md").write_text("\n".join(lines) + "\n")


def _norm(g, p):
    r = (g.get(p) or {}).get("ratio")
    q = (g.get("random_pairs_ratio_median_over_seeds") or {}).get(p)
    return r / q if (r is not None and q) else None


def readouts(l):
    r = [("linear (ridge)", l["linear"])]
    if l.get("mlp"):
        r.append(("MLP", l["mlp"]["ensemble"]))
    for m, b in (l.get("knn") or {}).items():
        r.append((f"kNN {m}", b))
    return r


# ============================================================================ main

def main(argv=None) -> int:
    a = parse_args(argv)
    from threadpoolctl import threadpool_limits
    threadpool_limits(_THREADS)
    t_start = time.time()
    out = a.out_dir
    out.mkdir(parents=True, exist_ok=True)
    if a.readme_only:
        R = json.loads((out / "mode_a_nonlinear_results.json").read_text())
        write_readme(out, R, a.reading.read_text() if a.reading and a.reading.exists() else None)
        say(f"rebuilt {out / 'README_nonlinear.md'} from the JSON")
        return 0
    say(f"{SCRIPT_VERSION}: threads={_THREADS} out={out}")
    LR = json.loads((a.linear_dir / "mode_a_results.json").read_text())
    LG = np.load(a.linear_dir / "mode_a_pair_gaps.npz")
    lin_layers = {l["name"]: l for l in LR["layers"]}
    lin_index = {str(n): i for i, n in enumerate(LG["layer_names"])}
    n_boot = a.n_boot or int(LR["config"]["n_boot"])
    cfg = LR["config"]

    need_input = "input" in set(a.mlp_layers) | set(a.knn_layers)
    meta = lp.load_meta(a.processed, need_rna=need_input)
    n_total = len(meta["split"])
    key = np.array([lp.ANNOTATION_MAP[c] or "" for c in meta["cell_types"]])
    ctype, split, names = meta["cell_types"], meta["split"], meta["adt_names"]
    col = {n: j for j, n in enumerate(names)}
    p95 = lp.p95_replicate(meta)
    p95v = np.array([p95[n] for n in names], dtype=np.float32)
    meas4 = np.clip(meta["adt"] / p95v, 0, 1)[:, [col[p] for p in PROTEINS]].astype(np.float32)

    emb = lp.Embedding(a.embed_dir, None, n_total)
    if emb.partial and not a.smoke:
        raise SystemExit("partial embedding: only with --smoke")
    cells = emb.cells
    keyp, ctp, spp, Y = key[cells], ctype[cells], split[cells], meas4[cells]
    rows = {"train": np.where(spp == "train")[0], "val": np.where(spp == "val")[0], "eval": np.where(spp == "test")[0]}
    tr = rows["train"]
    if a.max_train:
        tr = np.sort(np.random.default_rng(0).choice(tr, min(a.max_train, tr.size), replace=False))
    z_raw = np.asarray(emb.z, dtype=np.float32)
    pairs = lp.build_pairs(z_raw[rows["eval"]], rows["eval"], keyp, ctp, int(cfg["k"]), int(cfg["n_random_per_seed"]),
                           list(cfg["random_seeds"]))
    boot = {v: lp.cluster_boot_indices(pairs["neighbour"][v][0], n_boot, int(cfg["seed"])) for v in VARIANTS}
    same = all(np.array_equal(LG[f"nk_{v}"], cells[pairs["neighbour"][v][0]]) and
               np.array_equal(LG[f"t_{v}"], cells[pairs["neighbour"][v][1]]) for v in VARIANTS)
    if not same and not a.smoke:
        raise SystemExit("rebuilt look-alike pairs differ from the linear run's; refusing to compare")
    say(f"pairs: all {pairs['neighbour']['all'][0].size}, no_gdT158 {pairs['neighbour']['no_gdT158'][0].size}; identical to linear run: {same}")
    pos_eval = np.full(len(cells), -1, np.int64)
    pos_eval[rows["eval"]] = np.arange(rows["eval"].size)
    meas_gap = {v: np.abs(Y[pairs["neighbour"][v][0]] - Y[pairs["neighbour"][v][1]]) for v in VARIANTS}

    # head reference (per-pair gaps for the paired comparisons; must reproduce the linear run)
    hargs = argparse.Namespace(head_ckpt=a.head_ckpt, head_pred=None)
    pred, _hinfo = lp.head_predictions(hargs, emb, meta)
    pred4 = np.clip(pred / p95v, 0, 1)[:, [col[p] for p in PROTEINS]]
    hfull = np.full_like(pred4, np.nan)
    hfull[rows["eval"]] = pred4[rows["eval"]]
    hgap, hcache = lp.gap_block(hfull, Y, pairs, boot)
    Lh = LR["head_reference"]["gap"]
    hdiff = max(abs((hgap[v][p]["ratio"] or 0) - (Lh[v][p]["ratio"] or 0)) for v in VARIANTS for p in PROTEINS)
    say(f"head reference reproduced: max |ratio diff| = {hdiff:.4f}")
    checks = {"pairs_identical_to_linear_run": bool(same), "head_ratio_max_abs_diff_vs_linear_run": rnd(hdiff, 6),
              "n_boot": n_boot}
    head_ref = {"source": str(a.head_ckpt), "site4_pearson": LR["head_reference"]["site4_pearson"], "gap": hgap}

    depths = []
    for d in a.layers + [x for x in a.mlp_layers + a.knn_layers if x not in a.layers]:
        if d not in depths:
            depths.append(d)
    layers_out = []
    cnt = lambda r, c: int(np.sum(keyp[r] == c))  # noqa: E731

    def dump(status):
        R = {"smoke_test": bool(a.smoke), "smoke_note": a.smoke, "status": status, "created": time.strftime("%Y-%m-%d %H:%M:%S"),
             "config": {"script_version": SCRIPT_VERSION, "threads": _THREADS, "device": "cpu", "layers": a.layers,
                        "mlp_layers": a.mlp_layers, "knn_layers": a.knn_layers, "knn_metrics": a.knn_metrics, "knn_grid": KNN_GRID,
                        "mlp": {"hidden": a.mlp_hidden, "layers": 2, "activation": "ReLU", "dropout": a.mlp_dropout, "lr": a.mlp_lr,
                                "weight_decay": a.mlp_wd, "batch": a.mlp_batch, "max_epochs": a.mlp_max_epochs,
                                "patience": a.mlp_patience, "seeds": a.mlp_seeds, "loss": "MSE on the p95-scaled clipped target, 4 proteins jointly",
                                "early_stopping": "val split (training site), mean MSE of the clipped prediction"},
                        "reused_from_linear_run": {k: cfg[k] for k in ("k", "n_random_per_seed", "random_seeds", "seed", "target_scale", "pairs")},
                        "n_train_used": int(tr.size), "features": "each depth standardised on train cells (as the ridge); z = L2-normalised z_rna.npy"},
             "inputs": {"processed": str(a.processed), "embed_dir": str(a.embed_dir), "linear_dir": str(a.linear_dir),
                        "head_ckpt": str(a.head_ckpt), "ckpt": str(a.ckpt) if a.ckpt else None},
             "cells": {"n_train": int(rows["train"].size), "n_val": int(rows["val"].size), "n_eval_site4": int(rows["eval"].size),
                       "eval_NK": cnt(rows["eval"], "NK"), "eval_T": cnt(rows["eval"], "T")},
             "pairs": {v: {"n_pairs": int(pairs["neighbour"][v][0].size)} for v in VARIANTS},
             "checks": checks, "head_reference": head_ref, "layers": layers_out,
             "runtime_sec": round(time.time() - t_start, 1)}
        (out / "mode_a_nonlinear_results.json").write_text(json.dumps(R, indent=1))
        reading = a.reading.read_text() if a.reading and a.reading.exists() else None
        write_readme(out, R, reading)

    for d in depths:
        tl = time.time()
        name = {"input": "input", "z": "z"}.get(d, f"layer_{d}")
        label = {"input": "in"}.get(d, d)
        if name == "input":
            X = load_input_layer(a, meta, cells, emb)
        elif name == "z":
            X = emb.z_head_input()
        else:
            X = emb.layer(int(d) - 1)
        Xs = lp.standardise(X, tr)
        L0 = lin_layers[name]
        lin = {"site4_pearson": L0["protein_probe"]["site4_pearson"], "gap": L0["protein_probe"]["gap"], "alpha": L0["protein_probe"]["alpha"]}
        lin_pg = {v: LG[f"pred_gap_{v}"][lin_index[name]] for v in VARIANTS}
        entry = {"name": name, "label": label, "linear": lin, "comparisons": {}, "timing_sec": {}}
        ev = rows["eval"]
        if d in a.mlp_layers:
            preds, infos, per_seed = [], [], []
            for sd in a.mlp_seeds:
                p_, info = fit_mlp(Xs, Y, tr, rows["val"], ev, sd, a)
                full = np.full((len(cells), len(PROTEINS)), np.nan, np.float32)
                full[ev] = p_
                g, _ = lp.gap_block(full, Y, pairs, boot)
                info["site4_pearson"] = {p: rnd(pearson(p_[:, j], Y[ev, j])) for j, p in enumerate(PROTEINS)}
                info["gap"] = g
                per_seed.append(info)
                preds.append(p_)
                say(f"  {name} MLP seed {sd}: epochs {info['epochs_run']} (best {info['best_epoch']}) val {info['val_mse_mean']} "
                    f"r " + " ".join(f"{p}={info['site4_pearson'][p]}" for p in PROTEINS) + " | nogd ratio " +
                    " ".join(f"{p}={g['no_gdT158'][p]['ratio']}" for p in PROTEINS) + f" ({info['fit_sec']}s)")
            pe = np.mean(preds, 0).astype(np.float32)
            full = np.full((len(cells), len(PROTEINS)), np.nan, np.float32)
            full[ev] = pe
            g, gc = lp.gap_block(full, Y, pairs, boot)
            ens = {"site4_pearson": {p: rnd(pearson(pe[:, j], Y[ev, j])) for j, p in enumerate(PROTEINS)}, "gap": g}
            srange = {v: {p: [min(s["gap"][v][p]["ratio"] for s in per_seed), max(s["gap"][v][p]["ratio"] for s in per_seed)]
                          for p in PROTEINS} for v in VARIANTS}
            entry["mlp"] = {"ensemble": ens, "per_seed": per_seed, "seed_range": srange}
            entry["comparisons"]["MLP - linear"] = {v: paired(gc[v], lin_pg[v], meas_gap[v], boot[v]) for v in VARIANTS}
            entry["comparisons"]["MLP - head"] = {v: paired(gc[v], hcache[v], meas_gap[v], boot[v]) for v in VARIANTS}
            entry["timing_sec"]["mlp"] = round(sum(s["fit_sec"] for s in per_seed), 1)
        if d in a.knn_layers:
            entry["knn"] = {}
            for m in a.knn_metrics:
                p_, info, nb = knn_probe(X, Xs, m, Y, rows, tr)
                full = np.full((len(cells), len(PROTEINS)), np.nan, np.float32)
                full[ev] = p_
                g, gc = lp.gap_block(full, Y, pairs, boot)
                info["site4_pearson"] = {p: rnd(pearson(p_[:, j], Y[ev, j])) for j, p in enumerate(PROTEINS)}
                info["gap"] = g
                info["neighbour_overlap"] = neighbour_overlap(nb["nev"], nb["k"], pos_eval, pairs)
                entry["knn"][m] = info
                entry["comparisons"][f"kNN {m} - linear"] = {v: paired(gc[v], lin_pg[v], meas_gap[v], boot[v]) for v in VARIANTS}
                entry["comparisons"][f"kNN {m} - head"] = {v: paired(gc[v], hcache[v], meas_gap[v], boot[v]) for v in VARIANTS}
                say(f"  {name} kNN {m}: k {info['k']} r " + " ".join(f"{p}={info['site4_pearson'][p]}" for p in PROTEINS) +
                    " | nogd ratio " + " ".join(f"{p}={g['no_gdT158'][p]['ratio']}" for p in PROTEINS) + f" ({info['fit_sec']}s)")
        entry["timing_sec"]["total"] = round(time.time() - tl, 1)
        layers_out.append(entry)
        del X, Xs
        dump(f"partial ({len(layers_out)} of {len(depths)} depths)")
        say(f"{name} done ({time.time() - tl:.0f}s)")
    dump("complete")
    say(f"wrote {out / 'mode_a_nonlinear_results.json'} in {time.time() - t_start:.0f}s")
    return 0


def load_input_layer(a, meta, cells, emb) -> np.ndarray:
    if a.ckpt is None or a.medians is None:
        raise SystemExit("the input layer needs --ckpt and --medians")
    seq_len = int(emb.manifest.get("seq_len", 2048))
    norm = float(emb.manifest.get("normalize_total", 1e4))
    sig = hashlib.sha256(json.dumps([lp.file_sig(a.ckpt / "model.safetensors"), lp.file_sig(a.medians),
                                     lp.file_sig(a.processed / "cite_arrays.npz"), seq_len, norm, int(cells.size),
                                     int(cells.sum())], default=str).encode()).hexdigest()[:16]
    if a.cache_dir is not None:
        cp = a.cache_dir / f"input_layer_genemean_{sig}.npy"
        if cp.exists():
            say(f"input layer from cache {cp}")
            return np.load(cp)
    X, info = lp.input_layer_means(meta, cells, a.ckpt, a.medians, seq_len, norm, emb.ntokens)
    say(f"input layer recomputed: {info}")
    if a.cache_dir is not None:
        a.cache_dir.mkdir(parents=True, exist_ok=True)
        np.save(cp, X)
    return X


if __name__ == "__main__":
    sys.exit(main())
