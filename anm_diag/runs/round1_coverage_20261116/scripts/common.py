"""What the scripts of the round share: the head run from its weights in numpy, the proteins both studies have, scores."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RUN.parents[1]))            # the anm_diag package
sys.path.insert(0, str(Path(__file__).resolve().parent))

SF = 0.9266302585601807                            # the training-median size factor scripts/05_eval.py uses; no score here depends on it
ALPHAS = [float(a) for a in np.logspace(-1, 5, 13)]  # the ridge penalties of runs/scope_20261116


def silu(a):
    return a / (1.0 + np.exp(-a))


def head_forward(W, z_unit, chunk=20000):
    """The head as scored: unit-length state -> MLP -> decoder -> softplus, times the constant size factor. Any width."""
    out = []
    W = {k: np.asarray(v, np.float64) for k, v in W.items() if k != "adt_names"}
    for i in range(0, len(z_unit), chunk):
        h = np.asarray(z_unit[i:i + chunk], np.float64)
        for part in ("mlp__net", "dec__mu.net"):
            h = silu(h @ W[f"{part}.0.weight"].T + W[f"{part}.0.bias"])
            h = silu(h @ W[f"{part}.3.weight"].T + W[f"{part}.3.bias"])
            h = h @ W[f"{part}.6.weight"].T + W[f"{part}.6.bias"]
        out.append((SF * np.logaddexp(0.0, h)).astype(np.float32))
    return np.vstack(out)


def load_studies(bm_arrays, bm_z, pbmc_arrays, pbmc_z):
    """Both studies as the trainer reads them, and the columns of the proteins the second study measured."""
    import train_arm
    bm, pb = train_arm.load_data(bm_arrays, bm_z, pbmc_arrays, pbmc_z)
    shared = np.where(~np.isnan(pb["adt"]).all(axis=0))[0]
    if np.isnan(pb["adt"][:, shared]).any():
        raise SystemExit("a protein of the second study is measured in some cells only")
    return bm, pb, shared


def score(Y, P):
    """The head's ordinary score: the mean over proteins of the correlation between measured and predicted."""
    Y, P = np.asarray(Y, np.float64), np.asarray(P, np.float64)
    r = [np.corrcoef(Y[:, j], P[:, j])[0, 1] for j in range(Y.shape[1]) if Y[:, j].std() > 1e-8 and P[:, j].std() > 1e-8]
    return float(np.mean(r)) if r else float("nan")


def per_protein(Y, P):
    Y, P = np.asarray(Y, np.float64), np.asarray(P, np.float64)
    return np.array([np.corrcoef(Y[:, j], P[:, j])[0, 1] if Y[:, j].std() > 1e-8 and P[:, j].std() > 1e-8 else np.nan for j in range(Y.shape[1])])
