"""Shared cosine kNN for Mode B evidence-geometry probes."""
from __future__ import annotations

import numpy as np


def l2_normalize(z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, dtype=np.float32)
    norms = np.linalg.norm(z, axis=1, keepdims=True) + 1e-8
    return z / norms


def knn_cosine(z: np.ndarray, k: int, batch: int = 256) -> tuple[np.ndarray, np.ndarray]:
    zn = l2_normalize(z)
    n = zn.shape[0]
    k_eff = min(k, max(n - 1, 1))
    idx = np.zeros((n, k_eff), dtype=np.int64)
    sims = np.zeros((n, k_eff), dtype=np.float32)
    if n <= 1:
        return idx, sims
    for s in range(0, n, batch):
        e = min(s + batch, n)
        S = zn[s:e] @ zn.T
        for i in range(e - s):
            S[i, s + i] = -np.inf
        part = np.argpartition(-S, kth=k_eff - 1, axis=1)[:, :k_eff]
        row = np.take_along_axis(S, part, axis=1)
        order = np.argsort(-row, axis=1)
        idx[s:e] = np.take_along_axis(part, order, axis=1)
        sims[s:e] = np.take_along_axis(row, order, axis=1)
    return idx, sims
