"""E3 readouts on frozen TEDDY: R1 (MLP on the 12 gene-mean layer outputs) and R2 (one learned query
attending over layer-12 gene-token states). Architectures and training settings are the registered ones
(registration_v3.json experiments.E3.readouts) with the open details fixed in registration/addenda/E3.json.

R2 with its query at 0 attends uniformly, i.e. it starts as the gene-mean (TEDDY's own pooling) followed by
a linear map; training moves it away from mean pooling only where that lowers the validation error.
"""
from __future__ import annotations

import copy
import math
import time
from typing import Callable, Iterable

import numpy as np
import torch
import torch.nn as nn


class R1MLP(nn.Module):
    """Standardise (fixed training mean / SD), then 2 hidden layers of 256 (ReLU, dropout 0.1), linear output."""

    def __init__(self, d_in: int, d_out: int, hidden: int = 256, dropout: float = 0.1):
        super().__init__()
        self.register_buffer("mu", torch.zeros(d_in))
        self.register_buffer("sd", torch.ones(d_in))
        self.net = nn.Sequential(nn.Linear(d_in, hidden), nn.ReLU(), nn.Dropout(dropout),
                                 nn.Linear(hidden, hidden), nn.ReLU(), nn.Dropout(dropout),
                                 nn.Linear(hidden, d_out))

    def set_standardisation(self, mu: np.ndarray, sd: np.ndarray) -> None:
        self.mu.copy_(torch.as_tensor(mu, dtype=torch.float32))
        self.sd.copy_(torch.as_tensor(np.maximum(sd, 1e-6), dtype=torch.float32))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net((x - self.mu) / self.sd)


class R2AttnPool(nn.Module):
    """One learned query q over the token states of a cell; linear output.

    h~_t = (h_t - mu) / sd (fixed per-dimension affine from the training cells' gene-means; a
    reparametrisation, since a convex combination commutes with it and softmax ignores the constant
    shift), s_t = <q, h~_t> / sqrt(d) over real tokens, a = softmax(s), y = W (sum_t a_t h~_t) + b.
    q starts at 0 (uniform attention = gene-mean)."""

    def __init__(self, d: int, d_out: int):
        super().__init__()
        self.d = int(d)
        self.register_buffer("mu", torch.zeros(d))
        self.register_buffer("sd", torch.ones(d))
        self.q = nn.Parameter(torch.zeros(d))
        self.out = nn.Linear(d, d_out)

    def set_standardisation(self, mu: np.ndarray, sd: np.ndarray) -> None:
        self.mu.copy_(torch.as_tensor(mu, dtype=torch.float32))
        self.sd.copy_(torch.as_tensor(np.maximum(sd, 1e-6), dtype=torch.float32))

    def forward(self, h: torch.Tensor, mask: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        """h: [B, L, d] float; mask: [B, L] bool (True = real token). Returns (y [B, d_out], a [B, L])."""
        m = mask.bool()
        ht = ((h - self.mu) / self.sd).masked_fill(~m.unsqueeze(-1), 0.0)
        s = (ht @ self.q) / math.sqrt(self.d)
        s = s.masked_fill(~m, float("-inf"))
        a = torch.softmax(s, dim=1)
        a = torch.nan_to_num(a, nan=0.0)  # a cell without tokens (never in this data) gets a = 0
        pooled = (a.unsqueeze(-1) * ht).sum(dim=1)
        return self.out(pooled), a


def fit_with_early_stopping(model: nn.Module, train_batches: Callable[[int], Iterable[tuple]],
                            loss_on: Callable[[nn.Module, tuple], torch.Tensor],
                            val_mse: Callable[[nn.Module], float], *, lr: float, max_epochs: int, patience: int,
                            log: Callable[[str], None] = print) -> dict:
    """Adam(lr) on MSE; after each epoch the validation MSE; stop after ``patience`` epochs without a new
    best; restore the best weights. train_batches(epoch) yields batches for loss_on."""
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    best, best_state, best_epoch, bad = float("inf"), None, 0, 0
    hist = []
    for ep in range(1, int(max_epochs) + 1):
        t0 = time.time()
        model.train()
        tl, nb = 0.0, 0
        for batch in train_batches(ep):
            loss = loss_on(model, batch)
            opt.zero_grad()
            loss.backward()
            opt.step()
            tl += float(loss.detach())
            nb += 1
        model.eval()
        with torch.no_grad():
            v = float(val_mse(model))
        hist.append({"epoch": ep, "train_mse": tl / max(nb, 1), "val_mse": v, "sec": time.time() - t0})
        improved = v < best - 1e-12
        if improved:
            best, best_state, best_epoch, bad = v, copy.deepcopy(model.state_dict()), ep, 0
        else:
            bad += 1
        log(f"  epoch {ep}: train MSE {tl / max(nb, 1):.6f} val MSE {v:.6f}{' (best)' if improved else ''} "
            f"[{time.time() - t0:.0f}s]")
        if bad >= int(patience):
            break
    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    return {"history": hist, "best_epoch": best_epoch, "best_val_mse": best, "epochs_run": len(hist)}


def pad_states(chunks: list[np.ndarray]) -> tuple[torch.Tensor, torch.Tensor]:
    """Stack variable-length [L_i, d] float16 states into ([B, Lmax, d] float32, [B, Lmax] bool mask)."""
    lmax = max(1, max(c.shape[0] for c in chunks))
    d = chunks[0].shape[1]
    H = np.zeros((len(chunks), lmax, d), dtype=np.float32)
    M = np.zeros((len(chunks), lmax), dtype=bool)
    for i, c in enumerate(chunks):
        H[i, :c.shape[0]] = c
        M[i, :c.shape[0]] = True
    return torch.from_numpy(H), torch.from_numpy(M)
