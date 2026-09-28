#!/usr/bin/env python3
"""Re-score the saved phase-2 checkpoint on the site4 test cells, two ways.

The phase-2 training script (and the Block 2 export) score predict_adt(..., use_fm=True): one random
flow-matching sample per cell. The model also has a direct decode from the fused latent
(use_fm=False), which is the path its main loss trains. This prints the per-protein mean Pearson for
both, per input mode, next to the recorded metrics.

Usage (from the teddy_mm repo root; needs torch):
    python scripts/phase2_decode_check.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from teddy_mm.bidirectional import BidirectionalCite  # noqa: E402
from teddy_mm.data import load_prepared  # noqa: E402


def pearson_mean(a: np.ndarray, b: np.ndarray) -> float:
    s = [np.corrcoef(a[:, j], b[:, j])[0, 1] for j in range(a.shape[1]) if a[:, j].std() > 1e-8 and b[:, j].std() > 1e-8]
    return float(np.mean(s))


def main() -> None:
    dev = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    proc = ROOT / "data/processed/cite"
    pack = load_prepared(proc)
    z = np.load(proc / "z_rna.npy").astype(np.float32)
    z = z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)  # same L2 step as scripts/06_train_bidirectional.py
    m = pack["split"] == "test"
    zt = torch.from_numpy(z[m]).to(dev)
    y = torch.from_numpy(pack["adt"][m]).to(dev)
    sf = torch.from_numpy(pack["adt_size_factor"][m]).to(dev)
    ck = torch.load(ROOT / "outputs/cite_phase2/best.pt", map_location=dev, weights_only=False)
    model = BidirectionalCite(ck["z_dim"], ck["n_adt"]).to(dev)
    model.load_state_dict(ck["model"])
    model.eval()
    Y = y.cpu().numpy()
    out = {"test_cells": int(m.sum()), "recorded": json.loads((ROOT / "outputs/cite_phase2/metrics.json").read_text())["test"], "by_mode": {}}
    with torch.no_grad():
        for mode in ("rna_only", "adt_only", "joint"):
            torch.manual_seed(0)
            one = model.predict_adt(zt, y, sf, mode=mode, fm_steps=20, use_fm=True).cpu().numpy()
            five = np.mean([model.predict_adt(zt, y, sf, mode=mode, fm_steps=20, use_fm=True).cpu().numpy() for _ in range(5)], axis=0)
            direct = model.predict_adt(zt, y, sf, mode=mode, use_fm=False).cpu().numpy()
            out["by_mode"][mode] = {"fm_one_sample": round(pearson_mean(Y, one), 4), "fm_mean_of_5": round(pearson_mean(Y, five), 4),
                                    "direct_decode": round(pearson_mean(Y, direct), 4)}
    (ROOT / "outputs/cite_phase2/decode_check.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
