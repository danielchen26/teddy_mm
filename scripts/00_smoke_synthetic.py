#!/usr/bin/env python3
"""Build a tiny processed CITE pack + fake z_rna so training can be tested without GEO."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    out = ROOT / "data/processed/cite_smoke"
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    n, g, p, zdim = 240, 80, 16, 512
    rna = sparse.csr_matrix(rng.poisson(1.5, size=(n, g)).astype(np.float32))
    programs = rng.normal(size=(n, 4)).astype(np.float32)
    adt = np.clip(programs @ rng.normal(size=(4, p)).astype(np.float32) + rng.normal(0, 0.3, size=(n, p)), 0, None)
    adt = rng.poisson(adt + 0.2).astype(np.float32)
    split = np.array(["train"] * 160 + ["val"] * 40 + ["test"] * 40, dtype="U8")
    z = np.concatenate([programs, rng.normal(scale=0.1, size=(n, zdim - 4))], axis=1).astype(np.float32)
    np.savez_compressed(
        out / "cite_arrays.npz",
        rna_data=rna.data,
        rna_indices=rna.indices,
        rna_indptr=rna.indptr,
        rna_shape=np.array(rna.shape),
        adt=adt,
        token_ids=np.arange(g, dtype=np.int64),
        split=split,
        donors=np.array([f"d{i%6}" for i in range(n)], dtype="U32"),
        sites=np.array(["site1"] * 200 + ["site4"] * 40, dtype="U32"),
        cell_types=np.array(["t"] * n, dtype="U64"),
        rna_size_factor=np.ones(n, dtype=np.float32),
        adt_size_factor=np.ones(n, dtype=np.float32),
        adt_names=np.array([f"CD{i}" for i in range(p)]),
        rna_names=np.array([f"ENSG{i:011d}" for i in range(g)]),
    )
    np.save(out / "z_rna.npy", z)
    meta = {"n_cells": n, "n_adt": p, "n_train": 160, "n_val": 40, "n_test": 40, "smoke": True}
    (out / "meta.json").write_text(json.dumps(meta, indent=2))
    print("smoke pack ->", out)


if __name__ == "__main__":
    main()
