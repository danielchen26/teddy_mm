#!/usr/bin/env python3
"""RNA counts of the listed cells, on the owner's machine. Only the first step of
runs/batch2_forecast_20261125/scripts/states.py (its function `counts`), and nothing after it: no token, no state.

  python counts.py --repo TEDDY_MM --cells ifn_cells.tsv --out OUT_DIR [--states-script PATH]

-> OUT_DIR/counts_ifn.npz: rows, cols, vals, ids, genes, shape, as the count files of the batches.
Prints the number of cells and of entries. No value is printed and no protein is opened.
"""
import argparse
import hashlib
import importlib.util
import time
from pathlib import Path

import numpy as np

STUDY = "frangieh2021"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--states-script", type=Path, default=Path(__file__).resolve().parents[2] / "batch2_forecast_20261125" / "scripts" / "states.py")
    a = ap.parse_args()
    spec = importlib.util.spec_from_file_location("b2_states", a.states_script)
    st = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(st)
    rows = [l.rstrip("\n").split("\t") for l in a.cells.read_text().splitlines()[1:] if l.strip()]
    assert {r[0] for r in rows} == {STUDY}
    ids = [r[3] for r in rows]
    a.out.mkdir(parents=True, exist_ok=True)
    p = a.out / "counts_ifn.npz"
    if p.exists():
        raise SystemExit("counts_ifn.npz exists: it is not made again")
    R, C, V, genes, ncol = st.counts(a.repo / "data" / "raw" / "perturb" / st.FILES[STUDY][0], st.FILES[STUDY][1], ids, a.out / "parts_ifn", time.time() + 1e9)
    np.savez_compressed(a.out / "tmp_counts_ifn.npz", rows=R, cols=C, vals=V, ids=np.array(ids), genes=genes, shape=np.array([len(ids), ncol]))
    (a.out / "tmp_counts_ifn.npz").rename(p)
    print("counts:", len(ids), "cells,", len(V), "entries,", ncol, "genes | list sha256", hashlib.sha256(a.cells.read_bytes()).hexdigest()[:16], flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
