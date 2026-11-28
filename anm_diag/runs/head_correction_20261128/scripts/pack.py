#!/usr/bin/env python3
"""The measured protein of the 93 drawn units, on the owner's machine. Exported once, after the design, the lot and
the scripts that fit and judge are committed.

  python pack.py --repo TEDDY_MM --cells drawn_cells.tsv --units units.json --out OUT_DIR

The 24 antibody columns, as the log of one plus the counts, of the cells of the drawn units (role `drawn`) and of no
other cell: the control cells were exported long ago. The unit of every cell is derived again from its detected
guides and has to agree with the list.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve()
spec = importlib.util.spec_from_file_location("fa_pack", HERE.parents[2] / "final_audit_20261127" / "scripts" / "pack.py")
FA = importlib.util.module_from_spec(spec)
spec.loader.exec_module(FA)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--units", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    target = a.out / "drawn_pack.npz"
    assert not target.exists(), f"{target.name} exists: the drawn cells are exported once"
    want = set(json.loads(a.units.read_text())["drawn"]["units"])
    rows = [l.rstrip("\n").split("\t") for l in a.cells.read_text().splitlines()[1:] if l.strip()]
    f = h5py.File(a.repo / "data" / "raw" / "perturb" / "FrangiehIzar2021_protein.h5ad", "r")
    names = FA.column(f["obs"], "cell_name").astype(str)
    tg = np.array([FA.target_frangieh(g) for g in FA.column(f["obs"], "guide_id").astype(str)])
    cond = FA.column(f["obs"], "perturbation_2").astype(str)
    L = np.log1p(FA.dense(f))
    pos = {n: i for i, n in enumerate(names)}
    keep, idx, units = {k: [] for k in ("ids", "condition", "unit", "role")}, [], set()
    for study, c, unit, cell, role in rows:
        if role != "drawn":
            continue
        assert study == "frangieh2021"
        k = pos[cell]
        assert tg[k] == unit and cond[k] == c, (cell, tg[k], unit, cond[k], c)
        idx.append(k)
        units.add(f"{unit}@{study} | {c}")
        for key, v in (("ids", cell), ("condition", c), ("unit", unit), ("role", role)):
            keep[key].append(v)
    assert units == want, (len(units), len(want))
    out = {f"mel_{k}": np.array(v) for k, v in keep.items()}
    out.update(mel_log=L[np.array(idx)].astype(np.float32), mel_panel=FA.column(f["var"], f["var"].attrs["_index"]).astype(str),
               thp_panel=np.array(["CD86", "PDL1", "PDL2", "CD366"]), thp_ids=np.array([], dtype="<U8"), thp_log=np.zeros((0, 4), np.float32))
    a.out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out / "tmp_pack.npz", **out)
    (a.out / "tmp_pack.npz").rename(target)
    print(len(idx), "cells of", len(units), "drawn units |", len(out["mel_panel"]), "columns | list sha256",
          hashlib.sha256(a.cells.read_bytes()).hexdigest()[:16], "| pack sha256", hashlib.sha256(target.read_bytes()).hexdigest()[:16])


if __name__ == "__main__":
    main()
