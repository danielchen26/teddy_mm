#!/usr/bin/env python3
"""The measured protein of batch 4, on the owner's machine. Exported once, after the forecasts are committed.

  python pack.py --repo TEDDY_MM --cells batch4_cells.tsv --forecasts forecasts.json --out OUT_DIR

The 24 antibody columns of the melanoma screen, as the log of one plus the counts, for the cells of the 30 units of the
batch (role "event" in the list) and for no other cell. It refuses to run without the committed file of forecasts, and
stops if that file was not made for this list of cells. The unit of every cell is derived again from its detected
guides and has to agree with the list.
"""
import argparse
import hashlib
import json
import re
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

CONTROL_GUIDES = {"NO_SITE", "ONE_NON-GENE_SITE"}


def column(group, key):
    d = group[key]
    if isinstance(d, h5py.Group):
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in d["categories"][:]], dtype=object)
        codes = d["codes"][:]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], "nan")
    v = d[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v], dtype=object) if v.dtype.kind in "SO" else v


def target_frangieh(guides):
    if guides in ("nan", "", "None", "NaN"):
        return ""
    genes = set()
    for g in guides.split(";"):
        m = re.match(r"^(.+)_\d+$", g.strip())
        genes.add(m.group(1) if m else g.strip())
    genes -= CONTROL_GUIDES
    return "control" if not genes else (next(iter(genes)) if len(genes) == 1 else "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--forecasts", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    target = a.out / "batch4_pack.npz"
    assert not target.exists(), f"{target.name} exists: the batch is exported once"
    fc = json.loads(a.forecasts.read_text())
    assert fc["cells_file_sha256"] == sha(a.cells), "the forecasts were made for another list of cells"
    assert fc["pairs"] > 0
    rows = [l.rstrip("\n").split("\t") for l in a.cells.read_text().splitlines()[1:] if l.strip()]
    f = h5py.File(a.repo / "data" / "raw" / "perturb" / "FrangiehIzar2021_protein.h5ad", "r")
    names = column(f["obs"], "cell_name").astype(str)
    tg = np.array([target_frangieh(g) for g in column(f["obs"], "guide_id").astype(str)])
    cond = column(f["obs"], "perturbation_2").astype(str)
    X = f["X"]
    L = np.log1p(sp.csc_matrix((X["data"][:], X["indices"][:], X["indptr"][:]), shape=tuple(X.attrs["shape"])).toarray().astype(np.float64))
    pos = {n: i for i, n in enumerate(names)}
    keep, idx = {k: [] for k in ("ids", "condition", "unit", "role")}, []
    for study, c, unit, cell, role in rows:
        if role != "event":
            continue
        assert study == "frangieh2021", study
        k = pos[cell]
        assert tg[k] == unit and cond[k] == c, (cell, tg[k], unit, cond[k], c)
        idx.append(k)
        for key, v in (("ids", cell), ("condition", c), ("unit", unit), ("role", role)):
            keep[key].append(v)
    units = sorted({f"{u}@{c}" for u, c in zip(keep["unit"], keep["condition"])})
    assert len(units) == fc["units"], (len(units), fc["units"])
    out = {f"mel_{k}": np.array(v) for k, v in keep.items()}
    out.update(mel_log=L[np.array(idx)].astype(np.float32), mel_panel=column(f["var"], f["var"].attrs["_index"]).astype(str),
               thp_ids=np.array([], dtype=str), thp_unit=np.array([], dtype=str), thp_role=np.array([], dtype=str), thp_log=np.zeros((0, 4), np.float32), thp_panel=np.array(["CD86", "PDL1", "PDL2", "CD366"]))
    a.out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out / "tmp_pack.npz", **out)
    (a.out / "tmp_pack.npz").rename(target)
    print("batch 4 | cells", len(out["mel_ids"]), "| units", len(units), "| panel", len(out["mel_panel"]))
    print("sha256", sha(target))


if __name__ == "__main__":
    main()
