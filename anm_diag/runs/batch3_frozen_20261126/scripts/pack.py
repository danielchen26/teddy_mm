#!/usr/bin/env python3
"""The measured protein for batch 3, on the owner's machine, in two stages.

  python pack.py --repo TEDDY_MM --cells batch3_cells.tsv --library library_ids.txt --roles roles.json --stage A --out OUT_DIR
  python pack.py --repo TEDDY_MM --cells batch3_cells.tsv --library library_ids.txt --roles roles.json --stage B --forecasts forecasts.json --out OUT_DIR

Every antibody column of the two screens, as the log of one plus the counts (24 columns in the melanoma screen, 4 in the
THP-1 screen). The readouts are formed from these columns by the reading scripts.

Stage A, what the forecasts stand on: the control cells of every stratum; the cells of the library under interferon
(development; library_ids.txt); and the cells of the teachers of roles.json (recorded pathway events whose rows teach
their type on the readouts never read).
Stage B, what is forecast: every other cell of the list (the other recorded units, the units of batch 2, the seven
units of batch 3). It refuses to run without the committed file of forecasts, and stops if that file was not made for
this list of cells and for the stage-A file in OUT_DIR.

The unit of every listed cell is derived again from its detected guides and has to agree with the list.
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
CO, TH = "frangieh2021 | Co-culture", "papalexi2021 | all"


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


def target_papalexi(guide):
    m = re.match(r"^(.+?)g\d+$", guide)
    g = m.group(1) if m else guide
    return "control" if g in ("NT", "nan", "") else g


def dense(f):
    X = f["X"]
    return sp.csc_matrix((X["data"][:], X["indices"][:], X["indptr"][:]), shape=tuple(X.attrs["shape"])).toarray().astype(np.float64)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--library", type=Path, required=True)
    ap.add_argument("--roles", type=Path, required=True)
    ap.add_argument("--stage", choices=["A", "B"], required=True)
    ap.add_argument("--forecasts", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    target = a.out / f"batch3_pack_{a.stage}.npz"
    assert not target.exists(), f"{target.name} exists: a stage is exported once"
    if a.stage == "B":
        assert a.forecasts is not None, "stage B needs the committed file of forecasts"
        fc = json.loads(a.forecasts.read_text())
        assert fc["cells_file_sha256"] == sha(a.cells), "the forecasts were made for another list of cells"
        assert fc["pack_a_sha256"] == sha(a.out / "batch3_pack_A.npz"), "the forecasts were made from another stage-A file"
        assert fc["pairs"] > 0
    rows = [l.rstrip("\n").split("\t") for l in a.cells.read_text().splitlines()[1:] if l.strip()]
    teachers = json.loads(a.roles.read_text())["teachers"]
    lib = [l.strip() for l in a.library.read_text().splitlines() if l.strip()]

    def first(study, cond, unit, role):
        return role == "control" or (role == "recorded" and unit in teachers.get(f"{study} | {cond}", {}).get("teachers", []))

    raw, out = a.repo / "data" / "raw" / "perturb", {}
    # melanoma screen
    f = h5py.File(raw / "FrangiehIzar2021_protein.h5ad", "r")
    names = column(f["obs"], "cell_name").astype(str)
    tg = np.array([target_frangieh(g) for g in column(f["obs"], "guide_id").astype(str)])
    cond = column(f["obs"], "perturbation_2").astype(str)
    L = np.log1p(dense(f))
    pos = {n: i for i, n in enumerate(names)}
    keep = {k: [] for k in ("ids", "condition", "unit", "role")}
    idx = []
    for study, c, unit, cell, role in rows:
        if study != "frangieh2021" or first(study, c, unit, role) != (a.stage == "A"):
            continue
        k = pos[cell]
        assert tg[k] == unit and cond[k] == c, (cell, tg[k], unit, cond[k], c)
        idx.append(k)
        for key, v in (("ids", cell), ("condition", c), ("unit", unit), ("role", role)):
            keep[key].append(v)
    if a.stage == "A":
        listed = {r[3] for r in rows}
        for cell in lib:
            assert cell not in listed, "a cell of the library is in the list of the strata"
            k = pos[cell]
            assert cond[k] == "IFNγ", (cell, cond[k])
            idx.append(k)
            for key, v in (("ids", cell), ("condition", cond[k]), ("unit", tg[k]), ("role", "library")):
                keep[key].append(v)
    out.update({f"mel_{k}": np.array(v) for k, v in keep.items()}, mel_log=L[np.array(idx)].astype(np.float32), mel_panel=column(f["var"], f["var"].attrs["_index"]).astype(str))
    # THP-1 screen
    f = h5py.File(raw / "PapalexiSatija2021_eccite_protein.h5ad", "r")
    names = column(f["obs"], "_index").astype(str)
    tg = np.array([target_papalexi(g) for g in column(f["obs"], "guide_id").astype(str)])
    L = np.log1p(dense(f))
    pos = {n: i for i, n in enumerate(names)}
    keep, idx = {k: [] for k in ("ids", "unit", "role")}, []
    for study, c, unit, cell, role in rows:
        if study != "papalexi2021" or first(study, c, unit, role) != (a.stage == "A"):
            continue
        k = pos[cell]
        assert tg[k] == unit, (cell, tg[k], unit)
        idx.append(k)
        for key, v in (("ids", cell), ("unit", unit), ("role", role)):
            keep[key].append(v)
    out.update({f"thp_{k}": np.array(v) for k, v in keep.items()}, thp_log=L[np.array(idx)].astype(np.float32), thp_panel=column(f["var"], f["var"].attrs["_index"]).astype(str))
    a.out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out / f"tmp_pack_{a.stage}.npz", **out)
    (a.out / f"tmp_pack_{a.stage}.npz").rename(target)
    print("stage", a.stage, "| melanoma cells", len(out["mel_ids"]), "| THP-1 cells", len(out["thp_ids"]), "| panels", len(out["mel_panel"]), len(out["thp_panel"]))
    for name in ("mel", "thp"):
        u, n = np.unique(np.array([f"{r}" for r in out[f"{name}_role"]]), return_counts=True)
        print(" ", name, dict(zip(u.tolist(), n.tolist())))
    print("sha256", sha(target))


if __name__ == "__main__":
    main()
