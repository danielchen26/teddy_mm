#!/usr/bin/env python3
"""The units of the two screens, from the record only: for every study, condition and target gene, the number of cells
that carry guides for that one gene (plus any non-targeting guides), and the number of control cells. No measurement is
read: only the cell names, the detected guides and the fields of the design.

  python record_census.py FRANGIEH_PROTEIN.h5ad PAPALEXI_PROTEIN.h5ad OUT.json

The rule for a cell's target is the one of runs/teg_repair_test_20261119/scripts/cells.py (Frangieh) and of the
declarations of runs/split_20261111 (Papalexi: the guide's gene; NT is the control).
"""
import json
import re
import sys

import h5py
import numpy as np

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


def target_papalexi(guide):
    m = re.match(r"^(.+?)g\d+$", guide)
    g = m.group(1) if m else guide
    return "control" if g in ("NT", "nan", "") else g


def table(target, strata):
    out = {}
    for s in sorted(set(strata)):
        m = strata == s
        t, c = np.unique(target[m & (target != "")], return_counts=True)
        out[s] = {"cells": int(m.sum()), "cells_without_one_target": int((m & (target == "")).sum()), "units": {str(a): int(b) for a, b in zip(t, c)}}
    return out


f = h5py.File(sys.argv[1], "r")
o = f["obs"]
tg = np.array([target_frangieh(g) for g in column(o, "guide_id").astype(str)])
out = {"frangieh2021": {"fields": {"cell line": "melanoma (one line)", "condition": sorted(set(column(o, "perturbation_2").astype(str)))},
                        "by_condition": table(tg, column(o, "perturbation_2").astype(str))}}
f = h5py.File(sys.argv[2], "r")
o = f["obs"]
tg = np.array([target_papalexi(g) for g in column(o, "guide_id").astype(str)])
hto = column(o, "hto").astype(str)
out["papalexi2021"] = {"fields": {"cell line": "THP-1", "replicate and treatment (hto)": sorted(set(hto))}, "by_replicate": table(tg, hto),
                       "all": table(tg, np.full(len(tg), "all"))["all"]}
json.dump(out, open(sys.argv[3], "w"), indent=0)
print({k: {s: (v["cells"], len(v["units"])) for s, v in x.get("by_condition", x.get("by_replicate")).items()} for k, x in out.items()})
