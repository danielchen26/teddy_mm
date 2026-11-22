#!/usr/bin/env python3
"""The cells of batch 1: for which cells the system's states are to be computed next (census.json, batch_1).

  python batch_cells.py FRANGIEH_PROTEIN.h5ad PAPALEXI_PROTEIN.h5ad census.json coculture_control_cells_with_states.txt OUT.tsv

Record only: cell names, guides, design fields. For every unit of the batch, the cells that carry guides for that gene
only. For every stratum of the batch, the control cells that have no states yet and are not final
(sha256("final|<cell name>") not divisible by 3). THP-1 needs no new control cells: all of them have states.
"""
import hashlib
import json
import sys

import h5py
import numpy as np

import re

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


final = lambda cell: int(hashlib.sha256(f"final|{cell}".encode()).hexdigest(), 16) % 3 == 0
census = json.load(open(sys.argv[3]))
have = set(open(sys.argv[4]).read().split())
units = {}
for u in census["batch_1"]["units"]:
    gene, stratum = u.split("@")
    units.setdefault(stratum, []).append(gene)
rows = []
f = h5py.File(sys.argv[1], "r")
o = f["obs"]
names, cond = column(o, "cell_name").astype(str), column(o, "perturbation_2").astype(str)
tg = np.array([target_frangieh(g) for g in column(o, "guide_id").astype(str)])
for stratum, genes in units.items():
    study, c = [x.strip() for x in stratum.split("|")]
    if study != "frangieh2021":
        continue
    for g in genes:
        rows += [(study, c, g, n, "event") for n in names[(cond == c) & (tg == g)]]
    rows += [(study, c, "control", n, "control") for n in names[(cond == c) & (tg == "control")] if n not in have and not final(n)]
f = h5py.File(sys.argv[2], "r")
o = f["obs"]
names = column(o, "_index").astype(str)
tg = np.array([target_papalexi(g) for g in column(o, "guide_id").astype(str)])
for g in units.get("papalexi2021 | all", []):
    rows += [("papalexi2021", "all", g, n, "event") for n in names[tg == g]]
with open(sys.argv[5], "w") as out:
    out.write("study\tcondition\ttarget\tcell\trole\n")
    for r in rows:
        out.write("\t".join(r) + "\n")
count = {}
for r in rows:
    count[(r[0], r[1], r[2])] = count.get((r[0], r[1], r[2]), 0) + 1
print(len(rows), "cells")
for k, v in sorted(count.items()):
    print(" ", *k, v)
