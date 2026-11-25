#!/usr/bin/env python3
"""The cells of batch 2: the cells of its units, and the cells already read that its readings stand on. Record only.

  python cells.py FRANGIEH_PROTEIN.h5ad PAPALEXI_PROTEIN.h5ad roles.json known_cells.tsv OUT.tsv

Of the two protein files only the record of the cells is opened (names, guides, condition): no protein value. For every
unit of roles.json, the cells that carry guides for that gene only, in its condition (the rule of the census). known_cells.tsv
lists the cells earlier runs read (the development cells of the two strata and the cells of batch 1) with their role;
they are written through unchanged. The script stops if a cell of the batch is among them.
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


def main():
    roles = json.load(open(sys.argv[3]))
    known = [l.rstrip("\n").split("\t") for l in open(sys.argv[4]).read().splitlines()[1:] if l.strip()]
    seen = {r[3] for r in known}
    units = {}
    for u in roles["units"]:
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
    f = h5py.File(sys.argv[2], "r")
    o = f["obs"]
    names = column(o, "_index").astype(str)
    tg = np.array([target_papalexi(g) for g in column(o, "guide_id").astype(str)])
    for g in units.get("papalexi2021 | all", []):
        rows += [("papalexi2021", "all", g, n, "event") for n in names[tg == g]]
    new = [r[3] for r in rows]
    assert len(set(new)) == len(new), "a cell is listed twice"
    assert not (set(new) & seen), "a cell of the batch was read before"
    assert all(r[2] != "control" for r in rows)
    with open(sys.argv[5], "w") as out:
        out.write("study\tcondition\ttarget\tcell\trole\n")
        for r in known + rows:
            out.write("\t".join(r) + "\n")
    count = {}
    for r in rows:
        count[(r[0], r[1], r[2])] = count.get((r[0], r[1], r[2]), 0) + 1
    print(len(rows), "cells of the batch;", len(known), "cells read before")
    for k, v in sorted(count.items()):
        print(" ", *k, v)
    assert len(count) == len(roles["units"]) and min(count.values()) >= 100, "a unit of the batch has fewer than 100 cells or none"


if __name__ == "__main__":
    main()
