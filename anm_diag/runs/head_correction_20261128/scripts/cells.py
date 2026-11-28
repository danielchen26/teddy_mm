#!/usr/bin/env python3
"""The cells of the 93 drawn units, on the owner's machine. Record only: names, guides, condition. No protein value.

  python cells.py FRANGIEH_PROTEIN.h5ad units.json KNOWN_CELLS.tsv OUT.tsv

KNOWN_CELLS.tsv: the list of the final audit (every cell read so far, with its role). Added: the cells of the drawn
units, role `drawn`. Every control cell has been read already, so none is added; the units are read against the
judge half of them (DESIGN.md). The script stops if a drawn cell was read before or a unit has fewer than 100 cells.
"""
import importlib.util
import json
import sys
from pathlib import Path

import h5py
import numpy as np

spec = importlib.util.spec_from_file_location("cells2", Path(__file__).resolve().parents[2] / "batch2_forecast_20261125" / "scripts" / "cells.py")
C = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C)


def main():
    units = json.load(open(sys.argv[2]))["drawn"]["units"]
    known = [l.rstrip("\n").split("\t") for l in open(sys.argv[3]).read().splitlines()[1:] if l.strip()]
    assert not any(r[4] == "drawn" for r in known), "the known list already holds drawn cells"
    seen = {r[3] for r in known}
    by = {}
    for u in units:
        gene, stratum = u.split("@")
        by.setdefault(stratum, []).append(gene)
    assert set(by) <= {"frangieh2021 | Co-culture", "frangieh2021 | Control"}, sorted(by)
    o = h5py.File(sys.argv[1], "r")["obs"]
    names, cond = C.column(o, "cell_name").astype(str), C.column(o, "perturbation_2").astype(str)
    tg = np.array([C.target_frangieh(g) for g in C.column(o, "guide_id").astype(str)])
    rows, small = [], []
    for stratum, genes in sorted(by.items()):
        study, c = [x.strip() for x in stratum.split("|")]
        for g in sorted(genes):
            cells = names[(cond == c) & (tg == g)]
            if len(cells) < 100:
                small.append((c, g, len(cells)))
            rows += [(study, c, g, n, "drawn") for n in cells]
    new = [r[3] for r in rows]
    assert len(set(new)) == len(new), "a cell is listed twice"
    assert not (set(new) & seen), "a drawn cell was read before"
    assert not small, ("a drawn unit has fewer than 100 cells", small)
    with open(sys.argv[4], "w") as out:
        out.write("study\tcondition\ttarget\tcell\trole\n")
        for r in known + rows:
            out.write("\t".join(r) + "\n")
    n = {}
    for r in rows:
        n[r[1]] = n.get(r[1], 0) + 1
    print(len(rows), "cells of", len(units), "drawn units;", len(known), "read before |", n)


if __name__ == "__main__":
    main()
