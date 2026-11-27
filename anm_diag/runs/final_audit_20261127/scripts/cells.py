#!/usr/bin/env python3
"""The cells of the final audit, on the owner's machine. Record only: names, guides, condition. No protein value.

  python cells.py FRANGIEH_PROTEIN.h5ad PAPALEXI_PROTEIN.h5ad roles.json KNOWN_CELLS.tsv OUT.tsv

KNOWN_CELLS.tsv: the list of batch 4 (every cell read so far, with its role; the cells of batch 4 as `batch4`).
Added: the cells of the 133 final units (role `event`), by the rule of the census, and the final control cells of the
melanoma strata (role `reference`): the control cells no run has read. The script stops if one of them was read
before, or if a final control cell does not meet the rule of the census (sha256('final|<name>') divisible by 3).
"""
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import h5py
import numpy as np

spec = importlib.util.spec_from_file_location("cells2", Path(__file__).resolve().parents[2] / "batch2_forecast_20261125" / "scripts" / "cells.py")
C = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C)
final = lambda n: int(hashlib.sha256(("final|" + n).encode()).hexdigest(), 16) % 3 == 0


def main():
    roles = json.load(open(sys.argv[3]))
    known = [l.rstrip("\n").split("\t") for l in open(sys.argv[4]).read().splitlines()[1:] if l.strip()]
    assert not any(r[4] in ("event", "reference") for r in known), "the known list still holds cells marked as new"
    seen = {r[3] for r in known}
    units = {}
    for u in roles["units"]:
        gene, stratum = u.split("@")
        units.setdefault(stratum, []).append(gene)
    rows = []
    f = h5py.File(sys.argv[1], "r")
    o = f["obs"]
    names, cond = C.column(o, "cell_name").astype(str), C.column(o, "perturbation_2").astype(str)
    tg = np.array([C.target_frangieh(g) for g in C.column(o, "guide_id").astype(str)])
    for stratum, genes in units.items():
        study, c = [x.strip() for x in stratum.split("|")]
        if study != "frangieh2021":
            continue
        for g in genes:
            rows += [(study, c, g, n, "event") for n in names[(cond == c) & (tg == g)]]
        ref = [n for n in names[(cond == c) & (tg == "control")] if n not in seen]
        assert ref and all(final(n) for n in ref), "an unread control cell is not a final control cell by the rule of the census"
        rows += [(study, c, "control", n, "reference") for n in ref]
    f = h5py.File(sys.argv[2], "r")
    o = f["obs"]
    names = C.column(o, "_index").astype(str)
    tg = np.array([C.target_papalexi(g) for g in C.column(o, "guide_id").astype(str)])
    for g in units.get("papalexi2021 | all", []):
        rows += [("papalexi2021", "all", g, n, "event") for n in names[tg == g]]
    assert not [n for n in names[tg == "control"] if n not in seen], "THP-1 has control cells no run has read: they would be its reference"
    new = [r[3] for r in rows]
    assert len(set(new)) == len(new), "a cell is listed twice"
    assert not (set(new) & seen), "a cell of the final audit was read before"
    with open(sys.argv[5], "w") as out:
        out.write("study\tcondition\ttarget\tcell\trole\n")
        for r in known + rows:
            out.write("\t".join(r) + "\n")
    count = {}
    for r in rows:
        count[(r[0], r[1], r[4])] = count.get((r[0], r[1], r[4]), 0) + 1
    n_units = len({(r[0], r[1], r[2]) for r in rows if r[4] == "event"})
    small = sorted({(r[0], r[1], r[2]) for r in rows if r[4] == "event" and sum(1 for q in rows if q[:3] == r[:3]) < 100})
    print(len(rows), "cells of the final audit;", len(known), "cells read before;", n_units, "units")
    for k, v in sorted(count.items()):
        print(" ", *k, v)
    assert n_units == len(roles["units"]) and not small, ("a final unit has fewer than 100 cells or none", small)


if __name__ == "__main__":
    main()
