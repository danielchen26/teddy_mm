#!/usr/bin/env python3
"""The units of the run, before any value is read: names and counts only.

  python units.py --cells final_cells.tsv --library library.npz [--ifn-list OUT.tsv]

-> results/units.json: the knockouts of each composition, their cells per condition and half, the control cells.
--ifn-list: the IFN-γ cells whose RNA counts are exported on the owner's machine (the knockouts of the table and the
control cells), in the format of the lists of cells of the batches.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _ka as K  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--library", type=Path, required=True)
    ap.add_argument("--ifn-list", type=Path)
    a = ap.parse_args()
    C = K.cells(a.cells, a.library)
    ko = {t: K.knockouts(C, t) for t in K.TESTS}
    assert set(K.FIVE) <= set(ko["ifn"]) and set(K.FIVE) <= set(ko["co"])
    count = lambda c, g: [int(((c["unit"] == g) & c["A"]).sum()), int(((c["unit"] == g) & ~c["A"]).sum())]
    used = sorted(set(ko["ifn"]) | set(ko["co"]))
    res = {"cells_file_sha256": K.sha(a.cells), "library_ids_sha256": hashlib.sha256("\n".join(C["ifn"]["ids"].tolist()).encode()).hexdigest(),
           "halves": 'sha256("half|" + name), even: A', "the_five": list(K.FIVE), "declared_pairs_left_out": K.declared(),
           "knockouts": {t: {"conditions": [K.CONDITION[k] for k in K.TESTS[t]], "n": len(ko[t]), "names": ko[t]} for t in K.TESTS},
           "cells_in_halves_A_B": {K.CONDITION[k]: {"control": count(C[k], "control"), **{g: count(C[k], g) for g in used if (C[k]["unit"] == g).any()}} for k in C}}
    if a.ifn_list:
        m = np.isin(C["ifn"]["unit"], ko["ifn"] + ["control"])
        rows = ["study\tcondition\ttarget\tcell\trole"] + [f"{K.STUDY}\t{K.CONDITION['ifn']}\t{u}\t{i}\tlibrary" for i, u in zip(C["ifn"]["ids"][m], C["ifn"]["unit"][m])]
        a.ifn_list.write_text("\n".join(rows) + "\n")
        res["ifn_cells_exported"] = {"cells": int(m.sum()), "list_sha256": K.sha(a.ifn_list)}
    (K.RUN / "results" / "units.json").write_text(json.dumps(res, indent=1, ensure_ascii=False) + "\n")
    for t in K.TESTS:
        n = [sum(res["cells_in_halves_A_B"][K.CONDITION[k]][g]) for k in K.TESTS[t] for g in ko[t]]
        print(t, len(ko[t]), "knockouts | cells per knockout and condition", min(n), "to", max(n), "| smallest half", min(min(res["cells_in_halves_A_B"][K.CONDITION[k]][g]) for k in K.TESTS[t] for g in ko[t]))
    print("control cells (A, B):", {k: v["control"] for k, v in res["cells_in_halves_A_B"].items()}, "| IFN-γ cells to export:", res.get("ifn_cells_exported"))


if __name__ == "__main__":
    main()
