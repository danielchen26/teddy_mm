#!/usr/bin/env python3
"""How the answer key relates to the dataset's cell-type annotations (site4 held-out cells).

The answer key used everywhere is a 3-lineage rule on the measured proteins (B: CD19 CD72 CD22,
T: CD3 CD2 CD5, myeloid: CD16 CD11c CD36; top lineage if it leads by the question's margin).
It has no NK class. This script cross-tabulates that key (soft rule) against the annotated
coarse cell types, e.g. how many NK cells the key calls "myeloid" (NK/ILC coarse class, and NK cells
alone under each question: O0 soft, O1 strict, O2 key-marker).

Usage (from the teddy_mm repo root):
    python scripts/answer_key_vs_annotation.py
"""
from __future__ import annotations

import collections
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from bridge_anm.lib.lineage_panels import CRITERIA, expected_from_true  # noqa: E402
from bridge_anm.reverse_step1_must_separate import _coarse_from_cell_type  # noqa: E402

BRIDGE = ROOT / "outputs/anm_cite_bridge"


def main() -> None:
    p95 = {}
    with (BRIDGE / "cite_typed_events.jsonl").open() as f:
        for line in f:
            ev = json.loads(line)
            p95.setdefault(ev["protein"], ev["norm_p95_train"])
            if len(p95) == 9:
                break
    tab = collections.Counter()
    nk_by_question = collections.Counter()  # NK cells only (annotated "NK" or "NK CD158e1+"), keyed myeloid per question
    n_nk = 0
    with (BRIDGE / "cite_cells_meta.jsonl").open() as f:
        for line in f:
            cell = json.loads(line)
            if cell.get("site") != "site4":
                continue
            key = expected_from_true(cell["adt_true_panel_holdout"], p95, CRITERIA["O0"]) or "no clear answer"
            tab[(key, _coarse_from_cell_type(cell["cell_type"]))] += 1
            if cell["cell_type"] in ("NK", "NK CD158e1+"):
                n_nk += 1
                for q in ("O0", "O1", "O2"):
                    if expected_from_true(cell["adt_true_panel_holdout"], p95, CRITERIA[q]) == "myeloid":
                        nk_by_question[q] += 1
    keys = ["b_lineage", "t_lineage", "myeloid", "no clear answer"]
    coarse = sorted({c for _, c in tab})
    nk = {k: tab[(k, "nk_ilc")] for k in keys}
    keyed_myeloid = sum(tab[("myeloid", c)] for c in coarse)
    print(json.dumps({
        "crosstab_key_x_annotation": {k: {c: tab[(k, c)] for c in coarse} for k in keys},
        "nk_ilc_cells": sum(nk.values()),
        "nk_ilc_keyed_as": nk,
        "keyed_myeloid": keyed_myeloid,
        "keyed_myeloid_annotated_myeloid": tab[("myeloid", "myeloid")],
        "nk_cells_only": n_nk,
        "nk_cells_keyed_myeloid_by_question": dict(nk_by_question),
    }, indent=1))


if __name__ == "__main__":
    main()
