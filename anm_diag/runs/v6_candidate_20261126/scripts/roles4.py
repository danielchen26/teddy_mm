#!/usr/bin/env python3
"""The units a fourth batch would open, chosen by a written rule. results/batch4_roles.json. Nothing is opened here.

  python roles4.py

Until now the author chose the units of every batch. The rule:

1. Eligible: a unit that is unread, not final, and in none of the three batches read so far.
2. By an item on record: every eligible unit whose knocked-out gene (a) is a gene of a protein of the antibody panel
   (a declared edge), or (b) is the site of an edge on record that was proposed and has had no row to test it.
3. The rest by lot: every condition is filled to PER_CONDITION units from its other eligible units, taken in
   alphabetical order and permuted with seed SEED.

Names only: the census of runs/outer_step_20261122 and the roles of the batches. No value of any unit is read.
"""
import json
import sys
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
RUNS = RUN.parent
sys.path.insert(0, str(RUNS / "freeze_v5_20261126" / "scripts"))
import _v5 as V  # noqa: E402

PER_CONDITION, SEED = 15, 20261127


def main():
    census = json.loads((RUNS / "outer_step_20261122" / "results" / "census.json").read_text())
    final = set(census["roles"]["final"])
    before = set(census["batch_1"]["units"])
    for run in ("batch2_forecast_20261125", "batch3_frozen_20261126"):
        before |= set(json.loads((RUNS / run / "results" / "roles.json").read_text())["units"])
    eligible = sorted(u for u in census["roles"]["unread_not_final"] if u not in before)
    assert not set(eligible) & final
    own = {g for _r, (_ab, genes) in V.MELANOMA.items() for g in genes} | {g for genes in V.PARTNERS.values() for g in genes}
    tested = json.loads((RUNS / "batch3_frozen_20261126" / "results" / "reading.json").read_text())["edges"]
    waiting = {e["site"]: e["edge_id"] for e in tested["with_no_row_in_stage_B"]}
    out = {"rule": {"eligible": "unread, not final, in none of the batches read so far", "by_an_item_on_record": "the gene is a gene of a protein of the panel, or the site of an edge proposed with no row to test it yet",
                    "by_lot": f"each condition filled to {PER_CONDITION} units from its other eligible units, in alphabetical order, permuted with seed {SEED}"},
           "eligible": len(eligible), "conditions": {}, "units": []}
    rng = np.random.default_rng(SEED)
    for cond in sorted({u.split("@")[1] for u in eligible}):
        pool = [u for u in eligible if u.split("@")[1] == cond]
        item = {u: ("a declared edge" if u.split("@")[0] in own else f"the proposed edge {waiting[u.split('@')[0]]}") for u in pool if u.split("@")[0] in own or u.split("@")[0] in waiting}
        rest = [u for u in pool if u not in item]
        lot = sorted(np.array(rest)[rng.permutation(len(rest))[:PER_CONDITION - len(item)]].tolist())
        out["conditions"][cond] = {"eligible": len(pool), "by_an_item_on_record": item, "by_lot": lot, "left_unread_and_not_final": len(pool) - len(item) - len(lot)}
        out["units"] += sorted(item) + lot
    out["left_unread_and_not_final"] = len(eligible) - len(out["units"])
    (RUN / "results" / "batch4_roles.json").write_text(json.dumps(out, indent=1) + "\n")
    for cond, v in out["conditions"].items():
        print(cond, "| by an item:", v["by_an_item_on_record"], "| by lot:", [u.split("@")[0] for u in v["by_lot"]])
    print(len(out["units"]), "units;", out["left_unread_and_not_final"], "left unread and not final")


if __name__ == "__main__":
    main()
