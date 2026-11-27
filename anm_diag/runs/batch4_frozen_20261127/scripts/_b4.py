"""Batch 4: the cells, their states and their readouts, for the scripts of this run. One screen, two strata, 20 readouts."""
import json
import sys
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
RUNS = RUN.parent
sys.path.insert(0, str(RUNS / "v6_candidate_20261126" / "scripts"))
sys.path.insert(0, str(RUNS / "batch3_frozen_20261126" / "scripts"))
import _v6 as V6  # noqa: E402
from _v6 import V  # noqa: E402
import _b3 as B3  # noqa: E402  (the list of cells, the states and the antibody columns are read as in batch 3)

STUDY = "frangieh2021"
STRATA = {"mel_co": "Co-culture", "mel_ctl": "Control"}
READOUTS = list(V.MELANOMA)
ROLES = json.loads((RUN / "results" / "roles.json").read_text())
table, States, Packs = B3.table, B3.States, B3.Packs


def chosen():
    """{stratum: {unit: how it was chosen}} from roles.json."""
    out = {}
    for g, cond in STRATA.items():
        v = ROLES["conditions"][f"{STUDY} | {cond}"]
        out[g] = {**{u.split("@")[0]: "by an item on record" for u in v["by_an_item_on_record"]}, **{u.split("@")[0]: "by lot" for u in v["by_lot"]}}
    return out


def stratum(g, cells, states):
    """The cells of one stratum. kind: control; recorded (a unit read in batches 1 to 3 or before); new (a unit of this batch)."""
    m = (cells["study"] == STUDY) & (cells["condition"] == STRATA[g])
    ids, role, unit = cells["cell"][m], cells["role"][m], cells["target"][m]
    kind = np.where(role == "control", "control", np.where(role == "event", "new", "recorded"))
    assert set(unit[kind == "new"].tolist()) == set(chosen()[g]), "the units of the list are not those of roles.json"
    return {"ids": ids, "role": role, "unit": unit, "kind": kind, "raw": states.rows(ids), "Z": V.unit_length(states.rows(ids)), "states": states}


def edges_on_record():
    """The edges proposed from data so far, with what became of them: run 43, then batch 3."""
    r43 = [e for e in json.loads((RUNS / "one_loop_20261125" / "results" / "retro.json").read_text())["edges"] if e["status"] != "declared"]
    b3 = json.loads((RUNS / "batch3_frozen_20261126" / "results" / "reading.json").read_text())["edges"]
    return r43 + b3["tested"], b3["with_no_row_in_stage_B"]
