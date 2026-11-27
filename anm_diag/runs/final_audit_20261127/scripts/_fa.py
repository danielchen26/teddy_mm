"""The final audit: the cells, their states and their readouts, for the scripts of this run. Three strata."""
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
import _b3 as B3  # noqa: E402  (the list of cells, the states and the antibody columns are read as in batches 3 and 4)

STRATA = {"mel_co": ("frangieh2021", "Co-culture"), "mel_ctl": ("frangieh2021", "Control"), "thp1": ("papalexi2021", "all")}
MEL = ("mel_co", "mel_ctl")
READOUTS = {"mel_co": list(V.MELANOMA), "mel_ctl": list(V.MELANOMA), "thp1": list(V.THP1)}
table, States, Packs = B3.table, B3.States, B3.Packs


def final_units(g, roles=None):
    """The final units of a stratum, from results/roles.json (or the roles given)."""
    roles = roles or json.loads((RUN / "results" / "roles.json").read_text())
    return sorted(u.split("@")[0] for u in roles["units"] if u.split("@")[1] == " | ".join(STRATA[g]))


def stratum(g, cells, states, new=None):
    """The cells of one stratum. kind: control (a control cell read before: the laws stand on these); reference (a final
    control cell); new (a cell of a final unit); recorded (a unit read before). `ref`: the cells a new unit is read
    against, the final control cells where the stratum has any, else its control cells."""
    study, cond = STRATA[g]
    m = (cells["study"] == study) & (cells["condition"] == cond)
    ids, role, unit = cells["cell"][m], cells["role"][m], cells["target"][m]
    kind = np.where(role == "control", "control", np.where(role == "reference", "reference", np.where(role == "event", "new", "recorded")))
    if new is not None:
        assert set(unit[kind == "new"].tolist()) == set(new), "the units of the list are not the final units of the stratum"
    own = bool((kind == "reference").any())
    raw = states[study].rows(ids)
    return {"ids": ids, "role": role, "unit": np.where(kind == "reference", "reference", unit), "kind": kind, "raw": raw, "Z": V.unit_length(raw), "states": states[study],
            "ref": kind == ("reference" if own else "control"), "own_reference": own}


def relative(c, law, m):
    """The states of the cells m as the law is to read them: against the reference of the stratum and not against the
    control cells the law was fitted on. A law reads how far an event's cells stand from its centre; the centre is
    moved to the reference."""
    return c["Z"][m] - (c["Z"][c["ref"]].mean(0) - law._centre) if c["own_reference"] else c["Z"][m]


def against(rows, law, n_ref):
    """The spread of a forecast read against n_ref reference cells and not against the law's own control cells: the error
    of the reference is that of n_ref cells (`Law._new_noise` with n_ref in place of the law's count). The interval
    keeps its width in spreads."""
    for r in rows:
        q = law.r[r["readout"]]
        extra = (1 / n_ref - 1 / law._nc) * q["cv"] / q["sd"] ** 2
        k = (r["interval"][1] - r["interval"][0]) / 2 / r["sd"]
        r["sd"] = float(np.sqrt(r["sd"] ** 2 + extra))
        r["interval"] = [float(r["forecast"] - k * r["sd"]), float(r["forecast"] + k * r["sd"])]
        r["reference_cells"] = int(n_ref)
    return rows


def edges_on_record():
    """Every edge proposed from data so far, with what became of it and where: run 43, batch 3, batch 4."""
    out = [e for e in json.loads((RUNS / "one_loop_20261125" / "results" / "retro.json").read_text())["edges"] if e["status"] != "declared"]
    for run in ("batch3_frozen_20261126", "batch4_frozen_20261127"):
        out += json.loads((RUNS / run / "results" / "reading.json").read_text())["edges"]["tested"]
    return out
