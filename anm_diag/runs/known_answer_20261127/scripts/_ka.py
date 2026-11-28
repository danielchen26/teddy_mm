"""The known-answer check (run 51): constants, the halves and the lists of cells. DESIGN.md has the rule."""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
RUNS = RUN.parent
sys.path.insert(0, str(RUNS / "batch3_frozen_20261126" / "scripts"))
import _b3 as B3  # noqa: E402  (the list of cells, the states and the antibody columns are read as in batches 3 and 4)

STUDY = "frangieh2021"
FIVE = ("IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1")           # the known answer
CONDITION = {"ctl": "Control", "ifn": "IFNγ", "co": "Co-culture"}
TESTS = {"ifn": ("ifn", "ctl"), "co": ("co", "ctl")}           # the composition: (the condition with the signal, the one without)
SEED, RESAMPLES = 20261127, 2000
MIN_MEAN, EPS = 0.5, 0.1                                       # a gene: mean counts per 10,000 among control cells; what is added before the log
FOUR = 4.0                                                     # the four errors of every other rule
MIN_HALF = 20                                                  # cells of a knockout in a half of a condition, below which it is left out
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def half(name):
    """Half A (True) or B of a cell, by its name alone."""
    return int(hashlib.sha256(("half|" + str(name)).encode()).hexdigest(), 16) % 2 == 0


def declared():
    """The pairs the graph declares as direct: {knockout: [readouts]} (the frozen graph, as the final audit used it)."""
    out = {}
    for e in json.loads((RUNS / "final_audit_20261127" / "results" / "forecasts.json").read_text())["edges"]["declared"]:
        out.setdefault(e["site"], []).append(e["readout"])
    return out


def cells(cells_file, library):
    """Every cell of the run: {condition: {"ids", "unit"}} with unit "control" or the knockout. The cells with no signal
    given and of the co-culture come from the list of cells read so far; those under IFN-γ from the library."""
    t = B3.table(cells_file)
    out = {}
    for k in ("ctl", "co"):
        m = (t["study"] == STUDY) & (t["condition"] == CONDITION[k])
        out[k] = {"ids": t["cell"][m], "unit": t["target"][m]}
    L = np.load(library)
    assert set(L["stratum"].astype(str).tolist()) == {CONDITION["ifn"]}
    out["ifn"] = {"ids": L["ids"].astype(str), "unit": np.where(L["A"] == 0, "control", L["target"].astype(str)).astype(object).astype(str)}
    for c in out.values():
        assert len(set(c["ids"].tolist())) == len(c["ids"])
        c["A"] = np.array([half(i) for i in c["ids"]])
    return out


def knockouts(C, test):
    """The knockouts of a composition: cells in both of its conditions, at least MIN_HALF in every half."""
    a, b = TESTS[test]
    n = lambda c, g, h: int(((c["unit"] == g) & (c["A"] == h)).sum())
    both = (set(C[a]["unit"].tolist()) & set(C[b]["unit"].tolist())) - {"control"}
    return sorted(g for g in both if min(n(C[k], g, h) for k in (a, b) for h in (True, False)) >= MIN_HALF)
