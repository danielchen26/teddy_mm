"""What the scripts of batch 2 share: the strata, the edges known beforehand, and the cells of a stratum."""
import json
import sys
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
RUNS = RUN.parent
sys.path.insert(0, str(RUNS.parent))
sys.path.insert(0, str(RUNS / "teg_training_20261119" / "scripts"))
from anm_diag import account as AC  # noqa: E402,F401
from anm_diag import cycle  # noqa: E402,F401
from _common import SEED, unit  # noqa: E402,F401

ROLES = json.loads((RUN / "results" / "roles.json").read_text())
# stratum -> (study, condition, readouts)
STRATA = {"mel_co": ("frangieh2021", "Co-culture", ["PDL1", "MHC1", "HLAE"]), "mel_ctl": ("frangieh2021", "Control", ["PDL1", "MHC1", "HLAE"]), "thp1": ("papalexi2021", "all", ["PDL1"])}
# the edges known beforehand, as runs/account_20261125/scripts/library.py: the event's site is a gene of the measured protein
OWN = {"PDL1": ["CD274"], "MHC1": ["HLA-A", "HLA-B", "HLA-C", "B2M"], "HLAE": ["HLA-E", "B2M"]}
EDGES = {(e, r): "own" for r, genes in OWN.items() for e in genes}
MIN_ROWS = 8          # pairs without an edge a stratum needs to tell the spread of a direct part from its own rows


def table(path):
    rows = [l.rstrip("\n").split("\t") for l in Path(path).read_text().splitlines()[1:] if l.strip()]
    return {k: np.array([r[i] for r in rows]) for i, k in enumerate(("study", "condition", "target", "cell", "role"))}


def stratum(g, cells, states, read=None, new=None):
    """The cells of a stratum in the order of the list: ids, unit (the event of a recorded cell, 'control', or 'new:<gene>'
    for a cell of the batch), role, Z at unit length, the blocks of the typed state, and the readouts: of the cells read
    before from `read`, of the cells of the batch from `new` (the file made after the forecasts were committed), else NaN."""
    study, cond, readouts = STRATA[g]
    m = (cells["study"] == study) & (cells["condition"] == cond)
    ids, role, tg = cells["cell"][m], cells["role"][m], cells["target"][m]
    S = np.load(Path(states) / f"batch2_states_{study}.npz")
    pos = {i: k for k, i in enumerate(S["ids"].astype(str))}
    k = np.array([pos[i] for i in ids])
    c = dict(ids=ids, role=role, unit=np.where(role == "event", np.char.add("new:", tg), tg), raw=S["z_last"][k].astype(np.float64), Z=unit(S["z_last"][k]).astype(np.float64),
             sets=[str(s) for s in S["sets"]], blocks=S["blocks"][k], readouts=readouts)
    for r in readouts:
        y = np.full(len(ids), np.nan)
        for src in (read, new):
            if src is None:
                continue
            d = np.load(src)
            p = {i: j for j, i in enumerate(d["ids"].astype(str))}
            have = np.array([i in p for i in ids])
            y[have] = d[f"readout__{r}"][np.array([p[i] for i in ids[have]])]
        c[f"y_{r}"] = y
        assert np.isfinite(y[role != "event"]).all(), "a cell read before has no readout"
        assert (new is not None) or np.isnan(y[role == "event"]).all(), "a readout of a cell of the batch is here before its file was made"
    return c
