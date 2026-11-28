"""Correcting the head (run 52): constants, the halves, the cells and the head. DESIGN.md has the rule."""
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
RUNS = RUN.parent
sys.path.insert(0, str(RUNS / "batch3_frozen_20261126" / "scripts"))
import _b3 as B3  # noqa: E402

STUDY = "frangieh2021"
CONDITIONS = ("Co-culture", "Control")
OUTPUT = {"PDL1": "CD274", "MHC1": "HLA-A-B-C", "HLAE": "HLA-E"}      # the loop's readout -> the row of the head's last layer
READOUTS = tuple(OUTPUT)
LAST = "dec__mu.net.6"                                                # the only layer that moves
SEED, RESAMPLES = 20261128, 2000
ARMS = ("H0", "H1", "H2", "H3")
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def by_path(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def teach(name):
    """The teach half (True) or the judge half of a control cell, by its name alone."""
    return int(hashlib.sha256(("head|" + str(name)).encode()).hexdigest(), 16) % 2 == 0


def cells(cells_file):
    """Every melanoma cell of the run, by condition: ids, unit, role, and for a control cell its half."""
    t = B3.table(cells_file)
    out = {}
    for c in CONDITIONS:
        m = (t["study"] == STUDY) & (t["condition"] == c)
        ids, unit, role = t["cell"][m], t["target"][m], t["role"][m]
        out[c] = {"ids": ids, "unit": unit, "role": role, "control": unit == "control",
                  "drawn": role == "drawn", "teach": np.array([teach(i) for i in ids])}
        assert out[c]["drawn"].sum() and not (out[c]["drawn"] & out[c]["control"]).any()
    return out


def drawn_units(C):
    return {c: sorted(set(C[c]["unit"][C[c]["drawn"]].tolist())) for c in CONDITIONS}


def recorded_units(C, least=100):
    """The events recorded in these two conditions: read before this run, with at least `least` cells."""
    out = {}
    for c in CONDITIONS:
        g = C[c]
        keep = ~g["drawn"] & ~g["control"]
        out[c] = sorted(u for u in set(g["unit"][keep].tolist()) if int(((g["unit"] == u) & keep).sum()) >= least)
    return out


def head(states, ids):
    """The activation the last layer sees (512 per cell), the head's output as shipped, and the shipped rows."""
    TG = by_path("teg_common", RUNS / "teg_training_20261119" / "scripts" / "_common.py")
    P, names, _j = TG.weights()
    x = TG.unit(states)
    h = x.astype(np.float64)
    for name, act in zip(TG.LAYERS, TG.ACT):
        if name == LAST:
            break
        h = h @ P[name + ".weight"].T.astype(np.float64) + P[name + ".bias"]
        h = h * (1.0 / (1.0 + np.exp(-h))) if act else h
    rows = {r: (P[LAST + ".weight"][names.index(n)].astype(np.float64), float(P[LAST + ".bias"][names.index(n)])) for r, n in OUTPUT.items()}
    assert h.shape == (len(ids), 512)
    return h, rows


def out_of(a, w, b):
    """The head's output for one protein: log of softplus, as the head as released computes it."""
    return np.log(np.maximum(np.logaddexp(0.0, a @ w + b), 1e-6))
