"""Batch 3: the cells, their states and their readouts, for the scripts of this run."""
import json
import sys
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RUN.parent / "freeze_v5_20261126" / "scripts"))
import _v5 as V  # noqa: E402

ROLES = json.loads((RUN / "results" / "roles.json").read_text())
STRATA = {"mel_co": ("frangieh2021", "Co-culture"), "mel_ctl": ("frangieh2021", "Control"), "thp1": ("papalexi2021", "all")}
CONDITION = {"library": "IFNγ", "mel_co": "Co-culture", "mel_ctl": "Control"}
PANEL = {"frangieh2021": V.MELANOMA, "papalexi2021": V.THP1}
USED = {"frangieh2021": V.USED["melanoma"], "papalexi2021": V.USED["thp1"]}
NEW = {s: [r for r in PANEL[s] if r not in USED[s]] for s in PANEL}
NOTHING = "(nothing recorded)"                                # the placeholder event of a stratum with no recorded row: a copy of its control cells


def table(path):
    rows = [l.rstrip("\n").split("\t") for l in Path(path).read_text().splitlines()[1:] if l.strip()]
    return {k: np.array([r[i] for r in rows]) for i, k in enumerate(("study", "condition", "target", "cell", "role"))}


class States:
    """The states of the listed files, by cell: the first file that holds a cell is the one used."""

    def __init__(self, files):
        self.parts = [np.load(f) for f in files]
        self.pos = {}
        for k, p in enumerate(self.parts):
            for j, i in enumerate(p["ids"].astype(str)):
                self.pos.setdefault(i, (k, j))
        self.sets = [str(s) for s in self.parts[0]["sets"]]

    def rows(self, ids, key="z_last"):
        out = None
        where = np.array([self.pos[i] for i in ids])
        for k, p in enumerate(self.parts):
            m = where[:, 0] == k
            if m.any():
                v = p[key][where[m, 1]]
                if out is None:
                    out = np.zeros((len(ids),) + v.shape[1:], np.float64)
                out[m] = v
        return out


def teachers(g):
    return ROLES["teachers"].get(" | ".join(STRATA[g]), {}).get("teachers", [])


def kind(g, unit, role):
    """What a cell is in this batch: control; teacher (opened in stage A); read (a unit read before, forecast on the new
    readouts); new (one of the seven units)."""
    return np.where(role == "control", "control", np.where(role == "event", "new", np.where((role == "recorded") & np.isin(unit, teachers(g)), "teacher", "read")))


def stratum(g, cells, states):
    study, cond = STRATA[g]
    m = (cells["study"] == study) & (cells["condition"] == cond)
    ids, role, unit = cells["cell"][m], cells["role"][m], cells["target"][m]
    return {"ids": ids, "role": role, "unit": unit, "kind": kind(g, unit, role), "raw": states.rows(ids), "Z": V.unit_length(states.rows(ids)), "states": states}


class Packs:
    """The antibody columns of the exported stages, by cell, and the readouts formed from them."""

    def __init__(self, files):
        self.log = {"frangieh2021": {}, "papalexi2021": {}}
        self.panel = {}
        for f in files:
            P = np.load(f)
            for study, k in (("frangieh2021", "mel"), ("papalexi2021", "thp")):
                self.panel[study] = [str(x) for x in P[f"{k}_panel"]]
                for i, row in zip(P[f"{k}_ids"].astype(str), P[f"{k}_log"].astype(np.float64)):
                    assert i not in self.log[study], "a cell is in two stages"
                    self.log[study][i] = row

    def has(self, study, ids):
        return np.array([i in self.log[study] for i in ids])

    def readouts(self, study, ids, names=None):
        """{readout: values} for the listed cells, with the log count of every readout's antibody and the background where
        the panel has isotype antibodies."""
        L = np.array([self.log[study][i] for i in ids])
        panel = self.panel[study]
        col = {r: panel.index(ab) for r, (ab, _genes) in PANEL[study].items()}
        names = names or list(col)
        if study == "frangieh2021":
            assert len(panel) == 24
            clr = L - L.mean(1, keepdims=True)
            return {r: clr[:, col[r]] for r in names}, {r: L[:, col[r]] for r in names}, L[:, [panel.index(k) for k in V.ISOTYPES]].mean(1)
        return {r: L[:, col[r]] for r in names}, {r: L[:, col[r]] for r in names}, None


def used_values(read, pack2):
    """The readouts read before, by cell: {readout: {cell: value}} (the files of batch 2)."""
    out = {r: {} for r in ("PDL1", "MHC1", "HLAE")}
    for f in (read, pack2):
        d = np.load(f)
        ids = d["ids"].astype(str)
        for r in out:
            out[r].update(zip(ids.tolist(), d[f"readout__{r}"].tolist()))
    return out


def with_placeholder(x, Y, unit):
    """A stratum with no recorded event: the law is built with a placeholder event that is a copy of the control cells. Its
    change and its shift are zero exactly, it has no type: nothing teaches and nothing is taught."""
    c = unit == "control"
    return np.vstack([x, x[c]]), {r: np.concatenate([y, y[c]]) for r, y in Y.items()}, np.concatenate([unit, np.full(int(c.sum()), NOTHING)])
