"""Shared by the scripts of this run: where things are, the cells and their split, the head, and the reading.

Nothing here is machine-specific. Set these environment variables:

  ANM_FEATURES   exported cell features and head_weights.npz, as for runs/split_20261111 (frangieh2021/, papalexi2021/)
  ANM_BM_FIT     bone-marrow cells for the holding term: arrays z_train, z_val, adt_val, adt_names (as in runs/loop_guard_20261112)
  ANM_BM_TEST    bone-marrow test cells for the guard: arrays z, adt, adt_names
  ANM_TEG_OUT    where the heads of this run are written (default: features/ in this folder, not in git)
  ANM_GUARD_HEADS  optional: the heads of runs/loop_guard_20261112 (<data set>__<gene>__mlp_guarded.npz), to compare the refit of arm `data`
"""
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
ROOT = RUN.parents[1]
LOOP = ROOT / "runs" / "loop_20261112"
RES = RUN / "results"
FEATURES = Path(os.environ["ANM_FEATURES"]) if os.environ.get("ANM_FEATURES") else None
OUT = Path(os.environ.get("ANM_TEG_OUT", RUN / "features"))
BM_FIT, BM_TEST, GUARD_HEADS = os.environ.get("ANM_BM_FIT"), os.environ.get("ANM_BM_TEST"), os.environ.get("ANM_GUARD_HEADS")
sys.path.insert(0, str(ROOT))

TARGETS = {"frangieh2021": ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"], "papalexi2021": ["IFNGR1", "IFNGR2", "JAK2", "STAT1", "IRF1"]}
CELLS = {"frangieh2021": "melanoma", "papalexi2021": "THP-1"}
# the arms: what the supervision of the PD-L1 output carries (DESIGN.md). Everything else is the same in every arm.
ARMS = ["data", "class", "edge1", "edge05", "edge2", "graph", "judge", "rows", "data_rows"]
STANDS_FOR = {"as_trained": "the head as trained", "data": "the repair already made (measured protein, cell by cell)", "class": "an ontology: the class of each cell",
              "edge1": "a knowledge graph: a signed edge, margin 1 SD", "edge05": "a knowledge graph: a signed edge, margin 0.5 SD", "edge2": "a knowledge graph: a signed edge, margin 2 SD",
              "graph": "a KG embedding: the signed edge and the geometry of the graph", "judge": "a decision engine: the judgments of a learned model",
              "rows": "the TEG: its rows (tested comparisons)", "data_rows": "the repair already made, with the rows of the TEG added"}
MLP = ["mlp__net.0", "mlp__net.3", "mlp__net.6"]
DEC = ["dec__mu.net.0", "dec__mu.net.3", "dec__mu.net.6"]
LAYERS, ACT = MLP + DEC, [True, True, False, True, True, False]
LR, STEPS, EVERY, PATIENCE, KEEP, LAM, SEED = 3e-4, 400, 20, 3, 1.0, 1.0, 0      # as in runs/loop_guard_20261112 (its chosen weight of the holding term is 1)
MARGIN, FALL, B, ALPHA = 0.10, 0.02, 1000, 0.05                                 # the tool's margin, the guard of the step-4 rule, bootstrap draws
F = np.float32


def halves(groups, rng, share=0.5):
    first = np.zeros(len(groups), bool)
    for g in np.unique(groups):
        idx = rng.permutation(np.where(groups == g)[0])
        first[idx[: int(round(share * len(idx)))]] = True
    return first


def load(ds):
    """The cells and the split of runs/loop_20261112 and runs/loop_guard_20261112 (same code, same seed); the audit controls are
    checked against the record of the first."""
    d = np.load(FEATURES / ds / "dense.npz"); D = {k: d[k] for k in d.files}
    Z = (D.pop("node__Z") if "node__Z" in D else np.load(FEATURES / ds / "node__Z.npy")).astype(F)
    D = {k: v for k, v in D.items() if not k.startswith("node__")}
    x = (Z / (np.linalg.norm(Z, axis=1, keepdims=True) + 1e-6)).astype(F)
    A, tg = D["A"], D["target"].astype(str)
    st = D["stratum"].astype(str) if "stratum" in D else np.full(len(A), "all")
    rng = np.random.default_rng(SEED)
    ctl = A == 0
    audit_ctl = ctl & halves(np.where(ctl, st, "treated"), rng)
    grp = np.where(A == 1, np.char.add(np.char.add(D["guide"].astype(str), "|"), st), np.char.add("control|", st))
    val = halves(grp, rng, 0.2) & ~audit_ctl
    ctl_group = np.full(len(A), -1); idx = rng.permutation(np.where(audit_ctl)[0]); ctl_group[idx] = np.arange(len(idx)) % 5
    assert sorted(D["ids"].astype(str)[audit_ctl].tolist()) == json.load(open(LOOP / "results" / f"audit_controls_{ds}.json"))["audit_controls"]
    return dict(ids=D["ids"].astype(str), x=x, A=A, tg=tg, y=D["readout__PDL1"].astype(np.float64), audit_ctl=audit_ctl, train_ctl=ctl & ~audit_ctl, val=val, ctl_group=ctl_group)


def weights():
    W = dict(np.load(FEATURES / "head_weights.npz")); names = [str(n) for n in W["adt_names"]]
    return {l + s: W[l + s].astype(F) for l in LAYERS for s in (".weight", ".bias")}, names, names.index("CD274")


def unit(z):
    z = z.astype(np.float64)
    return (z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)).astype(F)


def forward_np(P, x):
    """The head in numpy: latent (output of the MLP) and output (log of the predicted mean of each protein)."""
    from scipy.special import expit
    h, lat = x.astype(np.float64), None
    for name, act in zip(LAYERS, ACT):
        h = h @ P[name + ".weight"].T.astype(np.float64) + P[name + ".bias"]
        h = h * expit(h) if act else h
        if name == "mlp__net.6":
            lat = h
    return lat, np.log(np.maximum(np.logaddexp(0.0, h), 1e-6))


def bm_score(P, x, y, names):
    """The head's own evaluation: per protein, the correlation of measured counts with the predicted mean; mean over proteins."""
    pred = np.exp(forward_np(P, x)[1])
    r = {n: float(np.corrcoef(y[:, c], pred[:, c])[0, 1]) for c, n in enumerate(names) if y[:, c].std() > 1e-8}
    return {"mean_pearson": float(np.mean(list(r.values()))), "pdl1_pearson": r.get("CD274"), "n_proteins": len(r)}


def own_output_u(y, pred, A, cluster=None, seed=SEED):
    """u of the head's own output for one event, as the tool computes it for the `frozen` scope (scope._u_shift) and as
    runs/scope_20261116/scripts/frozen_arms.py read it for the heads made before: 1 - a * shift(pred) / shift(y), with a the slope
    of y on pred among the control cells. Returns the draws (row 0 is the point value) and the slope. Treated cells are resampled
    by `cluster` when it is given, otherwise by cell; control cells by cell. The same seed gives the same resamples for every arm."""
    from anm_diag import scope
    ti, ci_ = np.where(A == 1)[0], np.where(A == 0)[0]
    t, c = scope._Side(y[ti, None], None if cluster is None else cluster[ti]), scope._Side(y[ci_, None], None)
    rng = np.random.default_rng(seed + 1)
    Ct, Cc = scope._counts(t.E, B, rng), scope._counts(c.E, B, rng)
    u, _, a = scope._u_shift(t, c, pred[ti, None], pred[ci_, None], Ct, Cc, np.ones(1))
    return u, float(np.atleast_1d(a)[0])


def interval(draws):
    from anm_diag.stats import ci
    return [float(v) for v in ci(draws, ALPHA)]


def reading(iv):
    """Two-sided, at the tool's margin: supported if the whole interval is inside (-margin, margin); rejected if it is entirely outside."""
    lo, hi = iv
    if -MARGIN < lo and hi < MARGIN:
        return "supported"
    if lo > MARGIN or hi < -MARGIN:
        return "rejected"
    return "unresolved"
