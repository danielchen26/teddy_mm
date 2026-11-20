#!/usr/bin/env python3
"""Arrays and SPEC files for `anm_diag scope`: three groups of events, the head as trained. No fit is made here.

  python build.py --work WORK

Environment: ANM_FEATURES, ANM_BM_FIT (as runs/teg_training_20261119) and ANM_TEST_CELLS (the folder written by
runs/teg_repair_test_20261119/scripts/cells.py: coculture.npz, library.npz).

A group is a cell type in a culture condition: the type every event of the group shares.
  mel_ifn   melanoma, IFN-gamma condition     PD-L1, MHC-I, HLA-E
  mel_co    melanoma, co-culture condition    PD-L1, MHC-I, HLA-E
  thp1      THP-1                             PD-L1 (the head has no PD-L2 output)
and one probe: the knockout of CD274, the gene of PD-L1 itself, in the IFN-gamma condition.

For every group and readout one SPEC: the stage is Z (the embedding the head reads); `frozen` is the head's own output on
the log scale; `train` is the head's bone-marrow training cells; `other` is another group (the other condition of the
same cell type, or for THP-1 the melanoma cells).
"""
import argparse
import copy
import os
import sys
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "teg_training_20261119" / "scripts"))
from _common import BM_FIT, FEATURES, TARGETS, forward_np, unit, weights  # noqa: E402

ALPHAS = [float(a) for a in np.logspace(-1, 5, 13)]
HGB = {"regressor": "HistGradientBoostingRegressor", "max_iter": 150, "learning_rate": 0.05, "max_leaf_nodes": 15}
DECISION = {"margin_u": 0.10, "step_margin": 0.10, "alpha": 0.05, "bootstrap": 1000, "seed": 0, "min_clusters": 10}
OUTPUT = {"PDL1": "CD274", "MHC1": "HLA-A-B-C", "HLAE": "HLA-E"}       # readout -> output of the head
GROUPS = {"mel_ifn": dict(cell_type="melanoma", stratum="IFN-gamma", study="frangieh2021", readouts=["PDL1", "MHC1", "HLAE"], other="mel_co"),
          "mel_co": dict(cell_type="melanoma", stratum="co-culture", study="frangieh2021", readouts=["PDL1", "MHC1", "HLAE"], other="mel_ifn"),
          "thp1": dict(cell_type="THP-1", stratum="stimulated", study="papalexi2021", readouts=["PDL1"], other="mel_ifn")}
PROBE = "CD274"


def group_arrays(name):
    cells = Path(os.environ["ANM_TEST_CELLS"])
    if name == "mel_co":
        d = np.load(cells / "coculture.npz")
        Z = d["Z"]
    else:
        ds = GROUPS[name]["study"]
        d = np.load(FEATURES / ds / "dense.npz")
        Z = d["node__Z"] if "node__Z" in d.files else np.load(FEATURES / ds / "node__Z.npy")
    A = d["A"].astype(int)
    tg = d["target"].astype(str)
    assert ((A == 0) | np.isin(tg, TARGETS[GROUPS[name]["study"]])).all(), name
    return dict(ids=d["ids"].astype(str), Z=Z.astype(np.float32), A=A, unit=np.where(A == 1, tg, "control"),
                **{f"y_{r}": d[f"readout__{r}"].astype(np.float64) for r in GROUPS[name]["readouts"]})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    a = ap.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    P0, names, _j = weights()
    col = {r: names.index(n) for r, n in OUTPUT.items()}
    head = lambda Z: forward_np(P0, unit(Z))[1]
    for name, g in GROUPS.items():
        K = group_arrays(name)
        out = head(K["Z"])
        for r in g["readouts"]:
            K[f"frozen_{r}"] = out[:, col[r]].astype(np.float64)
        np.savez(a.work / f"{name}.npz", **K)
        u, c = np.unique(K["unit"], return_counts=True)
        print(name, dict(zip(u.tolist(), c.tolist())))
    b = np.load(BM_FIT)
    assert [str(x) for x in b["adt_names"]] == names
    np.savez(a.work / "bm_train.npz", Z=b["z_train"].astype(np.float32), **{f"y_{r}": b["adt_train"][:, col[r]].astype(np.float64) for r in OUTPUT})
    # the probe: the IFN-gamma cells with the CD274 knockout cells of the library added as one more unit
    lib = np.load(Path(os.environ["ANM_TEST_CELLS"]) / "library.npz")
    K = dict(np.load(a.work / "mel_ifn.npz"))
    m = lib["target"].astype(str) == PROBE
    assert not np.isin(lib["ids"][m].astype(str), K["ids"]).any()
    P = dict(ids=np.r_[K["ids"], lib["ids"][m].astype(str)], Z=np.vstack([K["Z"], lib["Z"][m]]), A=np.r_[K["A"], np.ones(m.sum(), int)],
             unit=np.r_[K["unit"], np.full(m.sum(), PROBE)], y_PDL1=np.r_[K["y_PDL1"], lib["readout__PDL1"][m]])
    P["frozen_PDL1"] = head(P["Z"])[:, col["PDL1"]].astype(np.float64)
    np.savez(a.work / "mel_ifn_probe.npz", **P)
    print("probe", PROBE, int(m.sum()), "cells")
    base = dict(kind="perturbation", units_exclude=["control"], reads="Z", control_folds=5, decision=DECISION, label_basis="design", keep_draws=True,
                also_pooled=True, readout_classes={"ridge": {"regressor": "RidgeCV", "alphas": ALPHAS}, "hgb": HGB})
    n = 0
    for name, g in GROUPS.items():
        for r in g["readouts"]:
            s = {**copy.deepcopy(base), "name": f"{name}_{r}", "study": g["study"], "readout_name": r, "events": list(TARGETS[g["study"]]),
                 "keys": dict(readouts=[f"y_{r}"], stages={"Z": "Z"}, unit="unit", control={"key": "A", "equals": 0}, cluster="unit", frozen={"array": f"frozen_{r}"}),
                 "arrays": {"audit": f"{name}.npz", "train": "bm_train.npz"}}
            if r in GROUPS[g["other"]]["readouts"]:
                s["arrays"]["other"] = f"{g['other']}.npz"
            (a.work / f"{name}_{r}.yaml").write_text(yaml.safe_dump(s, sort_keys=False))
            n += 1
    s = {**copy.deepcopy(base), "name": "probe_PDL1", "study": "frangieh2021", "readout_name": "PDL1", "events": [PROBE], "also_pooled": False,
         "keys": dict(readouts=["y_PDL1"], stages={"Z": "Z"}, unit="unit", control={"key": "A", "equals": 0}, cluster="unit", frozen={"array": "frozen_PDL1"}),
         "arrays": {"audit": "mel_ifn_probe.npz", "train": "bm_train.npz"}}
    (a.work / "probe_PDL1.yaml").write_text(yaml.safe_dump(s, sort_keys=False))
    print("wrote", n + 1, "SPEC files")


if __name__ == "__main__":
    main()
