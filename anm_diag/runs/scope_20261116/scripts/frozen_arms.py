#!/usr/bin/env python3
"""An added check, not in DESIGN.md: the head's own PD-L1 output, before and after the repairs already made.

The design's forecast for G3 was written for the quantity the earlier audits measured (|u| at HEAD_OUT in the units
scope). The tool's own forecast (scope.forecast) is about the system's own output, the `frozen` scope. The outputs of
the repaired heads on their audit cells were kept, so that quantity can be measured for every arm without a new fit:

  python frozen_arms.py --loop-features DIR --sealed-melanoma FEATURES_DIR --sealed-thp1 FEATURES_DIR

DIR holds pooled_<study>_<arm>/ and transfer_<study>_<arm>/ (dense.npz, node__HEAD_OUT.npy) as written by
runs/loop_20261112. Writes ../results/frozen_arms.json.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RUN.parents[1]))
from anm_diag import scope  # noqa: E402
from anm_diag.stats import ci  # noqa: E402

J_PDL1, B, SEED, ALPHA = 1, 1000, 0, 0.05


def load(d: Path):
    z = np.load(d / "dense.npz", allow_pickle=False)
    out = z["node__HEAD_OUT"] if "node__HEAD_OUT" in z.files else np.load(d / "node__HEAD_OUT.npy")
    return dict(ids=z["ids"].astype(str), A=z["A"].astype(int), target=z["target"].astype(str), y=z["readout__PDL1"].astype(np.float64),
                pred=np.exp(out[:, J_PDL1].astype(np.float64)))


def frozen_u(x):
    ti, ci_ = np.where(x["A"] == 1)[0], np.where(x["A"] == 0)[0]
    t, c = scope._Side(x["y"][ti, None], x["target"][ti]), scope._Side(x["y"][ci_, None], None)
    rng = np.random.default_rng(SEED + 1)
    Ct, Cc = scope._counts(t.E, B, rng), scope._counts(c.E, B, rng)
    u, _, a = scope._u_shift(t, c, x["pred"][ti, None], x["pred"][ci_, None], Ct, Cc, np.ones(1))
    return u, float(np.atleast_1d(a)[0])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--loop-features", type=Path, required=True)
    p.add_argument("--sealed-melanoma", type=Path, required=True)
    p.add_argument("--sealed-thp1", type=Path, required=True)
    a = p.parse_args()
    res = {}
    for ct, ds, sealed in (("melanoma", "frangieh2021", a.sealed_melanoma), ("thp1", "papalexi2021", a.sealed_thp1)):
        sets = {"cells of the loop (half of the controls)": ("pooled", a.loop_features / f"pooled_{ds}_as_trained",
                                                              {"coverage": f"pooled_{ds}_mlp_only", "elsewhere": f"pooled_{ds}_decoder_only"}),
                "all cells of the sealed audit": ("transfer", sealed, {"other study": f"transfer_{ds}_mlp_only",
                                                                       "other study, elsewhere": f"transfer_{ds}_decoder_only"})}
        res[ct] = {}
        for label, (_, base_dir, arms) in sets.items():
            base = load(base_dir)
            ub, ab = frozen_u(base)
            entry = {"n_cells": int(len(base["A"])), "as_trained": {"u": float(ub[0]), "ci": ci(ub[1:], ALPHA), "slope": ab}, "arms": {}}
            for arm, folder in arms.items():
                d = a.loop_features / folder
                if not d.exists():
                    continue
                x = load(d)
                assert np.array_equal(x["ids"], base["ids"]) and np.array_equal(x["A"], base["A"]), f"{folder}: other cells than its baseline"
                u, sl = frozen_u(x)
                fall = np.abs(ub) - np.abs(u)
                entry["arms"][arm] = {"u": float(u[0]), "ci": ci(u[1:], ALPHA), "slope": sl, "fall": float(fall[0]), "fall_ci": ci(fall[1:], ALPHA),
                                      "run": f"loop_20261112/{folder}"}
            res[ct][label] = entry
            print(ct, "|", label, "| as trained", round(entry["as_trained"]["u"], 3),
                  {k: (round(v["u"], 3), "fall", round(v["fall"], 3), [round(q, 3) for q in v["fall_ci"]]) for k, v in entry["arms"].items()})
    (RUN / "results" / "frozen_arms.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
