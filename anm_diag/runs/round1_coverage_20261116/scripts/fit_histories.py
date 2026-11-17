#!/usr/bin/env python3
"""The validation curves of the six heads fitted for the test, copied from their metrics.json (the heads are not in
the repository). Written after the test, for the record of the head that is not admitted.

  python fit_histories.py --heads HEADS_DIR

Writes ../results/gate5_fit_histories.json.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import common
import plan


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--heads", type=Path, required=True)
    a = p.parse_args()
    out = {"what": "validation score on the bone-marrow validation donor after every epoch (the mean correlation over proteins), for the heads fitted at gate 5",
           "heads": {}}
    for s in plan.SEEDS:
        for arm in ("coverage_all", "labels_all"):
            n = f"s{s}_{arm}"
            m = json.loads((a.heads / n / "metrics.json").read_text())
            out["heads"][n] = {"checkpoint_epoch": m["checkpoint_epoch"], "validation_score": [round(float(h["val_mlp_pearson"]), 4) for h in m["history"]],
                               "train_loss": [round(float(h["train_loss"]), 4) for h in m["history"]]}
    (common.RUN / "results" / "gate5_fit_histories.json").write_text(json.dumps(out, indent=1) + "\n")
    for n, t in out["heads"].items():
        print(n, t["checkpoint_epoch"], t["validation_score"][-6:])


if __name__ == "__main__":
    main()
