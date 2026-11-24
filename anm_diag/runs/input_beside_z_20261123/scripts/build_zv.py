#!/usr/bin/env python3
"""The joint stage and the SPEC files of the audit. No readout and no label of an event enters the stage.

  python build_zv.py --work WORK

Environment as runs/teg_loop_20261120 (ANM_LOOP_WORK: the work folder of its build.py, with <stratum>.npz and <stratum>_rna.npz).
Writes in WORK, for every stratum, `<stratum>_zv.npz`: the cells of the stratum in the order of the earlier rounds, with
  Z      the state as the earlier SPEC files give it to the tool (the round puts it at unit length, as the head does)
  V      the 512 components of the stratum's RNA
  ZV     the state at unit length and V side by side, each divided by the square root of its total variance among the
         stratum's control cells
and one SPEC per target (DESIGN.md, Section 3) that reads Z and ZV in one run, with ridge and boosting.
"""
import argparse
import copy
import os
import sys
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
RUNS = HERE.parents[1]
sys.path.insert(0, str(RUNS / "teg_training_20261119" / "scripts"))
sys.path.insert(0, str(RUNS / "teg_loop_20261120" / "scripts"))
from _common import unit  # noqa: E402
from build import ALPHAS, DECISION, HGB  # noqa: E402

BOTH = {"ridge": {"regressor": "RidgeCV", "alphas": ALPHAS}, "hgb": HGB}
TYPE_A = ["IFNGR1", "IFNGR2", "JAK2", "STAT1"]
PATHWAY = ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"]
# target -> (stratum file, study, readout, events, pooled, stages)
TARGETS = {"T1": ("thp1", "papalexi2021", "PDL1", ["IRF1"], False, ["V", "Z", "ZV"]),
           "T2": ("mel_co", "frangieh2021", "MHC1", PATHWAY, True, ["V", "Z", "ZV"]),
           "C1": ("thp1a", "papalexi2021", "PDL1", TYPE_A, True, ["Z", "ZV"]),
           "C2": ("mel_co", "frangieh2021", "PDL1", PATHWAY, True, ["Z", "ZV"]),
           "C3": ("mel_ifn", "frangieh2021", "PDL1", PATHWAY, True, ["Z", "ZV"])}


def joint(z, v, ctl):
    """Z at unit length beside V, each block scaled to the same total variance among control cells."""
    z = unit(z).astype(np.float64)
    v = v.astype(np.float64)
    a, b = np.sqrt(z[ctl].var(0).sum()), np.sqrt(v[ctl].var(0).sum())
    return np.hstack([z / a, v / b]).astype(np.float32), float(a), float(b)


def spec(name, study, readout, events, pooled, arrays, stages):
    return dict(name=name, kind="perturbation", study=study, readout_name=readout, events=list(events), also_pooled=bool(pooled), units_exclude=["control"], reads="ZV",
                control_folds=5, decision=DECISION, label_basis="design", keep_draws=True, readout_classes=copy.deepcopy(BOTH),
                keys=dict(readouts=[f"y_{readout}"], stages={k: k for k in stages}, unit="unit", control={"key": "A", "equals": 0}, cluster="unit",
                          frozen={"array": f"frozen_{readout}"}), arrays={"audit": arrays})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    a = ap.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    src = Path(os.environ["ANM_LOOP_WORK"])
    for g in ("thp1", "mel_co", "mel_ifn"):
        K, R = dict(np.load(src / f"{g}.npz")), np.load(src / f"{g}_rna.npz")
        assert np.array_equal(K["ids"].astype(str), R["ids"].astype(str)), "the state and the RNA view are not of the same cells"
        ZV, sa, sb = joint(K["Z"], R["RNA512"], np.where(K["A"] == 0)[0])
        Z, V, rest = K["Z"].astype(np.float32), R["RNA512"].astype(np.float32), {k: v for k, v in K.items() if k != "Z"}
        assert not (a.work / f"{g}_zv.npz").is_symlink()
        np.savez(a.work / f"{g}_zv.npz", Z=Z, V=V, ZV=ZV, **rest)
        print(g, ZV.shape, "| scale of Z", round(sa, 4), "| scale of V", round(sb, 2))
        if g == "thp1":                                         # the four genes whose pool does not contain IRF1, as the first outer step
            m = np.isin(K["unit"].astype(str), TYPE_A + ["control"])
            np.savez(a.work / "thp1a_zv.npz", Z=Z[m], V=V[m], ZV=ZV[m], **{k: v[m] for k, v in rest.items()})
    for t, (g, study, readout, events, pooled, stages) in TARGETS.items():
        s = spec(f"{t}_{g}_{readout}", study, readout, events, pooled, f"{g}_zv.npz", stages)
        (a.work / f"{s['name']}.yaml").write_text(yaml.safe_dump(s, sort_keys=False))
        print("SPEC", s["name"], "| stages", stages, "| events", events, "| pooled", pooled)


if __name__ == "__main__":
    main()
