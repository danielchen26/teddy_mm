#!/usr/bin/env python3
"""The stage with the first annotation position beside Z, and the SPEC files of the audit. No readout enters a stage.

  python build_zc.py --thp1 THP1_ENC.npz --mel STRATA.npz --annot-thp1 ANNOT_THP1 --annot-mel ANNOT_MEL --ckpt TEDDY_G_70M --work WORK [--results RESULTS]

THP1_ENC.npz: thp1_enc.npz of runs/encoder_stages_20261124 (Z, V, the units and the readout of the THP-1 stratum).
STRATA.npz: the same for the strata mel_co and mel_ifn (`<stratum>__Z`, `__V`, ...), with a sample of the joint stage of
runs/input_beside_z_20261123 for a check. ANNOT_*: what runs/annotation_input_20261124/scripts/annot.py wrote (mode_P.npz:
the output of layer 12 at the four annotation positions, with the placeholders in front of the genes).

Writes in WORK, for every stratum, `<stratum>_zc.npz`: the cells in the order of the earlier rounds, with
  Z      the state, as the earlier SPEC files give it to the tool
  V      the 512 components of the stratum's RNA
  ZV     Z at unit length beside V, by the recipe of runs/input_beside_z_20261123 (the same function)
  C      the output of layer 12 at the first annotation position (the disease placeholder): what the official class
         returns as `cell_emb` when the placeholders are given
  ZC     Z at unit length beside C, by the same recipe: each block divided by the square root of its total variance
         among the stratum's control cells
and one SPEC per target (DESIGN.md, Section 3), with ridge and boosting.
"""
import argparse
import collections
import copy
import json
import sys
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
RUNS = HERE.parents[1]
sys.path.insert(0, str(RUNS / "teg_training_20261119" / "scripts"))
sys.path.insert(0, str(RUNS / "teg_loop_20261120" / "scripts"))
sys.path.insert(0, str(RUNS / "input_beside_z_20261123" / "scripts"))
from build import DECISION  # noqa: E402
from build_zv import BOTH, TARGETS, TYPE_A, joint  # noqa: E402

ORDER = ["disease", "tissue", "cell_type", "sex"]
STAGES = {t: (["V"] if "V" in st else []) + ["Z", "ZV", "C", "ZC"] for t, (_g, _s, _r, _e, _p, st) in TARGETS.items()}


def spec(name, study, readout, events, pooled, arrays, stages):
    return dict(name=name, kind="perturbation", study=study, readout_name=readout, events=list(events), also_pooled=bool(pooled), units_exclude=["control"], reads="ZC",
                control_folds=5, decision=DECISION, label_basis="design", keep_draws=True, readout_classes=copy.deepcopy(BOTH),
                keys=dict(readouts=[f"y_{readout}"], stages={k: k for k in stages}, unit="unit", control={"key": "A", "equals": 0}, cluster="unit",
                          frozen={"array": f"frozen_{readout}"}), arrays={"audit": arrays})


def said(top, lab, words, lo):
    out = {}
    for j, name in enumerate(ORDER):
        c = collections.Counter(words[int(i)] if int(i) < len(words) else f"id {int(i)}" for i in top[:, j])
        out[name] = {"ranked_first_share_of_cells": {k: round(v / len(top), 4) for k, v in c.most_common(4)}, "mean_probability_of_all_label_tokens": round(float(lab[:, j].sum(1).mean()), 4)}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--thp1", type=Path, required=True)
    ap.add_argument("--mel", type=Path, required=True)
    ap.add_argument("--annot-thp1", type=Path, required=True)
    ap.add_argument("--annot-mel", type=Path, required=True)
    ap.add_argument("--ckpt", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=HERE.parent / "results")
    a = ap.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    words = [line.rstrip("\n") for line in open(a.ckpt / "vocab.txt")]
    T, M = np.load(a.thp1), np.load(a.mel)
    PT, PM = np.load(a.annot_thp1 / "mode_P.npz"), np.load(a.annot_mel / "mode_P.npz")
    mt, mm = json.loads((a.annot_thp1 / "manifest.json").read_text()), json.loads((a.annot_mel / "manifest.json").read_text())
    assert mt["check_passed"] and mm["check_passed"], "a state without annotations does not repeat the record"
    lo = mm["label_tokens"][0]
    pos = {i: k for k, i in enumerate(PM["ids"].astype(str))}
    info = {"strata": {}, "state_without_annotations_against_the_recorded_state": {"thp1": mt["N_against_the_recorded_state"], "mel_co and mel_ifn": mm["N_against_the_recorded_state"]}}
    for g in ("thp1", "mel_co", "mel_ifn"):
        if g == "thp1":
            K = {k: T[k] for k in ("Z", "V", "ids", "A", "unit", "y_PDL1", "frozen_PDL1")}
            assert np.array_equal(PT["ids"].astype(str), K["ids"].astype(str)), "the annotation positions are of other cells"
            C, top, lab = PT["A"][:, 0], PT["top"], PT["lab"]
        else:
            K = {k[len(g) + 2:]: M[k] for k in M.files if k.startswith(g + "__") and not k.endswith("ZV_sample")}
            kk = np.array([pos[i] for i in K["ids"].astype(str)])
            C, top, lab = PM["A"][kk, 0], PM["top"][kk], PM["lab"][kk]
        ctl = np.where(K["A"] == 0)[0]
        Z, V, C = K["Z"].astype(np.float32), K["V"].astype(np.float32), C.astype(np.float32)
        ZV, _a, sv = joint(Z, V, ctl)
        ZC, sa, sc = joint(Z, C, ctl)
        d = float(np.abs(ZV[::97, ::31] - M[f"{g}__ZV_sample"]).max())
        assert d <= 1e-5, f"the joint stage of the run before is not rebuilt ({g}: {d})"
        rest = {k: v for k, v in K.items() if k not in ("Z", "V")}
        assert not (a.work / f"{g}_zc.npz").is_symlink()
        np.savez(a.work / f"{g}_zc.npz", Z=Z, V=V, ZV=ZV, C=C, ZC=ZC, **rest)
        info["strata"][g] = {"cells": int(len(Z)), "control_cells": int(len(ctl)), "scale_of_Z": sa, "scale_of_V": sv, "scale_of_C": sc, "joint_stage_of_the_run_before_rebuilt_to": d,
                             "the_model_says_at_the_placeholders": said(top, lab, words, lo), "the_same_among_control_cells": said(top[ctl], lab[ctl], words, lo)}
        print(g, ZC.shape, "| scale of Z", round(sa, 4), "| scale of C", round(sc, 3), "| ZV rebuilt to", d)
        if g == "thp1":
            m = np.isin(K["unit"].astype(str), TYPE_A + ["control"])
            np.savez(a.work / "thp1a_zc.npz", Z=Z[m], V=V[m], ZV=ZV[m], C=C[m], ZC=ZC[m], **{k: v[m] for k, v in rest.items()})
    for t, (g, study, readout, events, pooled, _st) in TARGETS.items():
        s = spec(f"{t}_{g}_{readout}", study, readout, events, pooled, f"{g}_zc.npz", STAGES[t])
        (a.work / f"{s['name']}.yaml").write_text(yaml.safe_dump(s, sort_keys=False))
        print("SPEC", s["name"], "| stages", STAGES[t], "| events", events, "| pooled", pooled)
    a.results.mkdir(parents=True, exist_ok=True)
    (a.results / "stages_manifest.json").write_text(json.dumps({**info, "stages": STAGES, "states": {"thp1": {k: mt[k] for k in ("device", "torch", "python", "seconds")}, "mel": {k: mm[k] for k in ("device", "torch", "python", "seconds")}}}, indent=1) + "\n")


if __name__ == "__main__":
    main()
