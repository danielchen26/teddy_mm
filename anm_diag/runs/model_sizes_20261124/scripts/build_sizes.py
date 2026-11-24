#!/usr/bin/env python3
"""The arrays and SPEC files for the three sizes of TEDDY-G, for `anm_diag scope`. No readout enters a stage.

  python build_sizes.py --enc WORK_ENC --annot70 ANNOT_70M --s160 STATES_160M --s400 STATES_400M --ckpt TEDDY_G_70M --work WORK [--results RESULTS]

WORK_ENC: the work folder of runs/encoder_stages_20261124 (thp1_enc.npz and its two SPEC files: the input view, the
readout, the units, the permuted view, and the stages of the 70M checkpoint). ANNOT_70M: what
runs/annotation_input_20261124/scripts/annot.py wrote (mode_P.npz). STATES_*: what scripts/sizes.py wrote.

Stages, all of the same 7,478 cells (DESIGN.md, Section 3). For a size m in 70, 160, 400:

  E<m>   what enters the first layer, averaged over the cell's gene tokens
  P<m>   the output of the layer before the last, averaged
  Z<m>   the output of the last layer, averaged: the state a head would read
  C<m>   the output of the last layer at the first annotation position, with the placeholders in front of the genes

and for the two larger sizes, whose states are wider than 512:

  Z<m>r, C<m>r   the same reduced to 512 components by the recipe of the input view (TruncatedSVD, 4 iterations, the
                 same seed), fitted on all cells: a state of the width of the small model's
  ZP<m>          the state with the placeholders in front of the genes

V is the input view and PERM the same with its rows permuted (a control: it must read nothing).
"""
import argparse
import collections
import json
from pathlib import Path

import numpy as np
import yaml
from sklearn.decomposition import TruncatedSVD

SEED = 20261120                                                  # the seed of the input view (teg_loop_20261120/scripts/build_rna.py)
TYPE_A = ["IFNGR1", "IFNGR2", "JAK2", "STAT1"]
ORDER = ["disease", "tissue", "cell_type", "sex"]
STAGES = ["V", "E70", "P70", "Z70", "C70"] + [f"{k}{m}" for m in (160, 400) for k in ("E", "P", "Z", "C")] + [f"{k}{m}r" for m in (160, 400) for k in ("Z", "C")] + ["ZP160", "ZP400", "PERM"]
comp = lambda X: TruncatedSVD(n_components=512, algorithm="randomized", n_iter=4, random_state=SEED).fit_transform(X).astype(np.float32)


def said(top, lab, words):
    out = {}
    for j, name in enumerate(ORDER):
        c = collections.Counter(words[int(i)] if int(i) < len(words) else f"id {int(i)}" for i in top[:, j])
        out[name] = {"ranked_first_share_of_cells": {k: round(v / len(top), 4) for k, v in c.most_common(4)}, "mean_probability_of_all_label_tokens": round(float(lab[:, j].sum(1).mean()), 4)}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--enc", type=Path, required=True)
    ap.add_argument("--annot70", type=Path, required=True)
    ap.add_argument("--s160", type=Path, required=True)
    ap.add_argument("--s400", type=Path, required=True)
    ap.add_argument("--ckpt", type=Path, required=True, help="any of the three checkpoints: for the words of the vocabulary")
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=Path(__file__).resolve().parent.parent / "results")
    a = ap.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    words = [line.rstrip("\n") for line in open(a.ckpt / "vocab.txt")]
    K = np.load(a.enc / "thp1_enc.npz")
    ids = K["ids"].astype(str)
    P70 = np.load(a.annot70 / "mode_P.npz")
    assert np.array_equal(P70["ids"].astype(str), ids)
    S = {"V": K["V"].astype(np.float32), "E70": K["E0"].astype(np.float32), "P70": K["L11"].astype(np.float32), "Z70": K["L12"].astype(np.float32), "C70": P70["A"][:, 0].astype(np.float32)}
    info = {"cells": int(len(ids)), "sizes": {"70": {"layers": 12, "d_model": 512, "the_model_says_at_the_placeholders": said(P70["top"], P70["lab"], words)}}}
    big = {}
    for m, d in ((160, a.s160), (400, a.s400)):
        X = np.load(d / "states.npz")
        assert np.array_equal(X["ids"].astype(str), ids), f"the states of the {m}M checkpoint are of other cells"
        man = json.loads((d / "manifest.json").read_text())
        chk = json.loads((d / "check.json").read_text()) if (d / "check.json").exists() else None
        for k, key in (("E", "E0"), ("P", "LP"), ("Z", "Z"), ("C", "C")):
            S[f"{k}{m}"] = X[key].astype(np.float32)
        big[m] = X
        info["sizes"][str(m)] = {**{k: man[k] for k in ("layers", "d_model", "heads", "d_hid", "weights", "device", "torch", "seconds")},
                                 "the_model_says_at_the_placeholders": said(X["top"], X["lab"], words),
                                 "code_check": None if chk is None else {k: chk[k] for k in ("cells", "hold", "official_class_against_this_code_on_cpu", "this_device_against_cpu", "order_of_the_annotations")}}
    for m in (160, 400):
        S[f"Z{m}r"], S[f"C{m}r"] = comp(big[m]["Z"]), comp(big[m]["C"])
    for m in (160, 400):
        S[f"ZP{m}"] = big[m]["ZP"].astype(np.float32)
    S["PERM"] = K["PERM"].astype(np.float32)
    assert list(S) == STAGES, list(S)
    info["width"] = {k: int(v.shape[1]) for k, v in S.items()}
    rest = {k: K[k] for k in ("ids", "A", "unit", "y_PDL1", "frozen_PDL1")}
    for f in ("thp1_sizes.npz", "thp1a_sizes.npz"):
        assert not (a.work / f).is_symlink()
    np.savez(a.work / "thp1_sizes.npz", **S, **rest)
    msk = np.isin(K["unit"].astype(str), TYPE_A + ["control"])
    np.savez(a.work / "thp1a_sizes.npz", **{k: v[msk] for k, v in S.items()}, **{k: v[msk] for k, v in rest.items()})
    for src, arrays in (("T_thp1_IRF1", "thp1_sizes.npz"), ("C_thp1a_pool", "thp1a_sizes.npz")):
        s = yaml.safe_load((a.enc / f"{src}.yaml").read_text())  # the same events, readers, decision and keys as the runs before; only the stages change
        s["keys"]["stages"] = {k: k for k in STAGES}
        s["reads"] = "Z70"
        s["arrays"] = {"audit": arrays}
        (a.work / f"{src}.yaml").write_text(yaml.safe_dump(s, sort_keys=False))
    a.results.mkdir(parents=True, exist_ok=True)
    (a.results / "sizes_manifest.json").write_text(json.dumps({**info, "stages": STAGES}, indent=1) + "\n")
    print(json.dumps({"cells": info["cells"], "width": info["width"]}, indent=1))


if __name__ == "__main__":
    main()
