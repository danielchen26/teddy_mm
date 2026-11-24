#!/usr/bin/env python3
"""The arrays and SPEC files of the three ways of running the encoder, for `anm_diag scope`. No readout enters a stage.

  python build_modes.py --enc WORK_ENC --annot ANNOT_OUT --ckpt TEDDY_G_70M --work WORK [--results RESULTS]

WORK_ENC: the work folder of runs/encoder_stages_20261124 (thp1_enc.npz and its two SPEC files; the input view, the
readout, the units and the permuted view are taken from it unchanged). ANNOT_OUT: what scripts/annot.py wrote.

Stages, all of the same 7,478 cells (DESIGN.md, Section 2):

  V     the input view of the earlier rounds, as it is
  ZN    the state without annotation positions (the record)
  ZP    the state with the four placeholders in front of the genes
  ZT    the state with the labels of the THP-1 line in front of the genes
  ZD    the state with the labels the official tokeniser writes when a study has no annotation columns
  ZW    the state with labels that are wrong for these cells
  AP    the output of the last layer at the four annotation positions, with the placeholders: 2,048 numbers reduced to
        512 components by the recipe of the input view (TruncatedSVD, 4 iterations, the same seed), fitted on all cells
  AT    the same with the labels of the line
  CP    the first position with the placeholders, as it is (what the official class returns as `cell_emb`)
  ZAP   ZP beside AP (1,024 numbers)
  PERM  V with its rows permuted (a control: it must read nothing)

Also written, to RESULTS/modes_manifest.json: how far the state moves between the ways of running, what the model's own
output layer says at the annotation positions, and how large the input embeddings of the label tokens are.
"""
import argparse
import collections
import json
from pathlib import Path

import numpy as np
import yaml
from sklearn.decomposition import TruncatedSVD

SEED = 20261120                                                  # the seed of the input view (teg_loop_20261120/scripts/build_rna.py)
Z_STAGES = ["ZN", "ZP", "ZT", "ZD", "ZW"]
STAGES = ["V"] + Z_STAGES + ["AP", "AT", "CP", "ZAP", "PERM"]
TYPE_A = ["IFNGR1", "IFNGR2", "JAK2", "STAT1"]
ORDER = ["disease", "tissue", "cell_type", "sex"]
comp = lambda X: TruncatedSVD(n_components=512, algorithm="randomized", n_iter=4, random_state=SEED).fit_transform(X).astype(np.float32)


def moved(Zm, Zn):
    """How far the state of a cell moves when the way of running changes, as a share of how far the cell lies from the
    mean of the stratum: as it is, and after the shift common to all cells is taken out (a reader with an intercept does
    not see that shift)."""
    cn, cm = Zn - Zn.mean(0), Zm - Zm.mean(0)
    far = np.linalg.norm(cn, axis=1)
    raw, cen = np.linalg.norm(Zm - Zn, axis=1) / far, np.linalg.norm(cm - cn, axis=1) / far
    q = lambda x: {"median": float(np.median(x)), "p10": float(np.quantile(x, 0.1)), "p90": float(np.quantile(x, 0.9))}
    return {"as_it_is": q(raw), "common_shift_taken_out": q(cen), "common_shift_over_median_distance_from_the_mean": float(np.linalg.norm(Zm.mean(0) - Zn.mean(0)) / np.median(far))}


def said(M, words, lo, truth):
    """What the model's output layer ranks first at each annotation position, over all cells, and what it gives the labels."""
    out = {}
    for j, name in enumerate(ORDER):
        top = collections.Counter(words[int(i)] if int(i) < len(words) else f"id {int(i)}" for i in M["top"][:, j])
        mean = M["lab"][:, j].mean(0)
        best = np.argsort(-mean)[:4]
        out[name] = {"ranked_first_share_of_cells": {k: round(v / len(M["top"]), 4) for k, v in top.most_common(5)},
                     "mean_probability": {words[lo + int(b)]: round(float(mean[b]), 4) for b in best},
                     "mean_probability_of_all_label_tokens": round(float(M["lab"][:, j].sum(1).mean()), 4),
                     "label_of_the_line": truth[j], "mean_probability_of_the_label_of_the_line": round(float(mean[words.index(truth[j]) - lo]), 4)}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--enc", type=Path, required=True)
    ap.add_argument("--annot", type=Path, required=True)
    ap.add_argument("--ckpt", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=Path(__file__).resolve().parent.parent / "results")
    a = ap.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    K = np.load(a.enc / "thp1_enc.npz")
    ids = K["ids"].astype(str)
    man = json.loads((a.annot / "manifest.json").read_text())
    assert man["check_passed"], "the state without annotations does not repeat the record"
    M = {k: np.load(a.annot / f"mode_{k}.npz") for k in "NPTDW"}
    for q in M.values():
        assert np.array_equal(q["ids"].astype(str), ids), "a file of the modes is for other cells"
    assert float(np.abs(M["N"]["Z"] - K["L12"]).max()) <= 1e-4
    n = len(ids)
    AP, AT = comp(M["P"]["A"].reshape(n, -1)), comp(M["T"]["A"].reshape(n, -1))
    S = {"V": K["V"].astype(np.float32), **{f"Z{k}": M[k]["Z"].astype(np.float32) for k in "NPTDW"}, "AP": AP, "AT": AT, "CP": M["P"]["A"][:, 0].astype(np.float32),
         "ZAP": np.concatenate([M["P"]["Z"].astype(np.float32), AP], axis=1), "PERM": K["PERM"].astype(np.float32)}
    assert list(S) == STAGES
    rest = {k: K[k] for k in ("ids", "A", "unit", "y_PDL1", "frozen_PDL1")}
    for f in ("thp1_ann.npz", "thp1a_ann.npz"):
        assert not (a.work / f).is_symlink()
    np.savez(a.work / "thp1_ann.npz", **S, **rest)
    m = np.isin(K["unit"].astype(str), TYPE_A + ["control"])
    np.savez(a.work / "thp1a_ann.npz", **{k: v[m] for k, v in S.items()}, **{k: v[m] for k, v in rest.items()})
    for src, arrays in (("T_thp1_IRF1", "thp1_ann.npz"), ("C_thp1a_pool", "thp1a_ann.npz")):
        s = yaml.safe_load((a.enc / f"{src}.yaml").read_text())  # the same events, readers, decision and keys as the run before; only the stages change
        s["keys"]["stages"] = {k: k for k in STAGES}
        s["reads"] = "ZN"
        s["arrays"] = {"audit": arrays}
        (a.work / f"{src}.yaml").write_text(yaml.safe_dump(s, sort_keys=False))
    # --- described, not judged
    words = [line.rstrip("\n") for line in open(a.ckpt / "vocab.txt")]
    lo, hi = man["label_tokens"][0], man["label_tokens"][1] + 1
    info = {"cells": int(n), "cells_of_the_control_file": int(m.sum()), "modes": man["modes"], "token_ids": man["token_ids"], "N_against_the_recorded_state": man["N_against_the_recorded_state"],
            "the_state_moves": {f"Z{k}": moved(M[k]["Z"], M["N"]["Z"]) for k in "PTDW"},
            "the_state_moves_between_the_labels_of_the_line_and_wrong_labels": moved(M["W"]["Z"], M["T"]["Z"]),
            "the_model_says_at_the_placeholders": said(M["P"], words, lo, [man["modes"]["T"][o] for o in ORDER]),
            "states": {k: man[k] for k in ("device", "torch", "python", "seconds")}}
    from safetensors import safe_open
    with safe_open(str(a.ckpt / "model.safetensors"), "np") as f:
        E, Dh = f.get_tensor("embeddings.weight"), f.get_tensor("decoder_head.weight")
    norm = np.linalg.norm(E, axis=1)
    special = [words.index(w) for w in ("<disease>", "<tissue>", "<cell>", "<sex>")]
    unused = np.arange(len(words) + 5, E.shape[0])              # ids past the five special tokens of the tokeniser: no token has them
    med = lambda x: float(np.median(x))
    info["input_embeddings"] = {"median_norm": {"label_tokens": med(norm[lo:hi]), "placeholders": med(norm[special]), "ids_no_token_has": med(norm[unused]), "all_gene_tokens": med(norm[:lo])},
                                "label_tokens_min_max": [float(norm[lo:hi].min()), float(norm[lo:hi].max())],
                                "input_table_equals_output_table": bool(np.array_equal(E, Dh)), "largest_difference_input_output_table": float(np.abs(E - Dh).max())}
    info["input_embeddings"]["label_tokens_look_untrained"] = bool(info["input_embeddings"]["median_norm"]["label_tokens"] < 2.0 * info["input_embeddings"]["median_norm"]["ids_no_token_has"])
    a.results.mkdir(parents=True, exist_ok=True)
    (a.results / "modes_manifest.json").write_text(json.dumps({**info, "stages": STAGES}, indent=1) + "\n")
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
