#!/usr/bin/env python3
"""A check added after the reading: its stages and SPEC files. No readout enters.

  python after_build.py --enc WORK_ENC --typed OUT --after OUT_AFTER --work WORK_AFTER [--results RESULTS]

  V, Z70, TA, EA            as in the reading (for the pairs)
  TAx, EAx                  after_states.py: the typing of TA without the knocked-out genes
  Zs, TAs, TRs, TAxs, EAs   "at one scale": what the round read of Z70, TA, TR, TAx and EA (the state divided by the
                            length of the mean state), every column divided by its spread among the control cells of
                            the stratum. No knockout cell and no protein value enters the scale.

Why the scale. The round's linear reader centres the columns and penalises all of them alike; the tool's reference
reader brings every column to one scale first. A block that is a small part of the mean has small columns, so the
round may not use it where the reference reader does. `Zs` is the control: the same scaling of the mean alone.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import yaml

TYPE_A = ["IFNGR1", "IFNGR2", "JAK2", "STAT1"]
STAGES = ["V", "Z70", "TA", "EA", "TAx", "EAx", "Zs", "TAs", "TRs", "TAxs", "EAs"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--enc", type=Path, required=True)
    ap.add_argument("--typed", type=Path, required=True)
    ap.add_argument("--after", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=Path(__file__).resolve().parent.parent / "results")
    a = ap.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    K, T, X = np.load(a.enc / "thp1_enc.npz"), np.load(a.typed / "typed.npz"), np.load(a.after / "after.npz")
    ids, un = K["ids"].astype(str), K["unit"].astype(str)
    assert np.array_equal(T["ids"].astype(str), ids) and np.array_equal(X["ids"].astype(str), ids), "states of other cells, or in another order"
    z = K["L12"].astype(np.float64)
    nz = np.linalg.norm(z, axis=1, keepdims=True) + 1e-6
    ctl = K["A"].astype(int) == 0

    def scaled(x):
        r = x.astype(np.float64) / nz                             # what the round read
        sd = r[ctl].std(0)
        return (r / np.where(sd > 1e-12, sd, 1.0)).astype(np.float32), float((sd <= 1e-12).mean())

    S = {"V": K["V"].astype(np.float32), "Z70": K["L12"].astype(np.float32), "TA": T["TA"], "EA": T["EA"], "TAx": X["TAx"], "EAx": X["EAx"]}
    flat = {}
    for name, src in (("Zs", z), ("TAs", T["TA"]), ("TRs", T["TR"]), ("TAxs", X["TAx"]), ("EAs", T["EA"])):
        S[name], flat[name] = scaled(src)
    assert list(S) == STAGES
    rest = {k: K[k] for k in ("ids", "A", "unit", "y_PDL1", "frozen_PDL1")}
    np.savez(a.work / "thp1_af.npz", **S, **rest)
    msk = np.isin(un, TYPE_A + ["control"])
    np.savez(a.work / "thp1a_af.npz", **{k: v[msk] for k, v in S.items()}, **{k: v[msk] for k, v in rest.items()})
    for src, arrays in (("T_thp1_IRF1", "thp1_af.npz"), ("C_thp1a_pool", "thp1a_af.npz")):
        s = yaml.safe_load((a.enc / f"{src}.yaml").read_text())
        s["keys"]["stages"] = {k: k for k in STAGES}
        s["reads"] = "Z70"
        s["arrays"] = {"audit": arrays}
        s["readout_classes"] = {"ridge": s["readout_classes"]["ridge"]}
        (a.work / f"{src}.yaml").write_text(yaml.safe_dump(s, sort_keys=False))
    # how large the block of the events' genes is beside the rest, in what the round read: the reason for the check
    r = T["TA"].astype(np.float64) / nz
    d = z.shape[1]
    sd1, sd2 = r[ctl][:, :d].std(0), r[ctl][:, d:].std(0)
    man = json.loads((a.after / "after_manifest.json").read_text())
    man["stages"] = {k: int(v.shape[1]) for k, v in S.items()}
    man["share_of_columns_without_spread_among_control_cells"] = flat
    man["TA_as_the_round_read_it"] = {"spread_of_a_column_among_control_cells_median": {"block_of_the_events_genes": float(np.median(sd1)), "the_rest": float(np.median(sd2))},
                                     "length_of_the_block_mean_over_cells": {"block_of_the_events_genes": float(np.linalg.norm(r[:, :d], axis=1).mean()), "the_rest": float(np.linalg.norm(r[:, d:], axis=1).mean())}}
    a.results.mkdir(parents=True, exist_ok=True)
    (a.results / "after_manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    print(json.dumps({"stages": man["stages"], "TA_as_the_round_read_it": man["TA_as_the_round_read_it"]}))


if __name__ == "__main__":
    main()
