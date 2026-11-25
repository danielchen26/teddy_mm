#!/usr/bin/env python3
"""The typed states as stages, and the SPEC files of the reading. No readout enters.

  python build_ty.py --enc WORK_ENC --typed OUT --work WORK [--results RESULTS]

WORK_ENC: the work folder of runs/encoder_stages_20261124 (thp1_enc.npz and the two SPEC files). OUT: the folder of
typed.py. Stages: V (the input view), Z70 (the state a head reads today, as recorded), TK8, TA, TR, G (typed.py),
ZG (Z70 beside G), EA (typed.py), PERM. The events, the decision and the keys of the SPEC files are those of the runs before; the
readers are the linear class only (the class the tool's rule has chosen at every stage so far), because three of
the stages are wide.
Writes WORK/thp1_ty.npz, thp1a_ty.npz, the two SPEC files and RESULTS/typed_manifest.json (a copy of typed.py's).
"""
import argparse
import json
from pathlib import Path

import numpy as np
import yaml

TYPE_A = ["IFNGR1", "IFNGR2", "JAK2", "STAT1"]
STAGES = ["V", "Z70", "TK8", "TA", "TR", "G", "ZG", "EA", "PERM"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--enc", type=Path, required=True)
    ap.add_argument("--typed", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=Path(__file__).resolve().parent.parent / "results")
    a = ap.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    K, T = np.load(a.enc / "thp1_enc.npz"), np.load(a.typed / "typed.npz")
    ids, un = K["ids"].astype(str), K["unit"].astype(str)
    assert np.array_equal(T["ids"].astype(str), ids), "the typed states are of other cells, or in another order"
    assert float(np.abs(T["Zc"] - K["L12"]).max()) <= 2e-3, "the mean computed with the typed states is not the recorded state"
    z = K["L12"].astype(np.float32)
    S = {"V": K["V"].astype(np.float32), "Z70": z, "TK8": T["TK8"], "TA": T["TA"], "TR": T["TR"], "G": T["G"], "ZG": np.concatenate([z, T["G"]], 1), "EA": T["EA"], "PERM": K["PERM"].astype(np.float32)}
    assert list(S) == STAGES
    rest = {k: K[k] for k in ("ids", "A", "unit", "y_PDL1", "frozen_PDL1")}
    np.savez(a.work / "thp1_ty.npz", **S, **rest)
    msk = np.isin(un, TYPE_A + ["control"])
    np.savez(a.work / "thp1a_ty.npz", **{k: v[msk] for k, v in S.items()}, **{k: v[msk] for k, v in rest.items()})
    for src, arrays in (("T_thp1_IRF1", "thp1_ty.npz"), ("C_thp1a_pool", "thp1a_ty.npz")):
        s = yaml.safe_load((a.enc / f"{src}.yaml").read_text())  # the same events, decision and keys as the runs before; the stages and the classes change
        s["keys"]["stages"] = {k: k for k in STAGES}
        s["reads"] = "Z70"
        s["arrays"] = {"audit": arrays}
        s["readout_classes"] = {"ridge": s["readout_classes"]["ridge"]}
        (a.work / f"{src}.yaml").write_text(yaml.safe_dump(s, sort_keys=False))
    a.results.mkdir(parents=True, exist_ok=True)
    man = json.loads((a.typed / "typed_manifest.json").read_text())
    man["stages"] = {k: int(v.shape[1]) for k, v in S.items()}
    (a.results / "typed_manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    print(json.dumps(man["stages"]))


if __name__ == "__main__":
    main()
