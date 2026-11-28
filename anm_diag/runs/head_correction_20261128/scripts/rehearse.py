#!/usr/bin/env python3
"""The rehearsal, in one command, on made-up values alone: whether the judge gives what it should when it is handed a
perfect output and a useless one, and whether the whole path runs and answers C2 both ways.

  python rehearse.py --work DIR [--worlds control_is_enough rows_teach]

-> results/rehearsal.json
"""
import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _hc as K  # noqa: E402
from anm_diag import account as AC  # noqa: E402

half = K.by_path("hc_read", HERE / "read.py").half                       # several runs have a script named read.py


def the_judge_itself(d):
    """Hand `account.reading` an output that is the readout itself, and one that is noise. It should give one and zero."""
    C = K.cells(d / "cells.tsv")
    St, Pk = K.B3.States([d / "states.npz"]), K.B3.Packs([d / "pack.npz"])
    rng = np.random.default_rng(0)
    out = {"perfect": [], "noise": [], "half_of_it": []}
    for c in K.CONDITIONS:
        g = C[c]
        keep = g["drawn"] | (g["control"] & ~g["teach"])
        ids = g["ids"][keep]
        Z = K.B3.V.unit_length(St.rows(ids))
        y = Pk.readouts(K.STUDY, ids, ["PDL1"])[0]["PDL1"]
        unit = np.where(g["control"][keep], "control", g["unit"][keep])
        first = np.array([half(i) for i in ids])
        for u in sorted(set(unit.tolist()) - {"control"})[:4]:
            m = (unit == "control") | (unit == u)
            for name, o in (("perfect", y), ("noise", rng.normal(size=len(y))), ("half_of_it", 0.5 * y + 0.5 * rng.normal(size=len(y)) * y.std())):
                r = AC.reading(Z[m], y[m], unit[m], first[m], u, o=o[m], seed=K.SEED, draws=100)
                if r["snr"] >= AC.MIN_SNR:
                    out[name].append(r["shares"]["as_trained"]["share"])
    return {k: {"n": len(v), "median": float(np.median(v)) if v else None, "range": [float(min(v)), float(max(v))] if v else None} for k, v in out.items()}


def world(work, name):
    d = work / name
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    run = lambda *a: subprocess.run([sys.executable, *map(str, a)], check=True, capture_output=True, text=True)
    run(HERE / "rehearsal.py", "--out", d, "--seed", 1, "--world", name)
    run(HERE / "fit.py", "--cells", d / "cells.tsv", "--units", d / "units.json", "--states", d / "states.npz", "--packs", d / "pack.npz", "--out", d / "rows.npz")
    (d / "results").mkdir()
    run(HERE / "read.py", "--cells", d / "cells.tsv", "--units", d / "units.json", "--rows", d / "rows.npz", "--states", d / "states.npz", "--packs", d / "pack.npz", "--results", d / "results", "--rehearsal", "--draws", 500)
    R = json.loads((d / "results" / "reading.json").read_text())
    planted = json.loads((d / "planted.json").read_text())
    verdict = lambda v: "holds" if v.get("holds") else ("refuted" if v.get("refuted") else "unresolved")
    return {"world": name, "planted": planted, "pairs": len(R["pairs"]), "the_state_carries": float(np.median([p["carried"] for p in R["pairs"]])) if R["pairs"] else None,
            "pooled": {k: (None if v is None else v["share"]) for k, v in R["pooled"].items()}, "H2-H1": R["differences"]["H2-H1"],
            "lines": {k: verdict(v) for k, v in R["lines"].items()}, "ridge_weight": json.loads((d / "rows_report.json").read_text())["lambda"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--worlds", nargs="*", default=["control_is_enough", "rows_teach"])
    a = ap.parse_args()
    res = {"made_up_values_only": True, "worlds": [world(a.work, w) for w in a.worlds]}
    res["the_judge_itself"] = the_judge_itself(a.work / a.worlds[0])
    (K.RUN / "results" / "rehearsal.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(json.dumps({"the_judge_itself": res["the_judge_itself"], "worlds": [{k: w[k] for k in ("world", "pairs", "the_state_carries", "pooled", "lines")} for w in res["worlds"]]}, indent=1, default=float))


if __name__ == "__main__":
    main()
