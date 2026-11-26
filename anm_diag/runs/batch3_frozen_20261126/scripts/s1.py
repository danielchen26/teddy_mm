#!/usr/bin/env python3
"""S1: the system as it forecast batch 2, run with the package as it was then. Called by forecast.py.

  python s1.py --package OLD --bundle BUNDLE.npz --out S1.json

OLD: a folder that holds the package `anm_diag` of commit 903532a (git archive 903532a anm_diag/anm_diag). The steps are
those of runs/batch2_forecast_20261125/scripts/forecast.py: one law of the library under interferon; per stratum a law of
its recorded rows, every recorded row of one type, with the spread of a direct part with no edge from its own rows where
it has at least eight recorded pairs with no edge, else from the library; one size per class of edge, pooled over the
recorded pairs of the class in the library and the strata of the same study. A stratum with no recorded row holds a
placeholder event (a copy of its control cells) and takes the library's spread.
"""
import argparse
import json
import sys

import numpy as np

MIN_ROWS = 8


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--package", required=True)
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sys.path.insert(0, a.package)
    from anm_diag import account as AC
    assert not hasattr(AC, "MEDIAN_Z2"), "this is not the package of batch 2"
    B = np.load(a.bundle)
    meta = json.loads(str(B["meta"]))
    E = {(e, r): "own" for e, r in meta["edges"]}
    lib = AC.Law(B["lib_x"], {r: B["lib_Y"][:, j] for j, r in enumerate(meta["lib_readouts"])}, B["lib_unit"], edges=E).fit()
    laws = {}
    for m in meta["laws"]:
        k = m["key"]
        ed = E if m["study"] == "frangieh2021" else {}
        law = AC.Law(B[k + "_x"], {r: B[k + "_Y"][:, j] for j, r in enumerate(m["readouts"])}, B[k + "_unit"], edges=ed)
        free = 0 if m["placeholder"] else sum(1 for e in law.events for r in m["readouts"] if (e, r) not in ed)
        law.fit(sigma0=None if free >= MIN_ROWS else lib.sigma0)
        laws[k] = (law, ed, m)
    pooled = [v for _e, _r, v in lib.class_values().get("own", [])]
    for k, (law, ed, m) in laws.items():
        if m["study"] == "frangieh2021":
            pooled += [v for _e, _r, v in law.class_values().get("own", [])]
    out = {"package": "903532a", "library": {"events": len(lib.events), "spread_of_a_direct_part_without_an_edge": lib.sigma0},
           "class_own_pooled": {"values": len(pooled), "mean": float(np.mean(pooled)), "spread_of_a_new_pair": float(np.std(pooled, ddof=1) * np.sqrt(1 + 1 / len(pooled)))}, "laws": {}, "rows": []}
    for k, (law, ed, m) in laws.items():
        s = law.summary()
        out["laws"][k] = {"recorded_events": law.events, "spread_of_a_direct_part_without_an_edge": law.sigma0, "from": "the library" if law.fixed else "own rows",
                          "widening_of_the_control_prior": s["widening_of_the_control_prior"]}
        for r in law.new(B[k + "_xn"], B[k + "_un"], edges=ed, classes={"own": pooled}):
            out["rows"].append({"law": k, "event": r["event"], "readout": r["readout"], "read": r["read"], "from_the_recorded_events": r["from_the_recorded_events"], "by_its_edge": r["by_its_edge"],
                                "forecast": r["forecast"], "sd": r["sd"], "interval": r["interval"], "no_forecast": r["no_forecast"]})
        print("S1", k, len(law.events), "recorded events; spread", round(law.sigma0, 3), out["laws"][k]["from"], flush=True)
    open(a.out, "w").write(json.dumps(out, default=float) + "\n")


if __name__ == "__main__":
    main()
