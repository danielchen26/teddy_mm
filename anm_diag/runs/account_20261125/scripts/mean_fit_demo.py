#!/usr/bin/env python3
"""Section 2.1: on the planted system of anm_diag/toy_account.py, the event `m` moves the state along a direction the
readout does not follow and changes the readout directly. No real cell. results/mean_fit_demo.json.

  python mean_fit_demo.py

For `m` and, beside it, for an event the state does carry (`c1`) and one that does not move the state (`n1`): u with no
record, u with its own record under the rule of run 40 (repair.Full: a readout made to follow the event's mean on one
half, read on the other), and what the account reads at the state.
"""
import json

import numpy as np

from _read import AC, RUN
from anm_diag import cycle, toy_account as T
from repair import cross, u_of


def main():
    res = {}
    for seed in range(5):
        d = T.planted(seed)
        y, un = d["y"], d["unit"]
        for e in ("c1", "m", "n1"):
            m = np.isin(un, ["control", e])
            tr = un[m] == e
            first = cycle.halves(un[m], seed)
            o, _inside, _fits = cross(d["Z"][m], y[m], y[m], tr, un[m], first, [e])
            r = AC.reading(d["Z"][m], y[m], un[m], first, e, draws=100, seed=seed)
            res.setdefault(e, []).append({"u_with_no_record": float(u_of(y[m], o["controls"], tr)[0][0]), "u_with_its_own_record_run_40": float(u_of(y[m], o["rows"], tr)[0][0]),
                                          "account": {k: v["share"] for k, v in r["shares"].items()}, "left_among_cells_whose_state_is_alike": 1 - r["shares"]["carried"]["share"]})
    out = {e: {k: float(np.mean([q[k] for q in v])) for k in ("u_with_no_record", "u_with_its_own_record_run_40", "left_among_cells_whose_state_is_alike")} for e, v in res.items()}
    for e, v in out.items():
        print(e, {k: round(x, 2) for k, x in v.items()})
    (RUN / "results" / "mean_fit_demo.json").write_text(json.dumps({"note": "a planted system; five seeds; means", "planted": {"c1": "through the state; control cells teach it", "m": "moves the state; its change is direct", "n1": "does not move the state; its change is direct"}, "mean": out, "by_seed": res}, indent=1) + "\n")


if __name__ == "__main__":
    main()
