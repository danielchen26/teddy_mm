"""After the gate was read (post hoc): does the picture of place A depend on the ridge weight? Point values, pooled.
The gate chose the weight by cross-validation on training cells (0.1 in every fit). Writes results/sensitivity.json."""
import json

import numpy as np

import gate as G
from gate import CELLS, FEATURES, GRID, TARGETS, load


def main():
    res = {"note": "post hoc; place A, pooled point values of u for each ridge weight (units: the mean variance of a coordinate)", "data_sets": {}}
    for ds, genes in TARGETS.items():
        d = load(ds); guide = np.load(FEATURES / ds / "dense.npz")["guide"].astype(str)
        x, y, A, tg, C = d["x"].astype(np.float64), d["y"], d["A"], d["tg"], d["train_ctl"]
        use = (A == 1) | d["audit_ctl"]
        xc = x[C].mean(0); F = np.hstack([x - xc, np.ones((len(x), 1))]); pen = np.r_[np.ones(x.shape[1]), 0.0] * x[C].var(0).mean()
        M, v = G.moments(F[C], y[C])
        rows = {}
        for lam in GRID:
            Minv = np.linalg.inv(M + lam * np.diag(pen)); th0 = Minv @ v
            resid = float(((F[C] @ th0 - y[C]) ** 2).mean())
            pooled = {k: np.zeros(len(A)) for k in ("controls", "rows", "cells")}
            for k_, g in enumerate(genes):
                own = ((A == 1) & (tg == g)) | (d["ctl_group"] == k_)
                ev = G.events(d, guide, [e for e in genes if e != g], False)
                D = np.array([np.r_[x[T].mean(0) - xc, 0.0] for _n, T in ev]); r = np.array([y[T].mean() - y[C].mean() for _n, T in ev])
                S = np.array([C.sum() * (y[T].var(ddof=1) / T.sum() + y[C].var(ddof=1) / C.sum()) / resid for _n, T in ev])
                th, _c, _g = G.update(th0, Minv, D, r, S)
                tr = C | ((A == 1) & (tg != g))
                M2, v2 = G.moments(F[tr], y[tr]); thc = np.linalg.solve(M2 + lam * np.diag(pen), v2)
                for name, t_ in (("controls", th0), ("rows", th), ("cells", thc)):
                    pooled[name][own] = (F @ t_)[own]
            rows[str(lam)] = {k: G.u_point(y[use], p[use], A[use])[0] for k, p in pooled.items()}
        res["data_sets"][ds] = rows
        print(CELLS[ds]); [print(f"  ridge weight {lam:>7s}: " + "  ".join(f"{k} {v:6.3f}" for k, v in r.items())) for lam, r in rows.items()]
    (G.RES / "sensitivity.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
