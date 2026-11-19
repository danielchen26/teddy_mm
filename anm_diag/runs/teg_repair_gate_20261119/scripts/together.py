"""After the gate was read (post hoc, development): both cell types repaired at once, in one head.

One gate with three ways (bone marrow, melanoma, THP-1), fitted on embeddings alone (the training control cells of each
data set and the bone-marrow training cells). The corrected log output of a cell is
    p_bm(x) * head(x) + p_mel(x) * readout_mel(x) + p_thp(x) * readout_thp(x)
with each readout the `rows` readout of place A, fitted for its own cell type on the head's scale. Read as in posthoc.py
(pooled, log scale), with the readout of the other cell type present in the same output.

  python together.py     (same environment variables as gate.py; writes results/together.json)
"""
import json

import numpy as np
from sklearn.linear_model import LogisticRegression

import gate as G
import posthoc as PH
from gate import BM_FIT, BM_TEST, CELLS, FEATURES, TARGETS, forward_np, load, unit, weights


def main():
    P0, names, j = weights()
    b = np.load(BM_FIT); t_ = np.load(BM_TEST)
    xb, xt, yt = unit(b["z_train"]).astype(np.float64), unit(t_["z"]).astype(np.float64), t_["adt"].astype(np.float64)
    out_t = forward_np(P0, xt)[1]
    data = {}
    for ds, genes in TARGETS.items():
        d = load(ds); guide = np.load(FEATURES / ds / "dense.npz")["guide"].astype(str)
        x, y, C = d["x"].astype(np.float64), d["y"], d["train_ctl"]
        o = forward_np(P0, d["x"])[1][:, j]
        xc = x[C].mean(0); scale = o[C].std() / y[C].std()
        target = o[C].mean() + scale * (y - y[C].mean())
        F = np.hstack([x - xc, np.ones((len(x), 1))])
        th = {k: {g: PH.fit_kinds(d, guide, g, genes, F, target, C, None, scale)[k] for g in genes} for k in ("rows",)}
        data[ds] = dict(d=d, x=x, o=o, xc=xc, th=th["rows"], genes=genes)
    order = list(TARGETS)                                                  # classes: 0 bone marrow, 1 and 2 the two cell types
    X = np.vstack([xb] + [data[ds]["x"][data[ds]["d"]["train_ctl"]] for ds in order])
    lab = np.r_[np.zeros(len(xb)), *[np.full(int(data[ds]["d"]["train_ctl"].sum()), i + 1) for i, ds in enumerate(order)]]
    clf = LogisticRegression(C=1.0, max_iter=3000).fit(X, lab)
    feat = lambda z, ds: np.hstack([z - data[ds]["xc"], np.ones((len(z), 1))])
    res = {"note": "post hoc, development on cells already read: both cell types repaired in one output", "data_sets": {}}
    def combined(z, head_out, pick):
        """pick: for each data set, which of its five fits to use for these cells."""
        p = clf.predict_proba(z)
        return p[:, 0] * head_out + sum(p[:, i + 1] * (feat(z, ds) @ data[ds]["th"][pick[ds]]) for i, ds in enumerate(order)), p
    for ds in order:
        D = data[ds]; other = [o_ for o_ in order if o_ != ds][0]
        logs = {g: combined(D["x"], D["o"], {ds: g, other: data[other]["genes"][0]})[0] for g in D["genes"]}
        r = PH.pooled(D["d"], D["genes"], logs)
        p = clf.predict_proba(D["x"]); A = D["d"]["A"]
        res["data_sets"][ds] = {"reading": r, "gate_on_knockout_cells": {("bone_marrow" if i == 0 else order[i - 1]): float(p[A == 1, i].mean()) for i in range(3)},
                                "gate_on_audit_controls": {("bone_marrow" if i == 0 else order[i - 1]): float(p[D["d"]["audit_ctl"], i].mean()) for i in range(3)}}
        print(f"{CELLS[ds]}: log {r['log']['u']:.3f} [{r['log']['ci'][0]:.3f}, {r['log']['ci'][1]:.3f}] {r['log']['reading']} {'ok' if r['log']['valid'] else 'NO'} | gate on knockout cells",
              {k: round(v, 3) for k, v in res["data_sets"][ds]["gate_on_knockout_cells"].items()})
    new_j, p = combined(xt, out_t[:, j], {ds: data[ds]["genes"][0] for ds in order})
    res["bone_marrow"] = dict(PH.bm(out_t, j, yt, new_j), gate={("bone_marrow" if i == 0 else order[i - 1]): float(p[:, i].mean()) for i in range(3)})
    print("bone marrow test cells:", {k: (round(v, 4) if not isinstance(v, dict) else {a: round(b_, 4) for a, b_ in v.items()}) for k, v in res["bone_marrow"].items()})
    (G.RES / "together.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
