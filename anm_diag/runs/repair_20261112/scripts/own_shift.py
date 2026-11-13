"""Side measure, outside the rule: how far the head's own outputs move for the knockouts of each held-out gene, next to
the measured protein. Both in standard deviations of the audit controls. "As trained" is the head before any repair;
"repaired" is the head of this round fitted with that gene held out (from results/fits_s0.json). "fitted_genes" gives the
same three numbers for the knockouts of the other genes of the data set, which that head was fitted on (mean over genes),
to tell a readout the fit did not learn from one it learned and did not carry over.
Writes results/own_shift.json."""
import json

import numpy as np
from scipy.special import expit

from _paths import FEATURES, OUT_FEATURES, RES, TARGETS, TRAINED

LAYERS, ACT = ["mlp__net.0", "mlp__net.3", "mlp__net.6", "dec__mu.net.0", "dec__mu.net.3", "dec__mu.net.6"], [True, True, False, True, True, False]


def head_out(P, x):
    h = x
    for name, act in zip(LAYERS, ACT):
        h = h @ P[name + ".weight"].T.astype(np.float64) + P[name + ".bias"]
        h = h * expit(h) if act else h
    return np.log(np.maximum(np.logaddexp(0.0, h), 1e-6))


W = dict(np.load(FEATURES / "head_weights.npz")); names = [str(n) for n in W["adt_names"]]
side = json.load(open(RES / "fits_s0.json"))["side"]; out = {}
for ds, genes in TARGETS.items():
    f = np.load(FEATURES / ds / "dense.npz"); Z = f["node__Z"] if "node__Z" in f.files else np.load(FEATURES / ds / "node__Z.npy")
    x = Z.astype(np.float64); x = x / (np.linalg.norm(x, axis=1, keepdims=True) + 1e-6); o0 = head_out(W, x)
    keep = set(json.load(open(RES / f"audit_controls_{ds}_s0.json"))["audit_controls"]); C = np.array([i in keep for i in f["ids"].astype(str)])
    A, tg = f["A"], f["target"].astype(str); out[ds] = {}
    for g in genes:
        T = (A == 1) & (tg == g); out[ds][g] = {}
        for rid, p in TRAINED[ds].items():
            y, j = f["readout__" + rid], names.index(p)
            out[ds][g][rid] = {"measured": float((y[T].mean() - y[C].mean()) / y[C].std()), "as_trained": float((o0[T, j].mean() - o0[C, j].mean()) / o0[C, j].std()),
                               "repaired": side[ds][g][rid]["own_output_shift_sd"]}
            assert abs(out[ds][g][rid]["measured"] - side[ds][g][rid]["measured_shift_sd"]) < 1e-6
        og = head_out(dict(np.load(OUT_FEATURES / "heads" / f"s0__{g}.npz")), x)      # the head fitted with g held out
        for rid, p in TRAINED[ds].items():
            y, j = f["readout__" + rid], names.index(p); rows = []
            assert abs((og[T, j].mean() - og[C, j].mean()) / og[C, j].std() - out[ds][g][rid]["repaired"]) < 1e-3
            for h in genes:
                if h != g:
                    Th = (A == 1) & (tg == h)
                    rows.append([(y[Th].mean() - y[C].mean()) / y[C].std(), (o0[Th, j].mean() - o0[C, j].mean()) / o0[C, j].std(), (og[Th, j].mean() - og[C, j].mean()) / og[C, j].std()])
            out[ds][g][rid]["fitted_genes"] = dict(zip(("measured", "as_trained", "repaired"), (float(v) for v in np.mean(rows, 0))))
        print(ds, g, " | ".join(f"{rid} measured {v['measured']:+.2f}, as trained {v['as_trained']:+.2f}, repaired {v['repaired']:+.2f} (fitted genes: {v['fitted_genes']['measured']:+.2f}, {v['fitted_genes']['as_trained']:+.2f}, {v['fitted_genes']['repaired']:+.2f})" for rid, v in out[ds][g].items()))
json.dump(out, open(RES / "own_shift.json", "w"), indent=1)
