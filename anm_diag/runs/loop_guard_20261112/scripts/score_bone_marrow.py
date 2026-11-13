"""The guard: the head's own evaluation on bone-marrow test cells (per protein, the correlation between measured counts and
the predicted mean; mean over proteins), for the head as trained and for the five guarded heads of each data set.
Writes results/bone_marrow.json."""
import json

import numpy as np
from scipy.special import expit

from _paths import ARM, BM_TEST, FEATURES, OUT_FEATURES, RES, TARGETS

LAYERS, ACT = ["mlp__net.0", "mlp__net.3", "mlp__net.6", "dec__mu.net.0", "dec__mu.net.3", "dec__mu.net.6"], [True, True, False, True, True, False]


def score(P, x, y, names):
    h = x
    for name, act in zip(LAYERS, ACT):
        h = h @ P[name + ".weight"].T.astype(np.float64) + P[name + ".bias"]
        h = h * expit(h) if act else h
    pred = np.logaddexp(0.0, h)
    r = {n: float(np.corrcoef(y[:, c], pred[:, c])[0, 1]) for c, n in enumerate(names) if y[:, c].std() > 1e-8}
    return {"mean_pearson": float(np.mean(list(r.values()))), "pdl1_pearson": r.get("CD274"), "n_proteins": len(r)}


d = np.load(BM_TEST)
z = d["z"].astype(np.float64); x = z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6); y = d["adt"].astype(np.float64); names = [str(n) for n in d["adt_names"]]
W = dict(np.load(FEATURES / "head_weights.npz")); assert [str(n) for n in W["adt_names"]] == names
out = {"n_cells": int(len(x)), "as_trained": score(W, x, y, names), "heads": {}}
for ds, genes in TARGETS.items():
    for g in genes:
        out["heads"][f"{ds}|{g}"] = score(dict(np.load(OUT_FEATURES / "heads" / f"{ds}__{g}__{ARM}.npz")), x, y, names)
    v = [out["heads"][f"{ds}|{g}"]["mean_pearson"] for g in genes]
    print(ds, f"mean over proteins {min(v):.3f} to {max(v):.3f}, mean of five heads {np.mean(v):.3f} (as trained {out['as_trained']['mean_pearson']:.3f})")
json.dump(out, open(RES / "bone_marrow.json", "w"), indent=1)
