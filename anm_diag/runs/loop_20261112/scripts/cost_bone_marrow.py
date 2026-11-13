"""What the repairs cost on the head's original task. Reads 5,000 test cells of the bone-marrow CITE-seq data the head was
trained on (TEDDY embedding and measured protein counts; file given by the environment variable ANM_BM_TEST, not in git)
and scores every head the way the head's own evaluation does: per protein, the Pearson correlation between measured counts
and the predicted mean, averaged over proteins."""
import json
import os

import numpy as np
from scipy.special import expit

from _paths import FEATURES, OUT_FEATURES, RES, TARGETS

LAYERS, ACT = ["mlp__net.0", "mlp__net.3", "mlp__net.6", "dec__mu.net.0", "dec__mu.net.3", "dec__mu.net.6"], [True, True, False, True, True, False]
silu = lambda a: a * expit(a)


def predict(P, x):
    h = x
    for name, act in zip(LAYERS, ACT):
        h = h @ P[name + ".weight"].T + P[name + ".bias"]
        h = silu(h) if act else h
    return np.logaddexp(0.0, h)                                      # the decoder's mean, up to the constant size factor


def score(P, x, y, names):
    pred = predict({k: v.astype(np.float64) for k, v in P.items() if k != "adt_names"}, x)
    r = {n: float(np.corrcoef(y[:, j], pred[:, j])[0, 1]) for j, n in enumerate(names) if y[:, j].std() > 1e-8}
    return {"mean_pearson": float(np.mean(list(r.values()))), "pdl1_pearson": r.get("CD274"), "n_proteins": len(r)}


d = np.load(os.environ["ANM_BM_TEST"])
z = d["z"].astype(np.float64); x = z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6); y = d["adt"].astype(np.float64); names = [str(n) for n in d["adt_names"]]
W = dict(np.load(FEATURES / "head_weights.npz")); assert [str(n) for n in W["adt_names"]] == names
out = {"n_cells": int(len(x)), "as_trained": score(W, x, y, names), "heads": {}}
print("as trained", out["as_trained"])
for ds, genes in TARGETS.items():
    for part in ("mlp_only", "decoder_only"):
        for held in genes + ["all"]:
            p = OUT_FEATURES / "heads" / f"{ds}__{held}__{part}.npz"
            if p.exists():
                out["heads"][f"{ds}|{held}|{part}"] = score(dict(np.load(p)), x, y, names)
        v = [s["mean_pearson"] for k, s in out["heads"].items() if k.startswith(ds) and k.endswith(part)]
        if v:
            print(ds, part, f"mean Pearson {min(v):.3f} to {max(v):.3f} over {len(v)} heads (as trained {out['as_trained']['mean_pearson']:.3f})")
json.dump(out, open(RES / "cost_bone_marrow.json", "w"), indent=1)
