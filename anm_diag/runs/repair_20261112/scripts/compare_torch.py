"""Compare heads fitted in plain NumPy with the ones scripts/07_repair_head.py (teddy_mm's own model code, PyTorch) makes
from the same inputs. Writes results/torch_check.json.

  ANM_TORCH_HEAD     the head of this round fitted on all genes (heads/s0__all.npz here), as made by the script
  ANM_TORCH_ROUND2   optional: a folder with round2_torch_<data set>_IFNGR1_weights.npz, the script run with one data set,
                     PD-L1 only and IFNGR1 held out, which is the repair of runs/loop_guard_20261112
  ANM_GUARD_FEATURES the folder of that run's heads (heads/<data set>__IFNGR1__mlp_guarded.npz)

All weights are npz files with the names of head_weights.npz."""
import json
import os

import numpy as np
from scipy.special import expit

from _paths import BM_TEST, FEATURES, OUT_FEATURES, RES, TARGETS

LAYERS, ACT = ["mlp__net.0", "mlp__net.3", "mlp__net.6", "dec__mu.net.0", "dec__mu.net.3", "dec__mu.net.6"], [True, True, False, True, True, False]


def out(P, x):
    h = x
    for name, act in zip(LAYERS, ACT):
        h = h @ P[name + ".weight"].T.astype(np.float64) + P[name + ".bias"]
        h = h * expit(h) if act else h
    sp = np.logaddexp(0.0, h)
    return np.log(np.maximum(sp, 1e-6)), sp


def unit(z):
    z = z.astype(np.float64)
    return z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)


W0 = dict(np.load(FEATURES / "head_weights.npz")); A = dict(np.load(OUT_FEATURES / "heads" / "s0__all.npz")); B = dict(np.load(os.environ["ANM_TORCH_HEAD"]))
res = {"decoder_identical_to_as_trained": {"numpy": all(np.array_equal(A[l + s], W0[l + s]) for l in LAYERS[3:] for s in (".weight", ".bias")),
                                           "torch": all(np.array_equal(B[l + s], W0[l + s]) for l in LAYERS[3:] for s in (".weight", ".bias"))},
       "mlp_weights": {"largest_change_from_as_trained": float(max(np.abs(A[l + ".weight"] - W0[l + ".weight"]).max() for l in LAYERS[:3])),
                       "largest_difference_numpy_torch": float(max(np.abs(A[l + s] - B[l + s]).max() for l in LAYERS[:3] for s in (".weight", ".bias")))}}
d = np.load(BM_TEST); x, y = unit(d["z"]), d["adt"].astype(np.float64)
score = lambda sp: float(np.mean([np.corrcoef(y[:, c], sp[:, c])[0, 1] for c in range(y.shape[1]) if y[:, c].std() > 1e-8]))
oa, spa = out(A, x); ob, spb = out(B, x); o0, sp0 = out(W0, x)
res["bone_marrow_test"] = {"n_cells": int(len(x)), "score_as_trained": score(sp0), "score_numpy": score(spa), "score_torch": score(spb),
                           "outputs_rms_change_from_as_trained": float(np.sqrt(((oa - o0) ** 2).mean())), "outputs_rms_difference_numpy_torch": float(np.sqrt(((oa - ob) ** 2).mean()))}
res["perturbation_cells"] = {}
for ds in TARGETS:
    f = np.load(FEATURES / ds / "dense.npz"); Z = f["node__Z"] if "node__Z" in f.files else np.load(FEATURES / ds / "node__Z.npy")
    xa = unit(Z); pa, _ = out(A, xa); pb, _ = out(B, xa); p0, _ = out(W0, xa)
    res["perturbation_cells"][ds] = {"outputs_rms_change_from_as_trained": float(np.sqrt(((pa - p0) ** 2).mean())), "outputs_rms_difference_numpy_torch": float(np.sqrt(((pa - pb) ** 2).mean()))}
if os.environ.get("ANM_TORCH_ROUND2"):
    from pathlib import Path
    res["round_2_recipe"] = {"what": "scripts/07_repair_head.py with one data set, PDL1=CD274 only, --hold-out IFNGR1, against heads/<data set>__IFNGR1__mlp_guarded.npz of runs/loop_guard_20261112"}
    for ds in TARGETS:
        A2 = dict(np.load(Path(os.environ["ANM_GUARD_FEATURES"]) / "heads" / f"{ds}__IFNGR1__mlp_guarded.npz")); B2 = dict(np.load(Path(os.environ["ANM_TORCH_ROUND2"]) / f"round2_torch_{ds}_IFNGR1_weights.npz"))
        f = np.load(FEATURES / ds / "dense.npz"); Z = f["node__Z"] if "node__Z" in f.files else np.load(FEATURES / ds / "node__Z.npy")
        held = (f["A"] == 1) & (f["target"].astype(str) == "IFNGR1"); xa = unit(Z)
        pa, _ = out(A2, xa); pb, _ = out(B2, xa); p0, _ = out(W0, xa); ba, sa = out(A2, x); bb, sb = out(B2, x)
        rms = lambda u, v, m=slice(None): float(np.sqrt(((u[m] - v[m]) ** 2).mean()))
        res["round_2_recipe"][ds] = {"decoder_identical_to_as_trained": all(np.array_equal(B2[l + s], W0[l + s]) for l in LAYERS[3:] for s in (".weight", ".bias")),
                                     "perturbation_cells": {"outputs_rms_change_from_as_trained": rms(pa, p0), "outputs_rms_difference_numpy_torch": rms(pa, pb)},
                                     "held_out_knockout_cells": {"n": int(held.sum()), "outputs_rms_change_from_as_trained": rms(pa, p0, held), "outputs_rms_difference_numpy_torch": rms(pa, pb, held)},
                                     "bone_marrow_test": {"score_numpy": score(sa), "score_torch": score(sb), "outputs_rms_difference_numpy_torch": rms(ba, bb)}}
json.dump(res, open(RES / "torch_check.json", "w"), indent=1)
print(json.dumps(res, indent=1))
