"""Extra numbers for the target comparison: the group shift at each stage relative to the spread of control cells, and the shift by batch."""
import json, numpy as np
from _paths import FEATURES, RES
F = str(FEATURES) + "/"
def load(sub):
    d = np.load(F + sub + "/dense.npz"); D = {k: d[k] for k in d.files}
    for k in ("Z", "HEAD_LAT", "HEAD_OUT"):
        if "node__" + k not in D: D["node__" + k] = np.load(F + sub + f"/node__{k}.npy")
    return D
out = {}
for sub, targets in (("papalexi2021", ["IFNGR1", "IFNGR2", "JAK2", "STAT1", "IRF1"]), ("frangieh2021", ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"])):
    D = load(sub); A = D["A"]; tg = D["target"].astype(str); C = A == 0; r = {}
    for g in targets:
        T = (A == 1) & (tg == g); r[g] = {}
        for st in ("Z", "HEAD_LAT", "HEAD_OUT"):
            X = D["node__" + st].astype(np.float64); dm = X[T].mean(0) - X[C].mean(0); spread = np.sqrt(X[C].var(0).sum())
            r[g][st] = float(np.linalg.norm(dm) / spread)
    out[sub] = r
    print(sub); [print("  ", g, {k: round(v, 3) for k, v in x.items()}, "kept from Z to head output:", round(x["HEAD_OUT"] / x["Z"], 2)) for g, x in r.items()]
    if "stratum" in D:
        st = D["stratum"].astype(str); Z = D["node__Z"].astype(np.float64); b = {}
        for s in sorted(set(st)):
            m = st == s
            if (m & (A == 1)).sum() < 100: continue
            d = Z[m & (A == 1)].mean(0) - Z[m & C].mean(0); b[s] = float(np.linalg.norm(d))
        d_all = Z[A == 1].mean(0) - Z[C].mean(0); out[sub + "_batch_shift"] = b; out[sub + "_pooled_shift"] = float(np.linalg.norm(d_all)); print("   batch |dZ|:", {k: round(v, 2) for k, v in b.items()}, "pooled", round(float(np.linalg.norm(d_all)), 2))
json.dump(out, open(RES / "why_gap_over_spread.json", "w"), indent=1)
