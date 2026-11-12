"""Two questions the target comparison left open. Exploratory, on the exported features and the head's weights.

1. Why does the head's last step shrink IRF1 alone? The decoder maps 512 latent dimensions to 134 proteins, so around any
   point at least 378 latent directions cannot change the output. For every target: the share of its latent shift that lies in
   the 134 directions the decoder can see, the share on the decoder's ten strongest directions, the gain of the linearised
   decoder along the shift, and the decoder's reaction to that shift against random latent shifts of equal length.
2. Why does THP-1 batch 1 have a smaller mean shift without a larger loss? If fewer of its cells are perturbed while a
   perturbed cell moves as far as in the other batches, the mean shrinks and what the state says about each cell does not.
   Cells are projected on the shared direction of the four upstream knockouts; a treated cell counts as shifted when it lies
   beyond the 95th percentile of its batch's control cells. The slope of PD-L1 on that projection says how far PD-L1 moves
   per unit of state shift in each batch."""
import json
import numpy as np
from scipy.special import expit

from _paths import FEATURES, RES

W = dict(np.load(FEATURES / "head_weights.npz"))
silu = lambda a: a * expit(a)


def mlp(p, x):
    h = silu(x @ W[p + "net.0.weight"].T + W[p + "net.0.bias"])
    h = silu(h @ W[p + "net.3.weight"].T + W[p + "net.3.bias"])
    return h @ W[p + "net.6.weight"].T + W[p + "net.6.bias"]


def latent(z):
    z = np.asarray(z, np.float64)
    return mlp("mlp__", z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6))


def decode(lat):
    return np.log(np.maximum(np.logaddexp(0.0, mlp("dec__mu.", lat)), 1e-6))


def load(sub):
    d = np.load(FEATURES / sub / "dense.npz")
    D = {k: d[k] for k in d.files}
    if "node__Z" not in D:
        D["node__Z"] = np.load(FEATURES / sub / "node__Z.npy")
    return D


rng = np.random.default_rng(11)
out = {"decoder": {}, "batches": {}}
for sub, targets in (("papalexi2021", ["IFNGR1", "IFNGR2", "JAK2", "STAT1", "IRF1"]), ("frangieh2021", ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"])):
    D = load(sub)
    A, tg, Z = D["A"], D["target"].astype(str), D["node__Z"].astype(np.float64)
    lat = latent(Z)
    C = np.where(A == 0)[0]
    Ls = lat[rng.choice(C, min(600, len(C)), replace=False)]
    o0 = decode(Ls)
    J = np.zeros((o0.shape[1], lat.shape[1]))
    h = 0.05 * float(lat[C].std())
    for j in range(lat.shape[1]):                       # Jacobian of the decoder, averaged over control cells
        e = np.zeros(lat.shape[1]); e[j] = h
        J[:, j] = (decode(Ls + e) - decode(Ls - e)).mean(0) / (2 * h)
    U, s, Vt = np.linalg.svd(J, full_matrices=False)
    res = {"singular_values": [float(s[0]), float(s[len(s) // 2]), float(s[-1])], "rank": int((s > 1e-8 * s[0]).sum()), "chance_share": Vt.shape[0] / Vt.shape[1],
           "effective_rank": float((s ** 2).sum() ** 2 / (s ** 4).sum()), "targets": {}}
    for g in targets:
        T = np.where((A == 1) & (tg == g))[0]
        d = lat[T].mean(0) - lat[C].mean(0); n = float(np.linalg.norm(d))
        vis = float(((Vt @ d) ** 2).sum() / n ** 2)
        top10 = float(((Vt[:10] @ d) ** 2).sum() / n ** 2)
        wvis = float((((Vt @ d) * s) ** 2).sum() / n ** 2) ** 0.5       # gain of the linearised decoder along the shift
        move = float(np.linalg.norm(decode(Ls + d) - o0, axis=1).mean())
        rnd = []
        for _ in range(200):
            v = rng.standard_normal(lat.shape[1]); v *= n / np.linalg.norm(v)
            rnd.append(float(np.linalg.norm(decode(Ls + v) - o0, axis=1).mean()))
        res["targets"][g] = dict(latent_shift=n, visible_share=vis, top10_share=top10, linear_gain=wvis, output_move=move, output_move_vs_random=move / float(np.median(rnd)))
    out["decoder"][sub] = res
    print(sub, "rank", res["rank"], "effective rank", round(res["effective_rank"], 1), "chance share", round(res["chance_share"], 3))
    for g, t in res["targets"].items():
        print(f"   {g:7s} latent shift {t['latent_shift']:.2f} | visible share {t['visible_share']:.2f} | on ten strongest {t['top10_share']:.2f} | gain {t['linear_gain']:.3f} | output move vs random {t['output_move_vs_random']:.2f}")

D = load("papalexi2021")
A, tg, st, Z = D["A"], D["target"].astype(str), D["stratum"].astype(str), D["node__Z"].astype(np.float64)
y = D["readout__PDL1"]
up = ["IFNGR1", "IFNGR2", "JAK2", "STAT1"]
isup = (A == 1) & np.isin(tg, up)
d = Z[isup].mean(0) - Z[A == 0].mean(0); d /= np.linalg.norm(d)
for b in ("rep1-tx", "rep3-tx", "rep4-tx"):
    m = st == b; c = m & (A == 0); t = m & isup
    p = (Z - Z[c].mean(0)) @ d
    thr = float(np.quantile(p[c], 0.95)); sh = t & (p > thr); ns = t & ~(p > thr)
    out["batches"][b] = dict(n_treated=int(t.sum()), n_control=int(c.sum()), mean_shift=float(p[t].mean()), share_shifted=float(sh.sum() / t.sum()),
                             shift_of_shifted=float(p[sh].mean()), pdl1_change_shifted=float(y[sh].mean() - y[c].mean()), pdl1_change_not_shifted=float(y[ns].mean() - y[c].mean()),
                             pdl1_change_all=float(y[t].mean() - y[c].mean()), control_pdl1_mean=float(y[c].mean()), control_pdl1_sd=float(y[c].std()),
                             treated_pdl1_p5=float(np.quantile(y[t], 0.05)), slope_pdl1_on_shift=float(np.polyfit(p[c | t], y[c | t], 1)[0]),
                             corr_pdl1_shift=float(np.corrcoef(p[c | t], y[c | t])[0, 1]),
                             share_shifted_by_target={g: float((p[m & (A == 1) & (tg == g)] > thr).mean()) for g in up + ["IRF1"]})
    x = out["batches"][b]
    print(f"   {b}: treated {x['n_treated']}, mean shift {x['mean_shift']:.2f}, shifted {100 * x['share_shifted']:.0f}% (5% expected by chance), their shift {x['shift_of_shifted']:.2f}, PD-L1 change: shifted {x['pdl1_change_shifted']:+.2f}, not shifted {x['pdl1_change_not_shifted']:+.2f}, all {x['pdl1_change_all']:+.2f}; PD-L1 per unit of shift {x['slope_pdl1_on_shift']:+.3f} (r {x['corr_pdl1_shift']:+.2f}); shifted by target { {g: round(v, 2) for g, v in x['share_shifted_by_target'].items()} }")
json.dump(out, open(RES / "why_last_step.json", "w"), indent=1)
